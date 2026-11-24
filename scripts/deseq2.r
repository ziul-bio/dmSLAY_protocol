cat(
  "\n****************************************\n",
  "        Starting DESeq2 analysis...      \n",
  "****************************************\n",
  sep = ""
)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, aliases = character(), default = NULL) {
  names <- c(name, aliases)
  for (argument in names) {
    prefix <- paste0("--", argument, "=")
    match <- grepl(paste0("^", prefix), args)
    if (any(match)) {
      return(sub(prefix, "", args[which(match)[1]], fixed = TRUE))
    }
    index <- which(args == paste0("--", argument))
    if (length(index) == 1 && index < length(args)) {
      return(args[index + 1])
    }
  }
  default
}

count_file <- get_arg("count-file", c("counts", "count"), "results/counts/peptide_read_counts.txt")
coldata_file <- get_arg("coldata-file", c("coldata", "samples"))
groups_arg <- get_arg("groups")
reference <- get_arg("reference", c("ref"))
treatment <- get_arg("treatment", c("contrast"))
output_dir <- get_arg("output-dir", c("output"), "results/deseq2")
alpha <- as.numeric(get_arg("alpha", default = "0.05"))
lfc_threshold <- as.numeric(get_arg("lfc-threshold", default = "1"))
minimum_count <- as.integer(get_arg("minimum-count", default = "30"))

if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
  stop("--alpha must be a number greater than 0 and less than 1")
}
if (is.na(lfc_threshold) || lfc_threshold < 0) {
  stop("--lfc-threshold must be zero or greater")
}
if (is.na(minimum_count) || minimum_count <= 0) {
  stop("--minimum-count must be greater than 0")
}
if (!file.exists(count_file)) {
  stop("Count table not found: ", count_file)
}

count_matrix <- read.delim(count_file, check.names = FALSE)
metadata_columns <- intersect(
  c("ID", "peptide", "row_sum", "peptide_length"),
  names(count_matrix)
)
if (!"ID" %in% names(count_matrix) || !"peptide" %in% names(count_matrix)) {
  stop("Count table must contain ID and peptide columns")
}

sample_columns <- setdiff(names(count_matrix), metadata_columns)
if (length(sample_columns) == 0) {
  stop("Count table contains no sample columns")
}

counts <- as.matrix(count_matrix[, sample_columns, drop = FALSE])
storage.mode(counts) <- "integer"
if (any(is.na(counts)) || any(counts < 0)) {
  stop("Sample count columns must contain non-negative integers")
}
rownames(counts) <- count_matrix$ID



# Parse comma- or semicolon-separated SAMPLE=GROUP (or SAMPLE:GROUP) entries.
# Returns a named vector aligned to sample_columns; missing samples remain NA.
# Stops if an entry is malformed or names a sample absent from the count table.
parse_groups <- function(value) {
  entries <- trimws(strsplit(value, "[,;]", perl = TRUE)[[1]])
  entries <- entries[nzchar(entries)]
  groups <- rep(NA_character_, length(sample_columns))
  names(groups) <- sample_columns

  for (entry in entries) {
    pair <- regexec("^([^:=]+)\\s*[:=]\\s*(.+)$", entry)
    match <- regmatches(entry, pair)[[1]]
    if (length(match) != 3) {
      stop("Invalid --groups entry '", entry, "'. Use SAMPLE=GROUP.")
    }
    sample <- trimws(match[2])
    group <- trimws(match[3])
    if (!sample %in% sample_columns) {
      stop("Unknown sample in --groups: ", sample)
    }
    groups[sample] <- group
  }
  groups
}


# if --groups is provided, parse it; otherwise, read the coldata file.
if (!is.null(groups_arg)) {
  coldata <- data.frame(
    row.names = sample_columns,
    group = parse_groups(groups_arg),
    check.names = FALSE
  )
} else if (!is.null(coldata_file)) {
  if (!file.exists(coldata_file)) {
    stop("Coldata file not found: ", coldata_file)
  }
  separator <- if (grepl("\\.(csv|txt)$", coldata_file, ignore.case = TRUE)) {
    if (grepl("\\.csv$", coldata_file, ignore.case = TRUE)) "," else "\t"
  } else {
    "\t"
  }

  # Building the coldate data frame from the supplied file. We try to guess which columns are the sample and group columns.
  supplied_coldata <- read.delim(coldata_file, sep = separator, check.names = FALSE)

  sample_column <- intersect(
    c("sample", "sample_id", "sample_name", "Sample", "SampleID", "name"),
    names(supplied_coldata)
  )[1]
  group_column <- intersect(
    c("group", "condition", "treatment", "phenotype"),
    names(supplied_coldata)
  )[1]
  if (is.na(sample_column) || is.na(group_column)) {
    stop("Coldata must contain a sample column and a group/condition column")
  }

  duplicate_samples <- supplied_coldata[[sample_column]][
    duplicated(supplied_coldata[[sample_column]])
  ]
  if (length(duplicate_samples) > 0) {
    stop("Duplicate sample in coldata: ", duplicate_samples[1])
  }

# building the coldata data frame with row names as sample names and a single group column. 
# We check for missing or extra samples compared to the count table.
  coldata <- data.frame(
    row.names = as.character(supplied_coldata[[sample_column]]),
    group = as.character(supplied_coldata[[group_column]]),
    check.names = FALSE
  )
  missing_samples <- setdiff(sample_columns, rownames(coldata))
  extra_samples <- setdiff(rownames(coldata), sample_columns)
  if (length(missing_samples) > 0) {
    stop("Samples missing from coldata: ", paste(missing_samples, collapse = ", "))
  }
  if (length(extra_samples) > 0) {
    stop("Coldata contains samples absent from the count table: ",
         paste(extra_samples, collapse = ", "))
  }
  coldata <- coldata[sample_columns, , drop = FALSE]
} else {
  stop("No sample groups provided. Use --coldata-file or --groups.")
}

if (anyNA(coldata$group) || any(trimws(coldata$group) == "")) {
  stop("Every sample must have a group")
}
coldata$group <- trimws(coldata$group)

groups <- unique(coldata$group)
if (length(groups) < 2) {
  stop("At least two groups are required for DESeq2")
}
if (is.null(reference) || !nzchar(reference)) {
  reference <- groups[1]
}
if (!reference %in% groups) {
  stop("Reference group '", reference, "' was not found in the sample groups")
}
coldata$group <- factor(coldata$group, levels = c(reference, setdiff(groups, reference)))


######################################################
#      Create DESeq2 dataset and run the analysis
######################################################
counts <- counts[, rownames(coldata), drop = FALSE]
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~group
)
keep <- rowSums(counts(dds)) >= minimum_count
dds <- dds[keep, , drop = FALSE]
if (nrow(dds) == 0) {
  stop("No peptides remain after the minimum-count filter")
}

# Main DESeq2 analysis: relevel the group factor to set the reference group, then run DESeq.
print(paste("Running DESeq2 with reference group:", reference))
dds$group <- droplevels(relevel(dds$group, ref = reference))
dds <- DESeq(dds)



# By default, analyze every non-reference group; otherwise parse the requested groups.
if (is.null(treatment) || !nzchar(treatment)) {
  treatments <- setdiff(levels(coldata$group), reference)
} else {
  treatments <- trimws(strsplit(treatment, "[,;]", perl = TRUE)[[1]])
  treatments <- treatments[nzchar(treatments)]
}



# Create the output directory and choose groups to compare against the reference.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through each treatment group and perform DESeq2 analysis against the reference group.
for (treatment_group in treatments) {
  if (!treatment_group %in% levels(coldata$group)) {
    stop("Treatment group '", treatment_group, "' was not found")
  }
  if (treatment_group == reference) {
    next
  }
  # Perform DESeq2 results extraction for the current treatment group compared to the reference group.
  cat("\nAnalyzing ", treatment_group, " vs ", reference, "\n", sep = "")
  result <- results(dds, alpha = alpha,
                    contrast = c("group", treatment_group, reference))
  
  result$DE <- ifelse(
    !is.na(result$padj) & result$padj < alpha &
      abs(result$log2FoldChange) >= lfc_threshold,
    ifelse(result$log2FoldChange > 0, "Upregulated", "Downregulated"),
    "Not significant"
  )
  result$DE <- factor(
    result$DE,
    levels = c("Upregulated", "Downregulated", "Not significant")
  )

  result_table <- as.data.frame(result)
  result_table$ID <- rownames(result_table)
  result_table <- merge(
    count_matrix,
    result_table,
    by = "ID",
    all.x = FALSE,
    sort = FALSE
  )

  # Save the result table to a TSV file in the output directory, with a filename that includes the treatment and reference group names.
  comparison <- paste0(treatment_group, "_vs_", reference)
  safe_comparison <- gsub("[^A-Za-z0-9._-]", "_", comparison)
  write.table(
    result_table,
    file = file.path(output_dir, paste0("deseq2_", safe_comparison, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # this will list only the top DE peptides with the negative log2FoldChange (downregulated) and padj < alpha, sorted by padj.
  top_de <- result_table[
    !is.na(result_table$padj) &
      result_table$padj < alpha &
      result_table$log2FoldChange < -(lfc_threshold),
  ]
  top_de <- top_de[order(top_de$padj, na.last = TRUE), , drop = FALSE]

  up <- sum(result_table$DE == "Upregulated")
  down <- sum(result_table$DE == "Downregulated")

  x_range <- range(result_table$log2FoldChange, finite = TRUE)
  y_range <- range(-log10(result_table$pvalue), finite = TRUE)
  annotation_x <- x_range[1] + c(0.2, 0.8) * diff(x_range)
  annotation_y <- y_range[1] + c(0.8, 0.75) * diff(y_range)

  # create a volcano plot for the current comparison
  volcano <- ggplot(
    result_table,
    aes(x = log2FoldChange, y = -log10(pvalue))
  ) +
    geom_point(aes(color = DE), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c(
        "Upregulated" = "blue",
        "Downregulated" = "red",
        "Not significant" = "gray70"
      )
    ) +
    labs(
      title = "DESeq2 volcano plot",
      subtitle = paste(treatment_group, "vs", reference),
      x = "log2 fold change",
      y = expression(-log[10](p-value)),
      color = "Differential expression"
    ) +
    annotate(
      geom = "text", x = annotation_x, y = rep(annotation_y[1], 2),
      label = c("Down", "Up"), size = 5
    ) +
    annotate(
      geom = "text", x = annotation_x, y = rep(annotation_y[2], 2),
      label = c(down, up), size = 5
    ) +
    theme_bw()

  ggsave(
    filename = file.path(output_dir, paste0("volcano_", safe_comparison, ".jpeg")),
    plot = volcano,
    height = 5,
    width = 8,
    dpi = 300
  )

  summary(result)
  cat("Wrote", safe_comparison, "results to", output_dir, "\n")
}





