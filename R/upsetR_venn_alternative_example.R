###############################################################################
# DESCRIPTION
# This R file, takes in CSV files of normalized counts, runs DESeq2 on them,
# and calculates significant upregulated genes and plots intersection diagram.
#
# AUTHOR
# Akshay Paropkari
#
# VERSION
# 0.0.2
###############################################################################

# library imports
library(DESeq2)
library(tools)
library(UpSetR)

# note that all normalized count files need to be in one folder
# rename normalized count files such as TF_normalized_count.csv
# i.e. bcr1_normalized_count.csv
norm.count.files <- list.files(path = "~/Desktop/input_folder",
                               pattern = "*_normalized_counts.csv",
                               full.names = T)

# initiate list to save all significantly upregulated ORFs
final.list <- list()

# iterate through all files in input folder, and calculate significantly
# upregulated ORFs
for (f in norm.count.files) {
  mutant <- unlist(strsplit(x = unlist(strsplit(x = f,
                                                split = "/"))[6],
                            split = "_"))[1]
  message(paste0("FORMATTING ", toTitleCase(mutant),"-KO DATA SET"))
  norm.data <- read.csv(file = f, header = T)
  norm.data <- norm.data[grep(pattern = "^orf*", x = norm.data$X), ]
  norm.data <- norm.data[rowMeans(x = norm.data[, 2:ncol(norm.data)]) > 10, ]

  column.names <- c("ORF.ids", "WT.SN250.1", "WT.SN250.2",
                    paste(mutant, "KO", "1", sep = "."),
                    paste(mutant, "KO", "2", sep = "."))
  colnames(norm.data) <- column.names
  row.names(norm.data) <- norm.data$ORF.ids
  norm.data$ORF.ids <- NULL
  norm.data <- round(x = norm.data)
  norm.data <- as.matrix(norm.data)
  norm.coldata <- data.frame(condition=c("WT", "WT", "TF.KO", "TF.KO"),
                             row.names = column.names[2:5])

  # run DESeq2 nbinom analysis for mutant
  message(paste0("RUNNING DESeq2 ON ", toTitleCase(mutant),"-KO DATA SET"))
  dds <- DESeqDataSetFromMatrix(countData = norm.data,
                                colData = norm.coldata,
                                design = ~ condition, )
  dds$condition <- factor(dds$condition, levels = c("WT","TF.KO"))
  dds <- DESeq(object = dds)
  res <- results(object = dds)
  resLFC <- lfcShrink(dds = dds,
                      coef = "condition_TF.KO_vs_WT",
                      type = "apeglm")

  # remove NAs from DESeq2 result DF
  res <- res[complete.cases(res), ]

  # diagnostic MA plots for visual inspection, user can comment out lines 47-50
  plotMA(object = resLFC, y = c(-5, 5),
         main = paste(mutant, "KO", "vs", "WT", sep = " "))
  abline(h = c(-1, 1), col = "dodgerblue", lwd=3)
  # plotDispEsts(object = dds, main = paste(mutant, "KO", "vs", "WT", sep = " "))

  # collect significant upregulated genes in mutant
  message(paste0("COLLECTING SIGNIFICANTLY UPREGULATED ORFs IN ",
                 toTitleCase(mutant),"-KO DATA SET"))
  orf.vec <- row.names(res)[which(res$padj < 0.05 & res$log2FoldChange > 0)]

  # save ORF list into a list
  message(paste0("SAVING ORF LIST FOR ", toTitleCase(mutant), "-KO DATA SET"))
  name <- paste0(toTitleCase(mutant), "-KO")
  final.list[[name]] <- orf.vec
  message(paste0("=========================================================="))
}

# plot upset plot and save as  JPEG file on the Desktop
jpeg(filename = "~/Desktop/sig_upreg_orfs_intersection_plot.jpeg",
     width = 21, height = 9, units = "in", res=300)
upset(data = fromList(final.list), nsets = 100, nintersects = NA,
      order.by = "freq", point.size = 4,
      mainbar.y.label = "Intersection set size",
      sets.x.label = "Set size",
      text.scale = c(2.5, 0, 2, 1.3, 2, 1.75))
dev.off()
