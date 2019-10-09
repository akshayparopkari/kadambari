###############################################################################
# DESCRIPTION
# This R file, reads in output from DESeq2 run and plots UpSet diagram showing
# significant and differentially expressed genes in multiple conditions.
#
# AUTHOR
# Akshay Paropkari
#
# VERSION
# 0.0.4
###############################################################################

# library imports
library(DESeq2)
library(tools)
library(UpSetR)

# note that all normalized count files need to be in one folder
# rename normalized count files such as TF_normalized_count.csv
# i.e. bcr1_normalized_count.csv
input.folder <- "~/Desktop/input_folder"  # CHANGE THIS TO YOUR INPUT FOLDER LOCATION
deseq.res.files <- list.files(path = input.folder,
                              pattern = "*_deseq2_output.tsv", full.names = T)

# initiate list to save all significant ORFs
final.list <- list()

# iterate through all files in input folder, and calculate significantly
# upregulated ORFs
for (f in deseq.res.files) {
  mutant <- unlist(strsplit(x = unlist(strsplit(x = f,
                                                split = "/"))[6],
                            split = "_"))[1]
  message(paste0("FORMATTING ", toTitleCase(mutant),"-KO DATA SET"))
  deseq.data <- read.csv(file = f, sep = "\t", header = T)

  # consider entries with ORF ids only
  deseq.data <- deseq.data[grep(pattern = "^orf*", x = rownames(deseq.data)), ]

  # remove ORFs with NA values in any column
  deseq.data <- deseq.data[complete.cases(deseq.data), ]

  # collect significant upregulated genes in mutant

  # significant and upregulated ORFs
  # message(paste0("COLLECTING SIGNIFICANTLY UPREGULATED ORFs IN ",
  #                toTitleCase(mutant),"-KO DATA SET"))
  # expr <- "upregulated"
  # orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05 & deseq.data$log2FoldChange > 0)]

  # significant and downregulated ORFs
  # message(paste0("COLLECTING SIGNIFICANTLY DOWNREGULATED ORFs IN ",
                 # toTitleCase(mutant),"-KO DATA SET"))
  # expr <- "downregulated"
  # orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05 & deseq.data$log2FoldChange < 0)]

  # significant ORFs
  message(paste0("COLLECTING SIGNIFICANT AND DIFFERENTIALLY EXPRESSED ORFs IN ",
                 toTitleCase(mutant),"-KO DATA SET"))
  expr <- "differentially_expressed"
  orf.vec <- row.names(deseq.data)[which(deseq.data$padj < 0.05)]

  # save ORF list into a list
  message(paste0("SAVING ORF LIST FOR ", toTitleCase(mutant), "-KO DATA SET"))
  name <- paste0(toTitleCase(mutant), "-KO")
  final.list[[name]] <- orf.vec
  message(paste0("=========================================================="))
}

# save ORF list to file
fnh <- paste0(input.folder, "/", expr, "_orf_list.csv")
# lapply(final.list,
#        function(x) write.table(x = data.frame(x),
#                                file = fnh,
#                                append= T, sep='\t' ))


# plot upset plot and save as  JPEG file on the Desktop
# NO SPACES in file name or output folder
image.fnh <- paste0(input.folder, "/", expr, "_orf_upset_plot.jpeg")
jpeg(filename = image.fnh, width = 21, height = 9, units = "in", res=300)
upset(data = fromList(final.list), nsets = 100, nintersects = NA,
      order.by = "freq", point.size = 4,
      mainbar.y.label = "Intersection set size",
      sets.x.label = "Set size",
      text.scale = c(2.5, 0, 2, 1.3, 2, 1.75))
dev.off()
