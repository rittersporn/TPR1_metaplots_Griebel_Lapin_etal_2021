#### script to select genes associated with TPR1 binding sites ####

# set working directory to the source directory for this file, this parameter will be adjusted later in the code

# define names of files with input peak information 
filename_input_peak_info <- c("TableS5_TPR1_Col_peak_annotation.txt",
                              "TableS7_TPR1_eds1_peak_annotation.txt")

# define names of output files with gene lists
filename_output_genesets <- c("TPR1_Col_bound_genes.txt",
                              "TPR1_eds1_bound_genes.txt")

# check where the working directory is
# and set it to directory two level up from the directory where this script is stored in downloaded repository
getwd()
setwd("./..")
getwd()

######## code below will extract geneIDs associated with Quest-called peaks #########
for (filename_order in 1:length(filename_input_peak_info)){
  # read TPR1 binding peaks
  peak_info <- read.delim(file.path("input_files",filename_input_peak_info[filename_order]),
                          sep = "\t", dec = ".", header = TRUE)
  
  # nr of peaks associated with unique annotated genes
  # using only Nearest Promoter ID or Annotation column will miss potential genes associated with the peak
  # the solution is to combine geneIDs from both columns
  # the difference between these two columns is observed when the peak has a TTS annotation
  peak_AGI_all <- unique(c(substr(peak_info$Nearest.PromoterID[peak_info$Distance.to.TSS >= -3000], start = 1, stop = 9),
                           substr(peak_info$Annotation[grep("AT", peak_info$Annotation)],
                                  start = regexpr("AT", peak_info$Annotation[grep("AT", peak_info$Annotation)]),
                                  stop = regexpr("AT", peak_info$Annotation[grep("AT", peak_info$Annotation)])+8)))
  peak_AGI_all <- peak_AGI_all[grep("ATM", peak_AGI_all, invert = TRUE)] # remove mitochondrial genes
  peak_AGI_all <- peak_AGI_all[grep("ATC", peak_AGI_all, invert = TRUE)] # remove chloroplast genes
  
  
  length(peak_AGI_all) # number of genes associated with TPR1
  
  # write geneIDs into the files
  write.table(peak_AGI_all[order(peak_AGI_all)], file.path("output_files", filename_output_genesets[filename_order]),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # clean-up to prevent carryover/crosscontamination of datasets
  rm(peak_info, peak_AGI_all)
}
