#--------------------------------------------------------------------------------------------------------
# Generate the bed files
# Target regions are 2kb upstream and downstream the TSS per gene --> we are plotting 4kb for all genes
#--------------------------------------------------------------------------------------------------------
rm(list = ls(all=TRUE))

# Working and output directories
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
all_bed <- file.path(wd, "bed_files")
in_bed = file.path(wd, "bed_files", "metagene_BED_gene_TAIR10.bed")
in_genes = file.path(wd, "gene_sets")
out_data = file.path(wd, "bed_files")

# empty directory "bed_files" (except "metagene_BED_gene_TAIR10.bed")
# to avoid mixing up outputs from subsequent runs
unlink(file.path(all_bed, list.files(all_bed)[list.files(all_bed) != "metagene_BED_gene_TAIR10.bed"])) # remove all but one bed files

#--------------------------------------------------------------------------------------------------------
# Generate the bed files for groups of genes in separate files in the directory "gene_sets";
# a for-loop is used
#--------------------------------------------------------------------------------------------------------

# Load the bed file for the entire genome
filename <- in_bed
bed <- read.table(filename,sep="\t",stringsAsFactors=FALSE, header=FALSE)
dim(bed) # 28770     6

# load gene sets
gene_sets <- list.files(path = in_genes)


# make bed files for each gene set
for (i in 1:length(gene_sets)){
  # load the gene set
  filename <- file.path(in_genes, gene_sets[i])
  gene_set_basename <- substr(gene_sets[i], start = 1, stop = regexpr("\\.txt", gene_sets[i])[1]-1)
  gene_set <- read.table(filename, sep="\t",stringsAsFactors=FALSE, header=FALSE)
  gene_set <- gene_set[grep("AT", gene_set[,1]),1] # select rows only with the TAIR gene code
  
  # substract gene codes without gene models (ATxXxxxxx vs ATxXxxxxx.x)
  for (k in 1:length(gene_set)){
    start_gene_code <- regexpr("AT", gene_set[k])
    end_gene_code <- start_gene_code + 9
    if (nchar(gene_set[k]) >= 9) {
      gene_set[k] <- substr(gene_set[k], start = start_gene_code, stop = end_gene_code)
    } else (print("check all gene codes in the provided set of genes are in the format ATxGxxxxx"))
  }
  
  print(paste("nr. genes in gene set", gene_set_basename, "-", length(gene_set)))
  gene_set <- as.vector(sort(gene_set))
  
  # prepare bed file for the gene set
  bed_gene_set <- bed[which(bed$V4 %in% gene_set),]
  print(paste("nr. genes found in the TAIR10 BED file -", nrow(bed_gene_set)))
  
  # Write bed files
  filename <- file.path(out_data, paste0(gene_set_basename, ".bed"))
  write.table(bed_gene_set, filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)
  
  # clean-up
  rm(filename, gene_set_basename, gene_set, bed_gene_set)
}

# make a bed file for a "random" subset of 2000 genes
# use of the whole is a heavy burden on the memory
bed_gene_set <- bed[sample(1:nrow(bed), size = 2000, replace = FALSE),]
filename <- file.path(out_data, paste0("TAIR_2000", ".bed"))
write.table(bed_gene_set, filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)
