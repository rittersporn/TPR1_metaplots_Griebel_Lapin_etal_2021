#--------------------------------------------------------------------------------------------------------
# Generate metaplots
# Target regions are 2kb upstream and downstream the TSS per gene --> we are plotting 4kb for all genes
#--------------------------------------------------------------------------------------------------------


#BiocManager::install("metagene")
library("metagene") # tested versions 2.18.0 and 2.22.0
library("BiocGenerics") # tested versions 0.32.0 and 0.36.1
library("stringr") # tested version 1.4.0
library("knitr") # tested version 1.30, 1.33
library("ggplot2") # tested versions 3.3.2 and 3.3.5


# Working and output directories
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
in_bed = file.path(wd, "bed_files")
in_bam = file.path(wd, "bam_bai", "TPR1")
out_data = file.path(wd, "metaplots")


# empty directory "metaplots"
# to avoid mixing up outputs from subsequent runs
unlink(file.path(out_data, "*")) # remove all files from the directory "metaplots"


#--------------------------------------------------------------------------------------------------------
# Metaplotting
#--------------------------------------------------------------------------------------------------------

#Load in regions in .bed format and remove the bed file for the entire CDS set - not enough memory
regions <- list.files(path=in_bed, full.names=FALSE, pattern=".bed$")
regions <- regions[grep("BED_gene", regions, invert = TRUE)]

# get short names of regions
regions_short <- regions
for (k in 1:length(regions)){
  regions_short[k] <- substr(regions[k], start = 1, stop = regexpr("\\.bed", regions[k])[1]-1)
}

# give full path for regions (gene sets) to display on the plot
regions_full <- list.files(path=in_bed, full.names=TRUE, pattern=".bed$")
regions_full <- regions_full[grep("BED_gene", regions_full, invert = TRUE)]



#--------------------------------------------------------------------------------------------------------
# Define mark name for plotting and BAM file
#--------------------------------------------------------------------------------------------------------
full_mark_name_and_location <- data.frame(mark_name = c("ChIP_TPR1_Col",
                                                        "input_TPR1_Col",
                                                        "ChIP_TPR1_eds1",
                                                        "input_TPR1_eds1"),
                                          bam_location = c("TPR1_Col_ChIP.bam$",
                                                           "TPR1_Col_input.bam$",
                                                           "TPR1_eds1_ChIP.bam$",
                                                           "TPR1_eds1_input.bam$"))

complete_genotype_set <- c("TPR1_Col", "TPR1_eds1")
#--------------------------------------------------------------------------------------------------------
# generate metaplots
#--------------------------------------------------------------------------------------------------------
for (genotype in complete_genotype_set){
  # select mark name and ChIP file names and  for the genotype of interest
  mark_name_and_location <- full_mark_name_and_location[grep(genotype, full_mark_name_and_location$mark_name), ]
  
  for (i in 1:nrow(mark_name_and_location)){
    #Name of the dataset to appear on plots
    mark_name <- as.character(mark_name_and_location$mark_name[i])
    bam_name <- as.character(mark_name_and_location$bam_location[i])
    
    print(paste0("metaplotting for ", mark_name))
    #Load BAM files in a vector format (should be in the same directory as the index files)
    bam_files<-list.files(path=in_bam, pattern=bam_name, full.names=TRUE)
    bam_files
    
    #mapping region coordinates against .BAM file and retrieving raw coverage
    mg <- metagene$new(regions = regions_full, bam_files = bam_files, verbose = TRUE)
    
    #Make sure to flip the genes that are on the -strand
    mg$produce_table(flip_regions = TRUE, normalization = 'RPM')
    mg$produce_data_frame(alpha = 0.01, sample_count = 1000, avoid_gaps = TRUE)
    df <- mg$get_data_frame()
    
    #rename groups to simplify legend on plots
    group_names = as.vector(df$group)
    for (k in 1:length(group_names)){
      for (n in 1: length(regions_short))
        group_names[k][grep(regions_short[n], group_names[k])] <- regions_short[n]
    }
    df$group <- as.factor(group_names)
    
    #plot the data
    my_plot = plot_metagene(df) + labs(y = "Coverage, RPM", x = "", title = mark_name) +
      coord_cartesian(ylim = c(0, 3)) +
      scale_x_continuous(breaks = c(1,50,100), labels = c("-2 kb","TSS","2 kb")) +
      theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
            axis.title.x = element_text(color="blue", size=14, face="bold"),
            axis.title.y = element_text(size=14))
    ggsave(file.path(out_data, paste0(mark_name,".pdf")))
    rm(my_plot)
    # save dataframes for making plots with basic normalization
    if(grepl("ChIP", mark_name)) df_ChIP <- df
    if(grepl("input", mark_name)) df_input <- df
  }
  
  
  # prepare a df with norm read data - very basic because normalization is performed for genes on average
  # not per gene
  df_norm <- df_ChIP[, c("region", "bin", "value")]
  df_norm$value <- df_ChIP$value - df_input$value
  
  
  #plot the normalized data
  my_plot <- ggplot(data = df_norm, aes(x = bin, y = value, group = region, colour = region)) +
    geom_line(size = 3) +
    labs(y = "Coverage, RPM, input normalized", x = "",
         title = paste0("Enrichment of TPR1-GFP at the gene sets in ", genotype)) +
    coord_cartesian(ylim = c(0, 3)) +
    scale_x_continuous(breaks = c(1,50,100), labels = c("-2 kb","TSS","2 kb")) +
    theme_bw() + 
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
          axis.title.x = element_text(color="blue", size=14, face="bold"),
          axis.title.y = element_text(size=14))
  ggsave(file.path(out_data, paste0(genotype, "_input_norm.pdf")))
  rm(my_plot)
}

# print info about used packages
sessionInfo()
