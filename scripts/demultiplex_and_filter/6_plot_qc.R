rm(list = ls())
library(ggplot2)

#From fastP multiqc
gc_content <- "0_sequence_processing_and_assembly/plate_info/fastp_reports/multiqc_data/fastp-seq-content-gc-plot_Read_1_After_filtering.txt"
# Read the file
gc_content <- read.csv(gc_content, sep = "\t", header = TRUE)
values<- do.call(cbind.data.frame, lapply(gc_content[,2:37], FUN=function(x) as.numeric(stringr::str_extract(x, "(?<=, ).*?(?=\\))"))))
gc_content[,2:37] <- values
gc_content<-reshape2::melt(gc_content,value.name = "gc_content")
names(gc_content)[2] <- "site"
gc_content$site <- as.numeric(gsub("X","",gc_content$site))

# View result
pdf("0_sequence_processing_and_assembly/plate_info/gc_content_per_plate_after__fasstp_filtering.pdf", height = 14, width = 8)
ggplot(gc_content)+
    geom_line(aes(x=site, y=gc_content, group=Sample), color="#0478c5ad")+
    facet_grid(plate ~ .)
dev.off()