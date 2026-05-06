
rm(list = ls())
#load plate map master and QC post filtering and normalization of orientation
plate_map <- read.table("0_sequence_processing_and_assembly/plate_info/Plate_maps_master.txt", header = T, sep = "\t")
reads <-read.table("0_sequence_processing_and_assembly/plate_info/multiqc_fastqc_after_strand_normalization_data/multiqc_fastqc.txt",
     header = T, sep = "\t")[,c(1,5,9,10)]
sequence_count <-read.table("0_sequence_processing_and_assembly/plate_info/multiqc_fastqc_after_strand_normalization_data/fastqc_sequence_counts_plot.txt",
    header = T, sep = "\t")

if(all(reads$Sample == sequence_count$Sample)){
    reads<-cbind.data.frame(reads,sequence_count[,-1])
}
#to simplyfy work will just remove things that werent "succesfully sequenced"
plate_map<- plate_map[plate_map$Sequence_ID %in% reads$Sample,]
#
matching_index<-match(plate_map$Sequence_ID, reads$Sample)
if(all(reads$Sequence_ID[matching_index] == plate_map$Sequence_ID)){
 plate_map<- cbind.data.frame(plate_map, reads[matching_index,-1])
}

#compute uniqueness
plate_map$Unique.Reads  + plate_map$Duplicate.Reads == plate_map$Total.Sequences
plate_map$uniqueness  <- plate_map$Unique.Reads/plate_map$Total.Sequences
plate_map$depth_proxy <- plate_map$Total.Sequences/plate_map$Unique.Reads


#merge technical replicates if they are from same platform otherwise move to low_reads
plate_map <- plate_map[-which(plate_map$Sequence_ID == "cubr_40644_1"),] #ALSO REMOVE THIS
plate_map$Project_ID <- gsub("_\\d$","",plate_map$Sequence_ID)
plate_map$include_analysis <- TRUE
dup_seqs <-aggregate(plate_map$Sequence_ID, by=list(plate_map$Project_ID), FUN=length)
dup_seqs<- dup_seqs[!dup_seqs$x == 1,]


select_commands<- c("mkdir -p unused_seqs", "rm cubr_40644_1.fastq.gz")
merged_samples<- c()

for(seq in 1:nrow(dup_seqs)){
    tmp<-plate_map[plate_map$Project_ID %in% dup_seqs$Group.1[seq],]
    if(length(unique(tmp$method))==1){
       select_commands<-c(select_commands,paste("cat" ,paste0(tmp$Sequence_ID,".fastq.gz", collapse = " "),">", paste0(tmp$Project_ID[1], "_8.fastq.gz"), collapse = " "))
       select_commands<-c(select_commands,paste("mv" ,paste0(tmp$Sequence_ID,".fastq.gz"), "unused_seqs"))
       merged_samples<- c(merged_samples,tmp$Project_ID[1])
       #update db
       new_index <- nrow(plate_map)+1
       plate_map$include_analysis[plate_map$Sequence_ID %in% tmp$Sequence_ID] <- FALSE
       plate_map[new_index,] <- tmp[1,]
       plate_map$Sequence_ID[new_index] <- paste0(tmp$Project_ID[1], "_8")
       plate_map$plate[new_index] <- "merged_duplicates"
       plate_map$Total.Sequences[new_index] <- sum(tmp$Total.Sequences)
       plate_map$include_analysis[new_index] <- TRUE
       plate_map[new_index, c("plate", "barcoded_PCR_oligo", "row", "column", "method",  "X.GC", "total_deduplicated_percentage", "Unique.Reads", "Duplicate.Reads", "uniqueness", "depth_proxy")]<- NA
    }
}
#run this command in the file path
# writeClipboard(select_commands)
sum(plate_map$include_analysis )

#exclude samples that will tend to have low depth
cuttoff_reads <- 20000
#an alternative filter to guide towards ~3x coverage
#plate_map[!plate_map$Unique.Reads*3 < plate_map$Total.Sequences, c("Sequence_ID","Total.Sequences","Unique.Reads","uniqueness")]
plate_map$include_analysis[plate_map$Total.Sequences < cuttoff_reads] <- FALSE
sum(plate_map$include_analysis)

#move samples not to include to another folder
plate_map[!plate_map$include_analysis, c("Sequence_ID","Total.Sequences","Unique.Reads","uniqueness")]

#writeClipboard(paste0("mv ",plate_map$Sequence_ID[!plate_map$include_analysis], ".fastq.gz unused_seqs/"))
#also move unwanted sequences like: "cubr_40644_1.fastq.gz"



plate_map_clean <- plate_map[plate_map$include_analysis,]
#THIS TIME NORMALIZATION IS NOT NEEDED FOR PLATE 16 But to improve assembly I will normalize all plates bases on quantile
boxplot(plate_map_clean$Total.Sequences)
plot(quantile(plate_map_clean$Total.Sequences, seq(0,1, 0.01)))

#bases on this will cap rather than normalize(since distribution is lognormal but with outliers)
# the idea is to attenuate the effect of outliers on loci discovery (paralog and sequencing error related)
cuttoff<- quantile(plate_map_clean$Total.Sequences, 0.95)
db_downsample <- plate_map_clean[plate_map_clean$Total.Sequences > cuttoff,]
db_downsample$fraction_to_keep <- round(cuttoff/db_downsample$Total.Sequences,3)
boxplot(db_downsample$fraction_to_keep ~ db_downsample$plate)

#make command to downsample these specimens and move them a folder for assembly
mv_cmd <- c("mkdir -p before_downsampling_reads",paste0("mv ",db_downsample$Sequence_ID,".fastq.gz before_downsampling_reads/"))
#writeClipboard(mv_cmd)

#cd /storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw/normalized/before_downsampling_reads
seqtk_command<- paste0("seqtk sample -s100 ",
    db_downsample$Sequence_ID,".fastq.gz ", db_downsample$fraction_to_keep,
    " | gzip > /storage1/fs1/christine.e.edwards/Active/weinmannia/data/all_raw/normalized/", db_downsample$Sequence_ID, ".fastq.gz")
#writeClipboard(seqtk_command)
#put on a script and run:
#bsub -G compute-christine.e.edwards -q general -a "docker(staphb/seqtk)" -R "rusage[mem=16GB]" ./normalize.sh
#NOW samples are ready for ustacks
#export update sequence info
plate_map$Normalized_Total.Sequences <- plate_map$Total.Sequences

plate_map$Normalized_Total.Sequences[plate_map$Sequence_ID %in% db_downsample$Sequence_ID] <- cuttoff


#visualize number of reads and GC content per plate
pdf("0_sequence_processing_and_assembly/plate_info/filtered_seq_stats.pdf", width = 14)
boxplot(plate_map$Total.Sequences~plate_map$plate)
boxplot(plate_map$Normalized_Total.Sequences~plate_map$plate)
boxplot(plate_map$X.GC ~plate_map$plate)
boxplot(plate_map$uniqueness ~plate_map$plate)
boxplot(plate_map$depth_proxy ~plate_map$plate)
boxplot(plate_map$X.GC~plate_map$plate)
dev.off()

write.table(plate_map,
    row.names = F, quote = F, sep = "\t",
    file="0_sequence_processing_and_assembly/plate_info/Plate_maps_master_post_cleaning.txt")
