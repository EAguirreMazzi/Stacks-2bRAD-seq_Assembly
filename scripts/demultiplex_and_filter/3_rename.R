rm(list = ls())
plate_map <-read.table("0_sequence_processing_and_assembly/plate_info/Plate_maps_master.txt", header = T, sep = "\t")
#check for duplicated in Sequence_ID and rename accordingly (these are replicates)
sum(duplicated(plate_map$Sequence_ID)) == 0

plate_map <- plate_map[plate_map$is_weinmannia,]
#filename after running the trim2bRAD_2barcodes_noAdap.pl
plate_map$demux_file <- paste0(plate_map$plate, ifelse(plate_map$bc_Label_file == "n/a","",paste0("_",plate_map$bc_Label_file)), "_row",plate_map$row,"_",plate_map$column_antiBC, ".fq")

renaming_cmd<- paste0("mv ", plate_map$demux_file, " ", plate_map$Sequence_ID, ".fastq")
writeClipboard(renaming_cmd)

#After running command there are some sequences not recovered (probably had less than 10000 reads and file was not created)
#We will figure out which were these when getting fastp reports