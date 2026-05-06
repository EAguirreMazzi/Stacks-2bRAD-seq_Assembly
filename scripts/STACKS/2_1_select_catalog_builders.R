rm(list = ls())
source("3_evolution/scripts/database/1_load_db.R")
db <- db[db$include_analysis == TRUE,]
#now is a good poing to export general popmap
write.table(cbind(db$Sequence_ID,rep("defaultpop",nrow(db))),
     row.names = F, col.names = F, quote = F, sep = "\t",
     file = paste0("0_sequence_processing_and_assembly/assembly_info/raw_full_popmap.txt"))

depth_stats<- c("loci_stacks", "mean_coverage", "sd_coverage", 
                  "max_coverage", "n_reads", "percentage_reads")

pdf("0_sequence_processing_and_assembly/assembly_info/ustacks_depths.pdf", width = 14)
for(stat in depth_stats){
    boxplot(db[,stat] ~ db$plate, xlab= "Plate", ylab =stat)
}

dev.off()

boxplot(db$max_coverage)
boxplot(db$mean_coverage)
summary(db$mean_coverage)

#select samples for catalog

# 1st criteria Number of Putative Loci
# Good indicator of:
# - How many BcgI sites were successfully sequenced
# - Library complexity
# - Whether regex filtering worked correctly

# Interpretation:
# Too few loci → poor library, low coverage, failed digestion
# Too many loci → repeat elements, collapsed paralogs,
#                 assembly artifacts, or regex too permissive
# SELECT: samples near the median-upper range
# EXCLUDE: extreme outliers in both directions
# samples with suspiciously HIGH locus counts may
# introduce artifactual loci into the catalog
ustacsk_loci<-db$loci_stacks
#remove outlier loci and compute
selected_by_loci_number <- which((ustacsk_loci >= quantile(ustacsk_loci,0.5) & ustacsk_loci <= quantile(ustacsk_loci,0.90)))

catalog_builders <- db[selected_by_loci_number,]

# 2nd criteria mean coverage
#High mean coverage ≠ high quality
#It may reflect:
#- PCR duplication (bad)
#- Few loci at deep coverage (bad)
#- Genuine good sequencing (good)
#Use mean coverage as a FLOOR not a ceiling:
#SELECT: samples above minimum viable mean (e.g. >5x)
#DO NOT: select purely highest mean coverage samples
#        → biases catalog toward NovaSeq/high-PCR samples
selected_by_mean_coverage <- which(catalog_builders$mean_coverage >=5)
catalog_builders <- catalog_builders[selected_by_mean_coverage,]
dim(catalog_builders)
#3rd criteris SD
# Low SD relative to mean (CV < 1):
# → Coverage relatively uniform across loci
# → Good library complexity
# → Reliable catalog builder
# High SD relative to mean (CV >> 1):
# → Few loci with extreme depth, many with low depth
# → Repeat elements or PCR artifacts inflating variance
# → Poor catalog builder despite possibly high mean
# Calculate coefficient of variation (CV = SD/mean):
# This is more informative than SD alone
#we will be more lenient than one in this filter just to remove disastrous samples
selected_by_CV <- which(catalog_builders$sd_coverage/catalog_builders$mean_coverage <= 3)
catalog_builders <- catalog_builders[selected_by_CV,]
dim(catalog_builders)
#4th criteria
#Highly informative for identifying problematic samples:
#Very high max coverage indicates:
# → Repeat/mitochondrial loci present
# → Collapsed paralogs
# → PCR duplication artifacts
# These loci will enter your catalog and create
# spurious matches across all 1000 samples
# EXCLUDE samples where:
# max_coverage > 10 × mean_coverage
# (suggests at least one severely inflated locus)
selected_by_maxCov <- which(catalog_builders$max_coverage <= quantile(catalog_builders$max_coverage, .90))
catalog_builders <- catalog_builders[selected_by_maxCov,]

#4th criteria make sure to include stuff from divergent lineages (we will use region as the broadest proxy)
#also include w. boliviensis (southern bolivia representative)
# outgroups and sparse representatives of other regions outside central and north andes and costa RIca which are well represented
divergent_country <-setdiff(unique(db$CountryName),unique(catalog_builders$CountryName))
divergent_region_specimens <- db[db$CountryName %in% divergent_country | db$Species == "boliviensis",]
#remove if they are really disastrous only
divergent_region_specimens[,c("Sequence_ID",depth_stats)]
#from these will remove those with very few loci (unusccesfull sequencing)
divergent_region_specimens <- divergent_region_specimens[divergent_region_specimens$loci_stacks >= 5000,]

#merge with catalog builders
catalog_builders <- rbind(catalog_builders, divergent_region_specimens)

#Explore how regions are sampled
counts_country <-aggregate(catalog_builders$Sequence_ID, by=list(catalog_builders$CountryName), FUN=length)
#I further downsample Bolivia (there are 155) to ~2* each species

catalog_builders_to_downsample <- catalog_builders[catalog_builders$CountryName == "Bolivia",]
set.seed(2021)
to_keep <- lapply(split(catalog_builders_to_downsample, 
                           catalog_builders_to_downsample$Species),
                     function(x) {
                       if(nrow(x) >= 2) {
                         x[sample(nrow(x), 2), ]
                       } else {
                         NULL  # Skip species with <2 samples
                       }
                     })
to_keep <- do.call(rbind, to_keep);rownames(to_keep) <- NULL
to_downsample<-setdiff(catalog_builders$Sequence_ID[catalog_builders$CountryName == "Bolivia"], to_keep$Sequence_ID)
catalog_builders <- catalog_builders[!catalog_builders$Sequence_ID %in% to_downsample,]
to_keep$Sequence_ID %in% catalog_builders$Sequence_ID #check approach worked

#export stacje metric for builders

pdf("0_sequence_processing_and_assembly/assembly_info/catalog_builders_ustacks_depths.pdf", width = 18)
for(stat in depth_stats){
    boxplot(catalog_builders[,stat] ~ catalog_builders$CountryName, xlab= "Country", ylab =stat)
}
dev.off()

dim(catalog_builders)
#also shuffle the catalog builders to void bias in how loci can be dicovered/grouped
set.seed(1502)
catalog_builders <- catalog_builders[permute::shuffle(catalog_builders),]

#SINCE IM USING ONLY 94 final catalog builder not staging this time required
#stages_1to3<-split(catalog_builders,cut(seq_along(catalog_builders),3, labels=FALSE))
#for( i in seq_along(stages_1to3)){
#    write.table(cbind(stages_1to3[[i]],rep("defaultpop",length(stages_1to3[[i]]))),
#     row.names = F, col.names = F, quote = F, sep = "\t",
#            file = paste0("3_evolution/scripts/assembly_prep/stacks/popmaps/catalog_builders_stage",i,"_popmap.txt"))
#}



write.table(cbind(catalog_builders$Sequence_ID,rep("defaultpop",length(catalog_builders$Sequence_ID))),
     row.names = F, col.names = F, quote = F, sep = "\t",
     file = paste0("0_sequence_processing_and_assembly/assembly_info/catalog_builders_popmap.txt"))

#also write string to add in cstacks commans with -i -o options (useful for branching different n values)
write.table(paste(paste0("-s ./", catalog_builders$Sequence_ID), collapse = " "),
     row.names = F, col.names = F, quote = F, sep = "\t",
     file = paste0("0_sequence_processing_and_assembly/assembly_info/catalog_builders_string.txt"))
