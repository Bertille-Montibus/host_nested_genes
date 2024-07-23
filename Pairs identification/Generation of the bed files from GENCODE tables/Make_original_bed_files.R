# Generation of the BED files to overlap to determine the list of Host and Nested genes #

##################
#### HUMAN #######
##################
library(dplyr)

##########
#Make a bed file for the host genes
##########
bed_host <- read.table('wgEncodeGencodeComprehensiveV36lift37.txt', header = F, sep = "\t")

#Remove the column with CDS_start and reorder the columns to make a bedtools
bed_host <- bed_host[,c(2,4,5,1,7,3,8)]

#Filtering out the biotypes associated with VDJ genes and to be experimentally confirmed (TEC)
#List of the biotypes
unique(bed_host$V8)
#[1] "protein_coding"          "processed_transcript"    "nonsense_mediated_decay"
#[4] "retained_intron"         "lncRNA"                  "lincRNA"
#[7] "sense_intronic"          "snRNA"                   "miRNA"
#[10] "misc_RNA"                "TEC"                     "snoRNA"
#[13] "non_stop_decay"          "polymorphic_pseudogene"  "antisense"
#[16] "sense_overlapping"       "rRNA"                    "scRNA"
#[19] "IG_V_gene"               "IG_C_gene"               "IG_J_gene"
#[22] "vault_RNA"               "TR_C_gene"               "TR_J_gene"
#[25] "TR_V_gene"               "TR_D_gene"               "Mt_tRNA"
#[28] "Mt_rRNA"                 "IG_D_gene"
#Use the index in the vector to make a filter
filter_biotypes <- unique(bed_host$V8)[c(1:10,12:18,22)]

#Filter the transcripts with the selected biotypes
bed_host = bed_host[bed_host$V8 %in% filter_biotypes,]

#Removal of the complex loci form the list
#Protocadherins
bed_host <- bed_host[!grepl("PCDH", bed_host$V7),]
#UGT (UDP-Glucuronosyltransferases)
bed_host <- bed_host[!grepl("UGT", bed_host$V7),]

#Save the bed file which will be used for the hosts
write.table(bed_host, "wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed", quote=F,col.names = F, row.names = F, sep = '\t')

##########
#Make a bed file for the nested genes
##########
#For the nested genes to make sure that all the transcripts of the nested gene will be included inside the host, recover the most extreme coordinates
#for each genes
Gene_names = unique(bed_host$V7)
New_bed = data.frame(chr = character(),
                     start = numeric(),
                     end = numeric(),
                     name = character())

for (i in Gene_names) {
  print(i)
  Gene_table = bed_host[bed_host$V7 == i,]
  print(Gene_table)
  start = min(Gene_table$V4)
  print(start)
  end = max(Gene_table$V5)
  table = data.frame(chr = Gene_table[1,1],
                     start = start,
                     end = end,
                     name = i)
  print(table)
  New_bed = rbind(New_bed,table)
}

#Save the bed file which will be used for the host
write.table(New_bed, "wgEncodeGencodeComprehensiveV36lift37_filefornested_filtered.bed", quote=F,col.names = F, row.names = F, sep = '\t')

##################
#### MOUSE #######
##################

##########
#Make a bed file for the host genes
##########
bed_host <- read.table('wgEncodeGencodeCompVM25.txt', header = F, sep = "\t")

#Remove the column with CDS_start and reorder the columns to make a bedtools
bed_host <- bed_host[,c(2,4,5,1,7,3,8)]

#Filtering out the biotypes associated with VDJ genes and to be experimentally confirmed (TEC)
#List of the biotypes
unique(bed_host$V8)
#[1] "protein_coding"                "lincRNA"
#[3] "sense_overlapping"             "antisense"
#[5] "processed_transcript"          "TEC"
#[7] "snRNA"                         "sense_intronic"
#[9] "miRNA"                         "snoRNA"
#[11] "misc_RNA"                      "rRNA"
#[13] "ribozyme"                      "scaRNA"
#[15] "polymorphic_pseudogene"        "bidirectional_promoter_lncRNA"
#[17] "macro_lncRNA"                  "3prime_overlapping_ncRNA"
#[19] "TR_V_gene"                     "TR_D_gene"
#[21] "TR_J_gene"                     "TR_C_gene"
#[23] "IG_LV_gene"                    "IG_V_gene"
#[25] "IG_J_gene"                     "IG_C_gene"
#[27] "sRNA"                          "scRNA"
#[29] "Mt_tRNA"                       "Mt_rRNA"
#[31] "IG_D_gene"
#Use the index in the vector to make a filter
filter_biotypes <- unique(bed_host$V8)[c(1:5,7:18,27:30)]

#Filter the transcripts with the selected biotypes
bed_host = bed_host[bed_host$V8 %in% filter_biotypes,]


bed_host <- bed_host[!grepl("Ugt", bed_host$V7),]

#Removal of the complex loci form the list
#Protocadherins
bed_host <- bed_host[!grepl("Pcdh", bed_host$V7),]
#UGT (UDP-Glucuronosyltransferases)
bed_host <- bed_host[!grepl("Ugt", bed_host$V7),]

#Save the bed file which will be used for the hosts
write.table(bed_host, "wgEncodeGencodeCompVM25_fileforhost_filtered.bed", quote=F,col.names = F, row.names = F, sep = '\t')

##########
#Make a bed file for the nested genes
##########
#For the nested genes to make sure that all the transcripts of the nested gene will be included inside the host, recover the most extreme coordinates
#for each genes
Gene_names = unique(bed_host$V7)
New_bed = data.frame(chr = character(),
                     start = numeric(),
                     end = numeric(),
                     name = character())

for (i in Gene_names) {
  print(i)
  Gene_table = bed_host[bed_host$V7 == i,]
  print(Gene_table)
  start = min(Gene_table$V4)
  print(start)
  end = max(Gene_table$V5)
  table = data.frame(chr = Gene_table[1,1],
                     start = start,
                     end = end,
                     name = i)
  print(table)
  New_bed = rbind(New_bed,table)
}

#Save the bed file which will be used for the host
write.table(New_bed, "wgEncodeGencodeCompVM25_filefornested_filtered.bed", quote=F,col.names = F, row.names = F, sep = '\t')
