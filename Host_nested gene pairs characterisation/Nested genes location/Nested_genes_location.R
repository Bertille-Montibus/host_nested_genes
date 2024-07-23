# Characterisation of the nested gene location inside their host #

library(stringr)
library(reshape2)
library("GenomicFeatures")
library("ChIPseeker")
library("ggupset")
library(dplyr)
library(viridis)

##################
#### HUMAN #######
##################

# ChIPseeker was used to annotate the location of the nested genes inside the host genes
# To do this it is necessary to make a gtf file with the host transcripts only

#Open the GTF file
GTF <- read.delim("gencode.v36lift37.annotation.gtf.gz", header=FALSE, comment.char="#")
#Open the list of host and nested genes
hn_gencode <- read.delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt", row.names = 1)

#Make a gtf file for the transcripts of the host genes
#Get the transcript ID separated from the the rest of the last column
GTF2 = GTF
GTF2$transcript = str_extract_all(GTF2$V9, "\\;[^()]+\\;")
#Remove the annotation "transcrip_id" and the ";" to recover only the ensembl transcript ID in the last column
GTF2$transcript = substring(GTF2$transcript, 17, nchar(GTF2$transcript)-1)
GTF2$transcript = sub('\\;.*', '', GTF2$transcript)

#Keep only the rows corresponding of transcripts of host genes
Keep1 = GTF2$transcript %in% hn_gencode$ENST_Host
GTF2_host = GTF2[Keep1,]
GTF_host = GTF[rownames(GTF) %in% rownames(GTF2_host),]

#Removing the host genes which are also nested in another gene to make the annotation easier at first
Unique = setdiff(hn_gencode$SYMBOL_Host, hn_gencode$SYMBOL_Nested)
Unique_ensembl = hn_gencode[hn_gencode$SYMBOL_Host %in% Unique,"ENST_Host"]
Keep2 = GTF2$transcript %in% Unique_ensembl
GTF2_host_only = GTF2[Keep2,]
GTF_host_only = GTF[rownames(GTF) %in% rownames(GTF2_host_only),]

#Save the new GTF files
write.table(GTF_host, file = "wgEncodeGencodeCompV36lift37_compHg19_host.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_host_only, file = "wgEncodeGencodeCompV36lift37_compHg19_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#For the annotation, we also need a bed file for the nested genes
#Make the bed file with all the nested genes
bed = hn_gencode[,c(1,2,3,4,6)]
colnames(bed)=NULL

#Save the table
write.table(bed, file = "wgEncodeGencodeCompV36lift37_compHg19_nested.bed", sep = "\t", row.names = F, quote = F)

#Make a bed file by removing the nested genes which are also host genes to simplify the annotation
Unique = setdiff(hn_gencode$SYMBOL_Nested, hn_gencode$SYMBOL_Host)
bed_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% Unique,c(1,2,3,4,6)]
colnames(bed_nestedonly)=NULL

#Save the table
write.table(bed, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedonly.bed", sep = "\t", row.names = F, quote = F)

#Annotation using ChIPseeker
#Make a txdb object from the personalised GTF file just generated
#Do the annotation with the host only GTF and the all nested bed and the all host gtf and the nested only bed to prevent
#misannotation and merge both at the end
#Start with host only and all nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nested.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_locationinhostonly = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                         TxDb=txdb,
                                         genomicAnnotationPriority = c("Exon","Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))

#Recover the results from the Annotation
Nested_location_table = as.data.frame(Nested_locationinhostonly)
Nested_location_table2 = Nested_locationinhostonly@detailGenomicAnnotation
Nested_location_table = cbind(Nested_location_table, Nested_location_table2)

#Second annotation with all host and nested genes only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_host.gtf.gz",
                       format="gtf",
                       dataSource= "GencodeCompV36lift37",
                       organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedonly.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                   TxDb=txdb,
                                   genomicAnnotationPriority = c( "Exon","Intron", "Promoter","5UTR","3UTR","Downstream", "Intergenic"))

#Recover the results from the Annotation
Nestedonly_location_table = as.data.frame(Nestedonly_location)
Nestedonly_location_table2 = Nestedonly_location@detailGenomicAnnotation
Nestedonly_location_table = cbind(Nestedonly_location_table, Nestedonly_location_table2)

#Pool the data from both annotations
Annotation = rbind(Nested_location_table, Nestedonly_location_table)
#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Few rows are duplicated this is due to the genes which are both host and nested but the annotation is the same
hn_gencode_position = hn_gencode_position[!duplicated(hn_gencode_position[,c(1,2)]),]

#Not all the pairs were annotated following this first round of annotation due to genes involved in multiple pairs and genes which are host but also nested
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation (n=1443)
reanno = hn_gencode_position[is.na(hn_gencode_position$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#2nd round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_3 = rbind(hn_gencode_position,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_3$exon_intron_total_number = as.numeric(hn_gencode_position_3$exon_intron_total_number)
hn_gencode_position_3 <- hn_gencode_position_3[order(hn_gencode_position_3$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_3 = hn_gencode_position_3[!duplicated(hn_gencode_position_3$pairid),]

#Still pairs left without annotation so aplly the reannotation for another round (n=377)
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation
reanno = hn_gencode_position_3[is.na(hn_gencode_position_3$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#3rd round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_4 = rbind(hn_gencode_position_3,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_4$exon_intron_total_number = as.numeric(hn_gencode_position_4$exon_intron_total_number)
hn_gencode_position_4 <- hn_gencode_position_4[order(hn_gencode_position_4$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_4 = hn_gencode_position_4[!duplicated(hn_gencode_position_4$pairid),]

#Still pairs left without annotation so aplly the reannotation for another round (n=243)
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation
reanno = hn_gencode_position_4[is.na(hn_gencode_position_4$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#4th round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_5 = rbind(hn_gencode_position_4,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_5$exon_intron_total_number = as.numeric(hn_gencode_position_5$exon_intron_total_number)
hn_gencode_position_5 <- hn_gencode_position_5[order(hn_gencode_position_5$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_5 = hn_gencode_position_5[!duplicated(hn_gencode_position_5$pairid),]

#Still pairs left without annotation so aplly the reannotation for another round (n=209)
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation
reanno = hn_gencode_position_5[is.na(hn_gencode_position_5$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#5th round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_6 = rbind(hn_gencode_position_5,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_6$exon_intron_total_number = as.numeric(hn_gencode_position_6$exon_intron_total_number)
hn_gencode_position_6 <- hn_gencode_position_6[order(hn_gencode_position_6$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_6 = hn_gencode_position_6[!duplicated(hn_gencode_position_6$pairid),]

#Still pairs left without annotation so aplly the reannotation for another round (n=206)
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation
reanno = hn_gencode_position_6[is.na(hn_gencode_position_6$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#6th round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "wgEncodeGencodeCompV36lift37_compHg19_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeCompV36lift37",
                        organism="Homo sapiens")
#Make a peak object from the bed file
Nested_genes = readPeakFile("wgEncodeGencodeCompV36lift37_compHg19_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_7 = rbind(hn_gencode_position_6,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_7$exon_intron_total_number = as.numeric(hn_gencode_position_7$exon_intron_total_number)
hn_gencode_position_7 = hn_gencode_position_7[order(hn_gencode_position_7$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_7 = hn_gencode_position_7[!duplicated(hn_gencode_position_7$pairid),]

#No improvement here anymore (still 206 pairs without annotation) --> manual annotation of these pairs left
#Save the table and complete the table manually
write.table(hn_gencode_position_11, file = "hn_gencode_h19_position.txt", col.names = T, sep = "\t")

#After manual annotation of the remaining pair, overlap the annotation table with the initiale host nested gene table
hn_gencode_position = read.delim("hn_gencode_h19_position.txt", header = T, sep = "\t", row.names = 1)

df = merge(hn_gencode, hn_gencode_position, by.x = 17, by.y = 17)
#Remove the column 18 to 39 which are already availble or not useful for the analysis
df = df[,c(-18:-39)]
#Change the column names back to the names in the original table
colnames(df)[2:17] = colnames(hn_gencode)[1:16]
#Name it again hn_gencode_position
hn_gencode_position = df

#Make graph from the table

#Add an annotation for nested genes when they overlap both exonic and intronic sequence as for the moment they are annotated as overlaping exon
for (i in 1:length(rownames(hn_gencode_position))){
  if (hn_gencode_position[i,"Exon"] == TRUE & hn_gencode_position[i,"Intron"] == TRUE) {
    hn_gencode_position[i,"position"] = "Exon/Intron"
  }else {
    hn_gencode_position[i,"position"] = hn_gencode_position[i,"annotation"]
  }
}

#Count the different categories
annotation_count <- hn_gencode_position %>% group_by(position) %>% count() %>% as.data.frame()
#Include information about the orientation
annotation_count_orientation <- hn_gencode_position %>% group_by(position,Strandness) %>% count() %>% as.data.frame()
annotation_count_orientation$Strandness = as.factor(annotation_count_orientation$Strandness)

#Boxplot with the proportion of pairs in each categories
ggplot(annotation_count_orientation, aes(x="", y=n, fill=position, alpha=Strandness))+
  geom_bar(width = 1, position = "stack", stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_viridis(discrete=TRUE,option="C")+
  scale_alpha_manual(values=c(0.5, 1)) +
  geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%", "\n","(", n,")")),
            position = position_stack(vjust = 0.5), size = 8)


#For the fully intronic and exonic look in more details if there is a bias according toward the position of the intron or exon containing the nested gene
#For the intronic one
#Intron
hn_gencode_position_intron= hn_gencode_position[hn_gencode_position$position == "Intron",]
for (i in 1:length(rownames(hn_gencode_position_intron))){
  if (hn_gencode_position_intron[i,"exon_intron_number"] == 1 & hn_gencode_position_intron[i,"exon_intron_total_number"] == 1) {
    hn_gencode_position_intron[i,"Intron_number"] = "First and only intron"
  }else {
    if (hn_gencode_position_intron[i,"exon_intron_number"] == 1 & hn_gencode_position_intron[i,"exon_intron_total_number"] != 1) {
      hn_gencode_position_intron[i,"Intron_number"] = "First intron"
    } else {
      if (hn_gencode_position_intron[i,"exon_intron_number"] == hn_gencode_position_intron[i,"exon_intron_total_number"]) {
        hn_gencode_position_intron[i,"Intron_number"] = "Last intron"
      }else {
        hn_gencode_position_intron[i,"Intron_number"] = "Other intron"
      }
    }
  }
}

annotation_count <- hn_gencode_position_intron %>% group_by(Intron_number,Strandness) %>% count() %>% as.data.frame()

ggplot(data=annotation_count, aes(x=Intron_number, y=n,alpha=Strandness)) +
  geom_bar(stat="identity")+
  scale_alpha_manual(values=c(0.5, 1))

#Exon
hn_gencode_position_exon= hn_gencode_position[hn_gencode_position$position == "Exon",]
for (i in 1:length(rownames(hn_gencode_position_exon))){
  if (hn_gencode_position_exon[i,"exon_intron_number"] == 1 & hn_gencode_position_exon[i,"exon_intron_total_number"] == 1) {
    hn_gencode_position_exon[i,"Exon_number"] = "First and only exon"
  }else {
    if (hn_gencode_position_exon[i,"exon_intron_number"] == 1 & hn_gencode_position_exon[i,"exon_intron_total_number"] != 1) {
      hn_gencode_position_exon[i,"Exon_number"] = "First exon"
    } else {
      if (hn_gencode_position_exon[i,"exon_intron_number"] == hn_gencode_position_exon[i,"exon_intron_total_number"]) {
        hn_gencode_position_exon[i,"Exon_number"] = "Last exon"
      }else {
        hn_gencode_position_exon[i,"Exon_number"] = "Other exon"
      }
    }
  }
}

annotation_count <- hn_gencode_position_exon %>% group_by(Exon_number, Strandness) %>% count() %>% as.data.frame()

ggplot(data=annotation_count, aes(x=Exon_number, y=n, alpha=Strandness)) +
  geom_bar(stat="identity") +
  scale_alpha_manual(values=c(0.5, 1))

#Graph more detailed to include more information for intronic and exonic location
#What is the distribution of the
count_intron_number = hn_gencode_position_intron %>% group_by(exon_intron_total_number) %>% count() %>% as.data.frame()
#For visualisation purpose selection of the gene with less than 20 introns (84.53% of the genes)
count_intron_number_less20 = count_intron_number[count_intron_number$exon_intron_total_number <= 20, ]
sum(count_intron_number_less20$n)/sum(count_intron_number$n)*100

#Table to make the plot
count_intron_number_location = hn_gencode_position_intron %>% group_by(exon_intron_number,exon_intron_total_number) %>% count() %>% as.data.frame()
#Select only the genes between 2 and 20 introns
count_intron_number_location_graph1 = count_intron_number_location[count_intron_number_location$exon_intron_total_number >1 & count_intron_number_location$exon_intron_total_number <= 20, ]
#Prepare the data for the graph
count_intron_number_location_graph1$exon_intron_number = as.factor(count_intron_number_location_graph1$exon_intron_number)
count_intron_number_location_graph1$exon_intron_total_number = as.factor(count_intron_number_location_graph1$exon_intron_total_number)

ggplot(count_intron_number_location_graph1, aes(fill=exon_intron_number, y=n, x=exon_intron_total_number)) +
  geom_bar(position="dodge", stat="identity")+
  scale_fill_viridis(discrete = T, option = "C", direction = -1)

#For the nested genes fully contained into an intron, recover the size of the intron to see if they are longer
#Keep only the information about exons
GTF3 = GTF[GTF$V3 == "exon",]

#Transform the last column into individual columns
GTF3 = cbind(GTF3,(colsplit(GTF3[,9], ";",
                            names=c("GeneId","TranscriptID",
                                    "Gene_type", "Gene_name",
                                    "Transcript_type", "Transcript_name",
                                    "Exon_number","Rest"))))

#Remove gene_id and transcript_id from the column to keep only the ensembl ID and be able to select them
GTF3$GeneId = gsub("gene_id ", "", GTF3$GeneId)
GTF3$TranscriptID = gsub(" transcript_id ", "", GTF3$TranscriptID)

#Adjust for some of the transcript which in addition to the normal column have gene_status and transcript_status
#Isolate them
Status = GTF3[(grepl("gene_status",GTF3$Gene_name) | grepl("transcript_status",GTF3$Exon_number)), ]
#Remove the additional columns
Status = Status[,c(-13,-16)]
#Reisolate the information from the last column
GTF3_status = cbind(Status,colsplit(Status[,15], ";",
                                    names=c("Transcript_name","Exon_number","Rest")))

GTF3_status = GTF3_status[,c(-15)]
#Change the column names to adjust to the new information
colnames(GTF3_status) = colnames(GTF3)

#Remove these from the GTF file to put them back
GTF3 = GTF3[!(grepl("gene_status",GTF3$Gene_name) | grepl("transcript_status",GTF3$Exon_number)), ]
GTF3 = rbind(GTF3, GTF3_status)

#Remove the transcripts which are duplicated because on the Y and X chr because it prevent removing the last intron which is not an intron
remove = grepl("PAR_Y",GTF3$TranscriptID)
GTF3 = GTF3[!remove,]

#Remove gene_name from the column to keep only the name and be able to select them
GTF3$Gene_name = gsub(" gene_name ", "", GTF3$Gene_name)
#Remove Gene_type, Transcript_type, transcript_name and Rest
GTF3 = GTF3[,c(-12,-14,-15,-17)]
#Remove exon_number from the column to keep only the number and be able to select them
GTF3$Exon_number = gsub(" exon_number ", "", GTF3$Exon_number)
#Remove the initial big column with all the information
GTF3 = GTF3[,-9]
GTF3$Exon_number = as.numeric(GTF3$Exon_number)
#Calculate the size of the exon
GTF3$Exon_size = GTF3$V5-GTF3$V4

#Comparison with introns into hosts and all the other genes 
#Arrange the transcripts by name and by exon number
GTF4 = GTF3
GTF4=GTF4 %>% arrange(TranscriptID, Exon_number)

#Calculate the intron size
#The calculation is different according to the orientation of the transcript
for (i in 1:length(rownames(GTF4))) {
  if (GTF4[i,"V7"] == "+") {
    GTF4[i,"Intron_size"] = GTF4[i+1,"V4"]-GTF4[i,"V5"]
  } else {
    GTF4[i,"Intron_size"] = GTF4[i,"V4"]-GTF4[i+1,"V5"]
  }
  print(i)
}

#Because all the genes are one after each other we need to remove the row which correspond to the last exon as there is no intron after
#Calculate the total number of introns
Exon_number = GTF4 %>% count(TranscriptID)
colnames(Exon_number)[2]="Exon_count"

#Merge this with the table where we calculated the intron size
GTF5 = merge(GTF4,Exon_number, by.x = 10, by.y = 1, all.x = T)
#Remove the rows where exon number is equal to the total number of exons
Remove = GTF5$Exon_number == GTF5$Exon_count
GTF6 = GTF5[!Remove,]

#Make the table with intron size and save it
All_Intron_size = GTF6[,c(1,11,12,14)]
colnames(All_Intron_size)[3]= "Intron_number"

#Save the table 
write.table(All_Intron_size, file = "gencode_v36lift37_annotation_intron_size.txt", sep = "\t", col.names = T)

#Open the table with all the transcripts use at the beginning
all <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")
#Select only the introns form these transcripts as a reference
All_Intron_size = All_Intron_size[All_Intron_size$TranscriptID %in% all$V4,]
All_Intron_size$legend = "All_genes"

#Recover the size of the nested genes containing introns by merging the size table with the intron table
#Merge for the size
hn_gencode_position_intron_size = merge.data.frame(hn_gencode_position_intron,All_Intron_size , by.x = c(12,28), by.y = c(1,3), all.x = T)
hn_gencode_position_intron_size$legend = "Nested_containing_intron"

#Recover all the introns of the host genes to compare
hn_gencode_host_intron_size = merge.data.frame(hn_gencode_position_intron,All_Intron_size , by.x = c(12), by.y = c(1), all.x = T)
hn_gencode_host_intron_size$legend = "Host_genes_intron"

#Merge the tables for the graph
graph = rbind(hn_gencode_position_intron_size[,33:34],All_Intron_size[,4:5],hn_gencode_host_intron_size[,34:35])

#Make the boxplot
ggplot(graph, aes(x=legend, y= Intron_size)) +
  geom_boxplot(fill='#A4A4A4', color="black") +
  scale_y_continuous(trans='log10')+
  labs(title="Comparaison with all introns",x="", y = "Intron size")

#Stats
##Human
t.test(log10(All_Intron_size$Intron_size),log10(hn_gencode_host_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(All_Intron_size$Intron_size) and log10(hn_gencode_host_intron_size$Intron_size)
#t = -122.98, df = 128884, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3151297 -0.3052427
#sample estimates:
#  mean of x mean of y
# 3.229999  3.540186

t.test(log10(All_Intron_size$Intron_size),log10(hn_gencode_position_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(All_Intron_size$Intron_size) and log10(hn_gencode_position_intron_size$Intron_size)
#t = -185.11, df = 10278, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -1.295433 -1.268285
#sample estimates:
#  mean of x mean of y
#  3.229999  4.511859


t.test(log10(hn_gencode_host_intron_size$Intron_size),log10(hn_gencode_position_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(hn_gencode_host_intron_size$Intron_size) and log10(hn_gencode_position_intron_size$Intron_size)
#t = -132.96, df = 12730, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#   -0.9859977 -0.9573483
#sample estimates:
#  mean of x mean of y
#3.540186  4.511859 


##################
#### MOUSE #######
##################
# ChIPseeker was used to annotate the location of the nested genes inside the host genes
# To do this it is necessary to make a gtf file with the host transcripts only

#Open the GTF file
GTF <- read.delim("gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz", header=FALSE, comment.char="#")
#Open the list of host and nested genes
hn_gencode <- read.delim("Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", sep = "\t")

#Make a gtf file for the transcripts of the host genes
#Get the transcript ID separated from the the rest of the last column
GTF2 = GTF
GTF2$transcript = str_extract_all(GTF2$V9, "\\;[^()]+\\;")
#Remove the annotation "transcrip_id" and the ";" to recover only the ensembl transcript ID in the last column
GTF2$transcript = substring(GTF2$transcript, 17, nchar(GTF2$transcript)-1)
GTF2$transcript = sub('\\;.*', '', GTF2$transcript)

#Change the name of the random chromosomes for later as this will be needed
GTF$V1 = gsub("GL456221.1", "chr1_GL456221_random", GTF$V1)
GTF$V1 = gsub("GL456233.1", "chrX_GL456233_random", GTF$V1)

#Keep only the rows corresponding of transcripts of host genes
Keep1 = GTF2$transcript %in% hn_gencode$ENST_Host
GTF2_host = GTF2[Keep1,]
GTF_host = GTF[rownames(GTF) %in% rownames(GTF2_host),]

#Removing the host genes which are also nested in another gene to make the annotation easier at first
Unique = setdiff(hn_gencode$SYMBOL_Host, hn_gencode$SYMBOL_Nested)
Unique_ensembl = hn_gencode[hn_gencode$SYMBOL_Host %in% Unique,"ENST_Host"]
Keep2 = GTF2$transcript %in% Unique_ensembl
GTF2_host_only = GTF2[Keep2,]
GTF_host_only = GTF[rownames(GTF) %in% rownames(GTF2_host_only),]

#Save the new GTF files
write.table(GTF_host, file = "GencodeVM25_comprehensive_mm10_host.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_host_only, file = "GencodeVM25_comprehensive_mm10_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)

#For the annotation, we also need a bed file for the nested genes
#Make the bed file with all the nested genes
bed = hn_gencode[,c(1,2,3,4,6)]
colnames(bed)=NULL

#Save the table
write.table(bed, file = "GencodeVM25_comprehensive_mm10_nested.bed", sep = "\t", row.names = F, quote = F)

#Make a bed file by removing the nested genes which are also host genes to simplify the annotation
Unique = setdiff(hn_gencode$SYMBOL_Nested, hn_gencode$SYMBOL_Host)
bed_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% Unique,c(1,2,3,4,6)]
colnames(bed_nestedonly)=NULL

#Save the table
write.table(bed, file = "GencodeVM25_comprehensive_mm10_nestedonly.bed", sep = "\t", row.names = F, quote = F)

#Annotation using ChIPseeker
#Make a txdb object from the personalised GTF file just generated
#Do the annotation with the host only GTF and the all nested bed and the all host gtf and the nested only bed to prevent
#misannotation and merge both at the end
#Start with host only and all nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes = readPeakFile("GencodeVM25_comprehensive_mm10_nested.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_locationinhostonly = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                         TxDb=txdb,
                                         genomicAnnotationPriority = c("Exon","Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))

#Recover the results from the Annotation
Nested_location_table = as.data.frame(Nested_locationinhostonly)
Nested_location_table2 = Nested_locationinhostonly@detailGenomicAnnotation
Nested_location_table = cbind(Nested_location_table, Nested_location_table2)

#Second annotation with all host and nested genes only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_host.gtf.gz",
                       format="gtf",
                       dataSource= "GencodeVM25",
                       organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("GencodeVM25_comprehensive_mm10_nestedonly.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                   TxDb=txdb,
                                   genomicAnnotationPriority = c( "Exon","Intron", "Promoter","5UTR","3UTR","Downstream", "Intergenic"))

#Recover the results from the Annotation
Nestedonly_location_table = as.data.frame(Nestedonly_location)
Nestedonly_location_table2 = Nestedonly_location@detailGenomicAnnotation
Nestedonly_location_table = cbind(Nestedonly_location_table, Nestedonly_location_table2)

#Pool the data from both annotations
Annotation = rbind(Nested_location_table, Nestedonly_location_table)
#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Few rows are duplicated this is due to the genes which are both host and nested but the annotation is the same
hn_gencode_position = hn_gencode_position[!duplicated(hn_gencode_position[,c(1,2)]),]

#Not all the pairs were annotated following this first round of annotation due to genes involved in multiple pairs and genes which are host but also nested
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation (n=455)
reanno = hn_gencode_position[is.na(hn_gencode_position$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "GencodeVM25_comprehensive_mm10_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)
#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#2nd round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_3 = rbind(hn_gencode_position,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_3$exon_intron_total_number = as.numeric(hn_gencode_position_3$exon_intron_total_number)
hn_gencode_position_3 <- hn_gencode_position_3[order(hn_gencode_position_3$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_3 = hn_gencode_position_3[!duplicated(hn_gencode_position_3$pairid),]

#Not all the pairs were annotated following this first round of annotation due to genes involved in multiple pairs and genes which are host but also nested
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation (n=88)
reanno = hn_gencode_position_3[is.na(hn_gencode_position_3$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "GencodeVM25_comprehensive_mm10_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)
#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#3rd round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_4 = rbind(hn_gencode_position_3,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_4$exon_intron_total_number = as.numeric(hn_gencode_position_4$exon_intron_total_number)
hn_gencode_position_4 <- hn_gencode_position_4[order(hn_gencode_position_4$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_4 = hn_gencode_position_4[!duplicated(hn_gencode_position_4$pairid),]

#Not all the pairs were annotated following this first round of annotation due to genes involved in multiple pairs and genes which are host but also nested
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation (n=50)
reanno = hn_gencode_position_4[is.na(hn_gencode_position_4$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "GencodeVM25_comprehensive_mm10_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)
#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#4th round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_5 = rbind(hn_gencode_position_4,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_5$exon_intron_total_number = as.numeric(hn_gencode_position_5$exon_intron_total_number)
hn_gencode_position_5 <- hn_gencode_position_5[order(hn_gencode_position_5$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_5 = hn_gencode_position_5[!duplicated(hn_gencode_position_5$pairid),]

#Not all the pairs were annotated following this first round of annotation due to genes involved in multiple pairs and genes which are host but also nested
#Make files with the remaining genes to make a new round of annotation
#Identified the pairs without annotation (n=48)
reanno = hn_gencode_position_5[is.na(hn_gencode_position_5$annotation),"pairid"]
reanno_table = hn_gencode[hn_gencode$pairid %in% reanno,]
#Select the nested only and the host only
reanno_nested_only = reanno_table[!(reanno_table$SYMBOL_Nested %in% reanno_table$SYMBOL_Host),]
reanno_host_only = reanno_table[!(reanno_table$SYMBOL_Host %in% reanno_table$SYMBOL_Nested),]

#Make a new GTF file for the reannotation with the remaining host or host only
Keep1 = GTF2$transcript %in% reanno_table$ENST_Host
Keep2 = GTF2$transcript %in% reanno_host_only$ENST_Host
GTF_reanno = GTF[Keep1,]
GTF_reanno_hostonly = GTF[Keep2,]

#Save the GTF for the next round of annotation
write.table(GTF_reanno, file = "GencodeVM25_comprehensive_mm10_reanno.gtf.gz", sep = "\t", row.names = F, quote = F)
write.table(GTF_reanno_hostonly, file = "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz", sep = "\t", row.names = F, quote = F)
#Make the bed file for the reannotation
bed_nestedreanno = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_table$SYMBOL_Nested,c(1,2,3,4,6)]
bed_nestedreanno_nestedonly = hn_gencode[hn_gencode$SYMBOL_Nested %in% reanno_nested_only$SYMBOL_Nested,c(1,2,3,4,6)]
colnames(bed_nestedreanno)=NULL
colnames(bed_nestedreanno_nestedonly)=NULL

#Save the bed files
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno.bed", sep = "\t", row.names = F, quote = F)
write.table(bed_nestedreanno, file = "GencodeVM25_comprehensive_mm10_nestedreanno_only.bed", sep = "\t", row.names = F, quote = F)

#5th round of annotation with ChIPseeker with the all the remaining host and nested only
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes_only = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno_only.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nestedonly_location_reanno = annotatePeak(Nested_genes_only, tssRegion=c(-100, 100),
                                          TxDb=txdb,
                                          genomicAnnotationPriority = c("Exon", "Intron","Promoter", "5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nestedonly_reanno = as.data.frame(Nestedonly_location_reanno)
Nestedonly_reanno_table2 = Nestedonly_location_reanno@detailGenomicAnnotation
Nestedonly_reanno = cbind(Nestedonly_reanno, Nestedonly_reanno_table2)

#ChIPseeker with the host only remaining and all the nested genes
#Make a txbd object from the gtf
txdb = makeTxDbFromGFF(file= "GencodeVM25_comprehensive_mm10_reanno_hostonly.gtf.gz",
                        format="gtf",
                        dataSource= "GencodeVM25",
                        organism="Mus musculus")
#Make a peak object from the bed file
Nested_genes = readPeakFile("GencodeVM25_comprehensive_mm10_nestedreanno.bed")

#Do the annotation by changing the priorities. We put Exon and Intron first so they should all be annotated as they are inside another gene
Nested_location_reanno = annotatePeak(Nested_genes, tssRegion=c(-100, 100),
                                      TxDb=txdb,
                                      genomicAnnotationPriority = c("Exon", "Intron","Promoter","5UTR","3UTR","Downstream", "Intergenic"))
#Recover the results
Nested_reanno = as.data.frame(Nested_location_reanno)
Nested_reanno_table2 = Nested_location_reanno@detailGenomicAnnotation
Nested_reanno = cbind(Nested_reanno, Nested_reanno_table2)

#Pool the data for this round of annotation
Annotation = rbind(Nested_reanno,Nestedonly_reanno)

#Extraction of the gene information from the results to get the exact location inside the gene
Annotation$Host = str_extract_all(Annotation$annotation,"(?<=\\().+?(?=\\))")
Annotation$annotation = gsub("\\s*\\([^\\)]+\\)","",Annotation$annotation)
Annotation$exon_intron_number = sub(".*,","",Annotation$Host)
Annotation$Host = gsub("/", "-", Annotation$Host)
Annotation$Host = sub("-.*", "", Annotation$Host)
Annotation$exon_intron_total_number = sub(".*of","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("\\of.*","",Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("intron", "", Annotation$exon_intron_number)
Annotation$exon_intron_number = gsub("exon", "", Annotation$exon_intron_number)

#Add the non-host gene name for the non intronic/exonic annotations
for (i in 1:nrow(Annotation)) {
  if (!Annotation[i,8] %in% c("Intron","Exon")) {
    Annotation[i,26] = Annotation[i,15]
  }
}

#Remove the information of the closest gene given by default by ChIPseeker
Annotation = Annotation[,-9:-16]
Annotation = Annotation[!duplicated(Annotation),]

#Add the information to the host/nested genes table by merging using the host and nested transcript IDs
hn_gencode_position_2 = merge.data.frame(hn_gencode,Annotation, by.x =c(4,11), by.y = c(6,18),all.x = T)
#Merge with the annotation
hn_gencode_position_6 = rbind(hn_gencode_position_5,hn_gencode_position_2)
#Remove the duplicated rows for pairids but make sure that the one removed are not with NA so used the exon_intron_total_number for that
hn_gencode_position_6$exon_intron_total_number = as.numeric(hn_gencode_position_6$exon_intron_total_number)
hn_gencode_position_6 <- hn_gencode_position_6[order(hn_gencode_position_6$exon_intron_total_number,decreasing=TRUE, na.last=TRUE),]
hn_gencode_position_6 = hn_gencode_position_6[!duplicated(hn_gencode_position_6$pairid),]

#No improvement here anymore (still 48 pairs without annotation) --> manual annotation of these pairs left
write.table(hn_gencode_position_6, file = "hn_gencode_mm10_position.txt", col.names = T, sep = "\t")


#After manual annotation of the remaining pair, overlap the annotation table with the initial host nested gene table
hn_gencode_position = read.delim("hn_gencode_mm10_position.txt", header = T, sep = "\t", row.names = 1)

df = merge(hn_gencode, hn_gencode_position, by.x = 17, by.y = 17)
#Remove the column 18 to 39 which are already availble or not useful for the analysis
df = df[,c(-18:-39)]
#Change the column names back to the names in the original table
colnames(df)[2:17] = colnames(hn_gencode)[1:16]
#Name it again hn_gencode_position
hn_gencode_position = df

#Make graph from the table

#Add an annotation for nested genes when they overlap both exonic and intronic sequence as for the moment they are annotated as overlaping exon
for (i in 1:length(rownames(hn_gencode_position))){
  if (hn_gencode_position[i,"Exon"] == TRUE & hn_gencode_position[i,"Intron"] == TRUE) {
    hn_gencode_position[i,"position"] = "Exon/Intron"
  }else {
    hn_gencode_position[i,"position"] = hn_gencode_position[i,"annotation"]
  }
}

#Count the different categories
annotation_count <- hn_gencode_position %>% group_by(position) %>% count() %>% as.data.frame()
#Include information about the orientation
annotation_count_orientation <- hn_gencode_position %>% group_by(position,Strandness) %>% count() %>% as.data.frame()
annotation_count_orientation$Strandness = as.factor(annotation_count_orientation$Strandness)

#Boxplot with the proportion of pairs in each categories
ggplot(annotation_count_orientation, aes(x="", y=n, fill=position, alpha=Strandness))+
  geom_bar(width = 1, position = "stack", stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_viridis(discrete=TRUE,option="C")+
  scale_alpha_manual(values=c(0.5, 1)) +
  geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%", "\n","(", n,")")),
            position = position_stack(vjust = 0.5), size = 8)


#For the fully intronic and exonic look in more details if there is a bias according toward the position of the intron or exon containing the nested gene
#For the intronic one
#Intron
hn_gencode_position_intron= hn_gencode_position[hn_gencode_position$position == "Intron",]
for (i in 1:length(rownames(hn_gencode_position_intron))){
  if (hn_gencode_position_intron[i,"exon_intron_number"] == 1 & hn_gencode_position_intron[i,"exon_intron_total_number"] == 1) {
    hn_gencode_position_intron[i,"Intron_number"] = "First and only intron"
  }else {
    if (hn_gencode_position_intron[i,"exon_intron_number"] == 1 & hn_gencode_position_intron[i,"exon_intron_total_number"] != 1) {
      hn_gencode_position_intron[i,"Intron_number"] = "First intron"
    } else {
      if (hn_gencode_position_intron[i,"exon_intron_number"] == hn_gencode_position_intron[i,"exon_intron_total_number"]) {
        hn_gencode_position_intron[i,"Intron_number"] = "Last intron"
      }else {
        hn_gencode_position_intron[i,"Intron_number"] = "Other intron"
      }
    }
  }
}

annotation_count <- hn_gencode_position_intron %>% group_by(Intron_number,Strandness) %>% count() %>% as.data.frame()

ggplot(data=annotation_count, aes(x=Intron_number, y=n,alpha=Strandness)) +
  geom_bar(stat="identity")+
  scale_alpha_manual(values=c(0.5, 1))

#Exon
hn_gencode_position_exon= hn_gencode_position[hn_gencode_position$position == "Exon",]
for (i in 1:length(rownames(hn_gencode_position_exon))){
  if (hn_gencode_position_exon[i,"exon_intron_number"] == 1 & hn_gencode_position_exon[i,"exon_intron_total_number"] == 1) {
    hn_gencode_position_exon[i,"Exon_number"] = "First and only exon"
  }else {
    if (hn_gencode_position_exon[i,"exon_intron_number"] == 1 & hn_gencode_position_exon[i,"exon_intron_total_number"] != 1) {
      hn_gencode_position_exon[i,"Exon_number"] = "First exon"
    } else {
      if (hn_gencode_position_exon[i,"exon_intron_number"] == hn_gencode_position_exon[i,"exon_intron_total_number"]) {
        hn_gencode_position_exon[i,"Exon_number"] = "Last exon"
      }else {
        hn_gencode_position_exon[i,"Exon_number"] = "Other exon"
      }
    }
  }
}

annotation_count <- hn_gencode_position_exon %>% group_by(Exon_number, Strandness) %>% count() %>% as.data.frame()

ggplot(data=annotation_count, aes(x=Exon_number, y=n, alpha=Strandness)) +
  geom_bar(stat="identity") +
  scale_alpha_manual(values=c(0.5, 1))

#Graph more detailed to include more information for intronic and exonic location
#What is the distribution of the
count_intron_number = hn_gencode_position_intron %>% group_by(exon_intron_total_number) %>% count() %>% as.data.frame()
#For visualisation purpose selection of the gene with less than 20 introns (81.84% of the genes)
count_intron_number_less20 = count_intron_number[count_intron_number$exon_intron_total_number <= 20, ]
sum(count_intron_number_less20$n)/sum(count_intron_number$n)*100

#Table to make the plot
count_intron_number_location = hn_gencode_position_intron %>% group_by(exon_intron_number,exon_intron_total_number) %>% count() %>% as.data.frame()
#Select only the genes between 2 and 20 introns
count_intron_number_location_graph1 = count_intron_number_location[count_intron_number_location$exon_intron_total_number >1 & count_intron_number_location$exon_intron_total_number <= 20, ]
#Prepare the data for the graph
count_intron_number_location_graph1$exon_intron_number = as.factor(count_intron_number_location_graph1$exon_intron_number)
count_intron_number_location_graph1$exon_intron_total_number = as.factor(count_intron_number_location_graph1$exon_intron_total_number)

ggplot(count_intron_number_location_graph1, aes(fill=exon_intron_number, y=n, x=exon_intron_total_number)) +
  geom_bar(position="dodge", stat="identity")+
  scale_fill_viridis(discrete = T, option = "C", direction = -1)


#For the nested genes fully contained into an intron, recover the size of the intron to see if they are longer
#Keep only the information about exons
GTF3 = GTF[GTF$V3 == "exon",]

#Transform the last column into individual columns
GTF3 = cbind(GTF3,(colsplit(GTF3[,9], ";",
                            names=c("GeneId","TranscriptID",
                                    "Gene_type", "Gene_name",
                                    "Transcript_type", "Transcript_name",
                                    "Exon_number","Rest"))))

#Remove gene_id and transcript_id from the column to keep only the ensembl ID and be able to select them
GTF3$GeneId = gsub("gene_id ", "", GTF3$GeneId)
GTF3$TranscriptID = gsub(" transcript_id ", "", GTF3$TranscriptID)

#Remove the transcripts which are duplicated because on the Y and X chr because it prevent removing the last intron which is not an intron
remove = grepl("PAR_Y",GTF3$TranscriptID)
GTF3 = GTF3[!remove,]

#Remove gene_name from the column to keep only the name and be able to select them
GTF3$Gene_name = gsub(" gene_name ", "", GTF3$Gene_name)
#Remove Gene_type, Transcript_type, transcript_name and Rest
GTF3 = GTF3[,c(-12,-14,-15,-17)]
#Remove exon_number from the column to keep only the number and be able to select them
GTF3$Exon_number = gsub(" exon_number ", "", GTF3$Exon_number)
#Remove the initial big column with all the information
GTF3 = GTF3[,-9]
GTF3$Exon_number = as.numeric(GTF3$Exon_number)
#Calculate the size of the exon
GTF3$Exon_size = GTF3$V5-GTF3$V4

#Comparison with introns into hosts and all the other genes 
#Arrange the transcripts by name and by exon number
GTF4 = GTF3
GTF4=GTF4 %>% arrange(TranscriptID, Exon_number)

#Calculate the intron size
#The calculation is different according to the orientation of the transcript
for (i in 1:length(rownames(GTF4))) {
  if (GTF4[i,"V7"] == "+") {
    GTF4[i,"Intron_size"] = GTF4[i+1,"V4"]-GTF4[i,"V5"]
  } else {
    GTF4[i,"Intron_size"] = GTF4[i,"V4"]-GTF4[i+1,"V5"]
  }
  print(i)
}

#Because all the genes are one after each other we need to remove the row which correspond to the last exon as there is no intron after
#Calculate the total number of introns
Exon_number = GTF4 %>% count(TranscriptID)
colnames(Exon_number)[2]="Exon_count"

#Merge this with the table where we calculated the intron size
GTF5 = merge(GTF4,Exon_number, by.x = 10, by.y = 1, all.x = T)
#Remove the rows where exon number is equal to the total number of exons
Remove = GTF5$Exon_number == GTF5$Exon_count
GTF6 = GTF5[!Remove,]

#Make the table with intron size and save it
All_Intron_size = GTF6[,c(1,11,12,14)]
colnames(All_Intron_size)[3]= "Intron_number"

#Save the table 
write.table(All_Intron_size, file = "gencode_vM25_annotation_intron_size.txt", sep = "\t", col.names = T)

#Open the table with all the transcripts use at the beginning
all <- read.table("wgEncodeGencodeCompVM25_fileforhost_filtered.bed")
#Select only the introns form these transcripts as a reference
All_Intron_size = All_Intron_size[All_Intron_size$TranscriptID %in% all$V4,]
All_Intron_size$legend = "All_genes"

#Recover the size of the nested genes containing introns by merging the size table with the intron table
#Merge for the size
hn_gencode_position_intron_size = merge.data.frame(hn_gencode_position_intron,All_Intron_size , by.x = c(12,28), by.y = c(1,3), all.x = T)
hn_gencode_position_intron_size$legend = "Nested_containing_intron"

#Recover all the introns of the host genes to compare
hn_gencode_host_intron_size = merge.data.frame(hn_gencode_position_intron,All_Intron_size , by.x = c(12), by.y = c(1), all.x = T)
hn_gencode_host_intron_size$legend = "Host_genes_intron"

#Merge the tables for the graph
graph = rbind(hn_gencode_position_intron_size[,33:34],All_Intron_size[,4:5],hn_gencode_host_intron_size[,34:35])

#Make the boxplot
ggplot(graph, aes(x=legend, y= Intron_size)) +
  geom_boxplot(fill='#A4A4A4', color="black") +
  scale_y_continuous(trans='log10')+
  labs(title="Comparaison with all introns",x="", y = "Intron size")

#Stats
t.test(log10(All_Intron_size$Intron_size),log10(hn_gencode_host_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(All_Intron_size$Intron_size) and log10(hn_gencode_host_intron_size$Intron_size)
#t = -90.768, df = 82438, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -0.2849661 -0.2729194
#sample estimates:
#  mean of x mean of y 
#3.132255  3.411198

t.test(log10(All_Intron_size$Intron_size),log10(hn_gencode_position_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(All_Intron_size$Intron_size) and log10(hn_gencode_position_intron_size$Intron_size)
#t = -102.96, df = 5752.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.277569 -1.229827
#sample estimates:
#  mean of x mean of y 
#3.132255  4.385953 

t.test(log10(hn_gencode_host_intron_size$Intron_size),log10(hn_gencode_position_intron_size$Intron_size))
#Welch Two Sample t-test

#data:  log10(hn_gencode_host_intron_size$Intron_size) and log10(hn_gencode_position_intron_size$Intron_size)
#t = -77.957, df = 6393.7, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.9992668 -0.9502436
#sample estimates:
#  mean of x mean of y 
#3.411198  4.385953 
