# Analysis of the conservation of the pairs #

#Download the table of the orthologs between mouse and human using the protocol described here :
#https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/

orthologs = read.csv("mart_export.txt", header = T)
#Open the table with host and nested genes for both species
hn_gencode_human = read.delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt", row.names = 1)
hn_gencode_mouse = read.delim("Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", row.names = 1)

#Find orthologs starting from the HUMAN list
#As stable ID version are needed for this, remove the numbers after the "_" for the transcripts
hn_gencode_human$ENST_Host = gsub("\\_.*","",hn_gencode_human$ENST_Host)
hn_gencode_human$ENST_Nested = gsub("\\_.*","",hn_gencode_human$ENST_Nested)

#Merge the tables with the orthologs tables to find which host and nested genes have an ortholog
hn_gencode_human_orthologs = merge(hn_gencode_human,orthologs[, c(4,6)],by.x = 11,by.y = 1, all.x = T)
colnames(hn_gencode_human_orthologs)[18] = "Host_mouse_ortholog"
hn_gencode_human_orthologs = merge(hn_gencode_human_orthologs,orthologs[, c(4,6)],by.x = 5,by.y = 1, all.x = T)
colnames(hn_gencode_human_orthologs)[19] = "Nested_mouse_ortholog"

#Calculate the number of host and nested genes in human with at least one ortholog in mouse
#Host
hn_gencode_human_orthologs_host = hn_gencode_human_orthologs[!is.na(hn_gencode_human_orthologs$Host_mouse_ortholog),]
length(unique(hn_gencode_human_orthologs_host$ENST_Host))
#5092 host transcripts with ortholog
length(unique(hn_gencode_human_orthologs_host$SYMBOL_Host))
length(unique(hn_gencode_human$SYMBOL_Host))
#4904 host genes with ortholog out of 7661 --> 64%
#Nested
hn_gencode_human_orthologs_nested = hn_gencode_human_orthologs[!is.na(hn_gencode_human_orthologs$Nested_mouse_ortholog),]
length(unique(hn_gencode_human_orthologs_nested$ENST_Nested))
#1835 transcripts nested with ortholog
length(unique(hn_gencode_human_orthologs_nested$SYMBOL_Nested))
length(unique(hn_gencode_human$SYMBOL_Nested))
#1835 nested genes with ortholog out of 11753 --> 15.61%

#Filter the table to keep the pairs where both the host and the nested genes have an ortholog in mouse
hn_gencode_human_orthologs = hn_gencode_human_orthologs[!is.na(hn_gencode_human_orthologs$Host_mouse_ortholog),]
hn_gencode_human_orthologs = hn_gencode_human_orthologs[!is.na(hn_gencode_human_orthologs$Nested_mouse_ortholog),]
#Number of pair conserved
number = unique(hn_gencode_human_orthologs$pairid)
length(unique(hn_gencode_human$pairid))
#1137 out of 13088 --> 8.6% of the pairs had an ortholog for both the host and the nested gene
write.table(hn_gencode_human_orthologs, "hn_gencode_human_host_and_nestedwithortholog.txt", col.names = T, sep = "\t")


#How many of these pairs of gene are also in a host nested gene organisation in mouse
mouse = merge(hn_gencode_mouse, hn_gencode_human_orthologs[,17:19], by.x = c(12,5), by.y = c(2,3), all.x = T)
mouse_conserved = mouse[!is.na(mouse$pairid.y),]
colnames(mouse_conserved)[17]= "pairid_mouse"
colnames(mouse_conserved)[18]= "pairid_human"
#349 pairs conserved out of 1137
number = unique(mouse_conserved[,17:18])
349/1137 #30.7%

#Find orthologs starting from the MOUSE list to verify and to recover the number of conserved genes
#Use the gene names for mouse as it is what we have

#Merge the tables with the orthologs tables to find which host and nested genes have an ortholog
hn_gencode_mouse_orthologs = merge(hn_gencode_mouse,orthologs[, c(6,4)],by.x = 12,by.y = 1, all.x = T)
colnames(hn_gencode_mouse_orthologs)[18] = "Host_human_ortholog"
hn_gencode_mouse_orthologs = merge(hn_gencode_mouse_orthologs,orthologs[, c(6,4)],by.x = 6,by.y = 1, all.x = T)
colnames(hn_gencode_mouse_orthologs)[19] = "Nested_human_ortholog"

#Calculate the number of host and nested genes in human with at least one ortholog in human
#Host
hn_gencode_mouse_orthologs_host = hn_gencode_mouse_orthologs[!is.na(hn_gencode_mouse_orthologs$Host_human_ortholog),]
length(unique(hn_gencode_mouse_orthologs_host$ENST_Host))
#3653 transcripts host with ortholog
length(unique(hn_gencode_mouse_orthologs_host$SYMBOL_Host))
length(unique(hn_gencode_mouse$SYMBOL_Host))
#3601 host genes with ortholog out of 4801 --> 75%
#Nested
hn_gencode_mouse_orthologs_nested = hn_gencode_mouse_orthologs[!is.na(hn_gencode_mouse_orthologs$Nested_human_ortholog),]
length(unique(hn_gencode_mouse_orthologs_nested$ENST_Nested))
#1483 transcripts nested with ortholog
length(unique(hn_gencode_mouse_orthologs_nested$SYMBOL_Nested))
length(unique(hn_gencode_mouse$SYMBOL_Nested))
#1483 nested genes with ortholog out of 7064 --> 20.99%

#Filter the table to keep the pairs where both the host and the nested genes have an ortholog in human
hn_gencode_mouse_orthologs = hn_gencode_mouse_orthologs[!is.na(hn_gencode_mouse_orthologs$Host_human_ortholog),]
hn_gencode_mouse_orthologs = hn_gencode_mouse_orthologs[!is.na(hn_gencode_mouse_orthologs$Nested_human_ortholog),]

#Number of pair conserved
number = unique(hn_gencode_mouse_orthologs$pairid)
length(unique(hn_gencode_mouse$pairid))
#948 out of 7560 --> 12.5% of the pairs had an ortholog for both the host and the nested gene
write.table(hn_gencode_mouse_orthologs, "hn_gencode_mouse_host_and_nestedwithortholog.txt", col.names = T, sep = "\t")

#How many of these pairs of gene are also in a host nested gene organisation in human
human = merge(hn_gencode_human, hn_gencode_mouse_orthologs[,c(17,18,19)], by.x = c(11,4), by.y = c(2,3), all.x = T)
human_conserved = human[!is.na(human$pairid.y),]
colnames(human_conserved)[18]= "pairid_mouse"
colnames(human_conserved)[17]= "pairid_human"
#349 pairs conserved out of 1357
349/948 #36%

#Verify that the pairs are the same and make a table with the orthologous pairs
all = merge(mouse_conserved, human_conserved, by.x = c(17,18), by.y = c(18,17), all.x = T, all.y = T)
colnames(all)[3] = "SYMBOL_Host_mouse"
colnames(all)[4] = "SYMBOL_Nested_mouse"
colnames(all)[5] = "geneChr_Nested_mouse"
colnames(all)[6] = "geneStart_Nested_mouse"
colnames(all)[7] = "geneEnd_Nested_mouse"
colnames(all)[8] = "ENST_Nested_mouse"
colnames(all)[9] =  "Strand_nested_mouse"
colnames(all)[10] =  "biotype_nested_mouse"
colnames(all)[11] = "geneChr_Host_mouse"
colnames(all)[12] =  "geneStart_Host_mouse"
colnames(all)[13] = "geneEnd_Host_mouse"
colnames(all)[14] = "ENST_Host_mouse"
colnames(all)[15] = "Strand_host_mouse"
colnames(all)[16] = "biotype_host_mouse"
colnames(all)[17] = "Overlap_size_mouse"
colnames(all)[17] = "Overlap_size_mouse"
colnames(all)[18] = "Strandness_mouse"
colnames(all)[19] = "ENST_Host_Human"
colnames(all)[20] = "ENST_Nested_Human"
colnames(all)[21] = "geneChr_Nested_Human"
colnames(all)[22] = "geneStart_Nested_Human"
colnames(all)[23] = "geneEnd_Nested_Human"
colnames(all)[24] = "SYMBOL_Nested_Human"
colnames(all)[25] = "Strand_nested_Human"
colnames(all)[26] = "biotype_nested_Human"
colnames(all)[27] = "geneChr_Host_Human"
colnames(all)[28] = "geneStart_Host_Human"
colnames(all)[29] = "geneEnd_Host_Human"
colnames(all)[30] = "SYMBOL_Host_Human"
colnames(all)[31] =  "Strand_host_Human"
colnames(all)[32] = "biotype_host_Human"
colnames(all)[33] = "Overlap_size_Human"
colnames(all)[34] = "Strandness_human"

#Save the table
write.table(all, file = "Conserved_pairs_mouse_human.txt", col.names = T, sep = "\t")
