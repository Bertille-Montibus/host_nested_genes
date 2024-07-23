#Filtering of the results of the overlap to detect the Host/nested gene pairs#
#Start from the tables generated after overlap#

##################
#### HUMAN #######
##################

#Open the table generated after the overlap
overlaps <- read.delim("V36lift37_overlap.bed", header = F, sep = "\t")

#Add information about the size (geneA and B are names given to the two genes after overlap)
overlaps$Gene_A_size <- overlaps$V3 - overlaps$V2
overlaps$Gene_B_size <- overlaps$V10 - overlaps$V9
overlaps$smallest <- pmin(overlaps$Gene_A_size, overlaps$Gene_B_size)

#Test to see if GeneB is the nested gene using the size
#Is GeneB the smallest ?
Test <- overlaps$Gene_B_size == overlaps$smallest
summary(Test)
# All TRUE so confirm that GeneB is the nested gene

#Add names to the columns
colnames <- c("geneChr_Host","geneStart_Host", "geneEnd_Host", "ENST_Host","SYMBOL_Host","Strand_Host", "biotype_Host",
              "geneChr_Nested","geneStart_Nested", "geneEnd_Nested","SYMBOL_Nested",
              "Overlap_size", "Gene_A_size", "Gene_B_size", "Smallest_size")
colnames(overlaps) <- colnames

#Remove the pairs where the 2 overlapping transcripts are from the same gene
same <- overlaps$SYMBOL_Host == overlaps$SYMBOL_Nested
df <- overlaps[!same,]

#When multiple pairs have the same host and nested gene keep the longest transcript of the host
df <- df[with(df, order(-Gene_A_size)), ]
same <- duplicated(df[,c(5,11)])
df <- df[!same,]

#Remove the genes with same size as they can't be host and nested genes but are probably the same but misannotated
Identical_genes <- df$Gene_A_size == df$Gene_B_size
df <- df[!Identical_genes,]

#Remove the gene with the same start and the same end as they are not proper host and nested genes
same_start <- df$geneStart_Nested == df$geneStart_Host
summary(same_start)
df = df[!same_start,]
same_end <- df$geneEnd_Nested == df$geneEnd_Host
summary(same_end)
df = df[!same_end,]
nestedgene_df = df


#Put back the information for the nested gene to be the same as the information for the host
#Take the file with the transcripts after initial filtering for the biotypes
bed <- read.table('wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed', header = F, sep = "\t")
bed$Gene_size = bed$V3 - bed$V2

#Keep the longest transcript as we did for the host
bed <- bed[with(bed, order(-Gene_size)), ]
same <- duplicated(bed[,5])
bed <- bed[!same,]

#Merge the tables to have complete information for the nested gene
#Merge the table to obtain a complete table as before
hn_table = merge.data.frame(df, bed, by.x = 11, by.y = 5)
hn_table = hn_table[,-c(9:11)]
hn_table = hn_table[,c(13:16,1,17:18,2:9)]

colnames(hn_table)[1] = "geneChr_Nested"
colnames(hn_table)[2] = "geneStart_Nested"
colnames(hn_table)[3] = "geneEnd_Nested"
colnames(hn_table)[4] = "ENST_Nested"
colnames(hn_table)[6] = "Strand_Nested"
colnames(hn_table)[7] = "biotype_Nested"

#Save the initiale table for host and nested genes
write.table(hn_table, file= "Host_nested_genes_GENCODEV36lift37.txt", col.names = T, row.names = F, sep = '\t')


##################
#### MOUSE #######
##################

#Open the table generated after the overlap
overlaps <- read.delim("VM25_overlap.bed", header = F, sep = "\t")

#Add information about the size (geneA and B are names given to the two genes after overlap)
overlaps$Gene_A_size <- overlaps$V3 - overlaps$V2
overlaps$Gene_B_size <- overlaps$V10 - overlaps$V9
overlaps$smallest <- pmin(overlaps$Gene_A_size, overlaps$Gene_B_size)

#Test to see if GeneB is the nested gene using the size
#Is GeneB the smallest ?
Test <- overlaps$Gene_B_size == overlaps$smallest
summary(Test)
# All TRUE so confirm that GeneB is the nested gene

#Add names to the columns
colnames <- c("geneChr_Host","geneStart_Host", "geneEnd_Host", "ENST_Host","SYMBOL_Host","Strand_Host", "biotype_Host",
              "geneChr_Nested","geneStart_Nested", "geneEnd_Nested","SYMBOL_Nested",
              "Overlap_size", "Gene_A_size", "Gene_B_size", "Smallest_size")
colnames(overlaps) <- colnames

#Remove the pairs where the 2 overlapping transcripts are from the same gene
same <- overlaps$SYMBOL_Host == overlaps$SYMBOL_Nested
df <- overlaps[!same,]

#When multiple pairs have the same host and nested gene keep the longest transcript of the host
df <- df[with(df, order(-Gene_A_size)), ]
same <- duplicated(df[,c(5,11)])
df <- df[!same,]

#Remove the genes with same size as they can't be host and nested genes but are probably the same but misannotated
Identical_genes <- df$Gene_A_size == df$Gene_B_size
df <- df[!Identical_genes,]

#Remove the gene with the same start and the same end as they are not proper host and nested genes
same_start <- df$geneStart_Nested == df$geneStart_Host
summary(same_start)
df = df[!same_start,]
same_end <- df$geneEnd_Nested == df$geneEnd_Host
summary(same_end)
df = df[!same_end,]
nestedgene_df = df


#Put back the information for the nested gene to be the same as the information for the host
#Take the file with the transcripts after initial filtering for the biotypes
bed <- read.table('wgEncodeGencodeCompVM25_fileforhost_filtered.bed', header = F, sep = "\t")
bed$Gene_size = bed$V3 - bed$V2

#Keep the longest transcript as we did for the host
bed <- bed[with(bed, order(-Gene_size)), ]
same <- duplicated(bed[,5])
bed <- bed[!same,]

#Merge the tables to have complete information for the nested gene
#Merge the table to obtain a complete table as before
hn_table = merge.data.frame(df, bed, by.x = 11, by.y = 5)
hn_table = hn_table[,-c(9:11)]
hn_table = hn_table[,c(13:16,1,17:18,2:9)]

colnames(hn_table)[1] = "geneChr_Nested"
colnames(hn_table)[2] = "geneStart_Nested"
colnames(hn_table)[3] = "geneEnd_Nested"
colnames(hn_table)[4] = "ENST_Nested"
colnames(hn_table)[6] = "Strand_Nested"
colnames(hn_table)[7] = "biotype_Nested"

#Save the initiale table for host and nested genes
write.table(hn_table, file= "Host_nested_genes_wgEncodeGencodeCompVM25.txt", col.names = T, row.names = F, sep = '\t')
