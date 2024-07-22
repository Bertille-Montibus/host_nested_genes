# Ortholog Analysis # 
# Author: James Cain 
# Paper: Global identification of mammalian host and nested gene pairs reveal tissue-specific transcriptional interplay


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Specific analysis on orthologous pair for heatmap generation # 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

library(dplyr)

# Set new working directory
setwd("~/Local_Documents/ENCODE_file/ENCODE_2/ortholog_info")

# Read files in that you generated from the separate human mouse scripts genocode_human_clean.R and genocode_mouse_clean.R
con_human <- read.csv("human_analysis_conserved_pairids.csv")
con_human_cor <- read.csv("conserved_pair_human_spearman_correlation.csv")
con_mouse <- read.csv("mouse_analysis_conserved_pairids.csv")
con_mouse_cor <- read.csv("conserved_pair_mouse_spearman_correlation.csv")

mouse_ids <- intersect(con_human$pairid_mouse, con_mouse$pairid_mouse)
human_ids <-intersect(con_human$pairid_human, con_mouse$pairid_human)

# make data frame to merge for conserved pair ids 
x <- con_human %>% select(pairid_mouse, conserved_pair_id) %>% left_join(con_mouse)
pair_id_meta <- x %>% dplyr::select(conserved_pair_id, pairid_mouse, pairid_human)

# Merge so can make use pheatmap package  
a <- con_human_cor %>% select(pairid_human=pair_id, cor) %>% left_join(pair_id_meta) %>% 
  rename(cor_human=cor) %>% 
  dplyr::select(conserved_pair_id, cor_human)
b <- con_mouse_cor %>% select(pairid_mouse=pair_id, cor) %>% left_join(pair_id_meta) %>% 
  rename(cor_mouse=cor) %>% 
  dplyr::select(conserved_pair_id, cor_mouse)

# Merge and remove where pair correlation isn't detected in RNA seq data 
cor_human_mouse_conserved <- merge(a,b, by="conserved_pair_id") %>% na.omit()

# pheatmap
library(pheatmap)
rownames(cor_human_mouse_conserved) <- cor_human_mouse_conserved$conserved_pair_id
pheatmap(cor_human_mouse_conserved[,-1], cluster_cols = F)
write.csv(cor_human_mouse_conserved, "cor_human_mouse_conserved.csv")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Specific analysis on non-orthologous org pair for heatmap generation # 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Assessment of correlation differences between non_conserved pairs 
non_con_human_cor <- read.csv("non_conserved_pair_human_spearman_correlation.csv")
non_con_mouse_cor <- read.csv("non_conserved_pair_mouse_spearman_correlation.csv")

# read in tables for non conserved pairs, find pairs ids which AREN'T CONSERVED in organisation 
ortho_pair_org <- read.table("Conserved_pairs_mouse_human.txt")
ortho_pair_any_human <- read.table("hn_gencode_human_host_and_nestedwithortholog.txt")
ortho_pair_any_mouse <- read.table("hn_gencode_mouse_host_and_nestedwithortholog.txt")

# genes which are host nested in human, but NOT in mouse 
non_ortho_org_pair_human <- ortho_pair_any_human %>% filter(!pairid %in% ortho_pair_org$pairid_human)

# genes which are host nested in mouse, but NOT in human 
non_ortho_org_pair_mouse <- ortho_pair_any_mouse %>% filter(!pairid %in% ortho_pair_org$pairid_human)

# Import ID conversion table 
ID_conversion <- read_tsv("ID_conversion_table_hg19")
colnames(ID_conversion)
ID_conversion$hg19.wgEncodeGencodeAttrsV36lift37.transcriptId <-  gsub("_.*", "", ID_conversion$hg19.wgEncodeGencodeAttrsV36lift37.transcriptId)

# Extracting the host gene names 
h <- non_ortho_org_pair_mouse %>% pull(Host_human_ortholog) %>% as.data.frame()
conversion <- merge(
  data.frame(hg19.wgEncodeGencodeAttrsV36lift37.transcriptId=h$.),
  ID_conversion,
  by = "hg19.wgEncodeGencodeAttrsV36lift37.transcriptId") %>% 
  dplyr::select(Host_human_GeneName=hg19.wgEncodeGencodeCompV36lift37.name2, 
                Host_human_ortholog=hg19.wgEncodeGencodeAttrsV36lift37.transcriptId) %>% 
  distinct(Host_human_GeneName,Host_human_ortholog, .keep_all = T)
non_ortho_org_pair_mouse <- left_join(non_ortho_org_pair_mouse, conversion)

# Extracting the nested gene names 
n <- non_ortho_org_pair_mouse %>% pull(Nested_human_ortholog) %>% as.data.frame()
conversion <- merge(
  data.frame(hg19.wgEncodeGencodeAttrsV36lift37.transcriptId=n$.),
  ID_conversion,
  by = "hg19.wgEncodeGencodeAttrsV36lift37.transcriptId") %>% 
  dplyr::select(Nested_human_GeneName=hg19.wgEncodeGencodeCompV36lift37.name2, 
                Nested_human_ortholog=hg19.wgEncodeGencodeAttrsV36lift37.transcriptId) %>% 
  distinct(Nested_human_GeneName,Nested_human_ortholog, .keep_all = T)
non_ortho_org_pair_mouse <- left_join(non_ortho_org_pair_mouse, conversion)

# Filter distinct pairs 
non_ortho_org_pair_mouse <-non_ortho_org_pair_mouse %>% distinct(SYMBOL_Nested,SYMBOL_Host, 
                                                                 Host_human_GeneName, Nested_human_GeneName,
                                                                 .keep_all = T)
non_ortho_org_pair_human <- non_ortho_org_pair_human %>% select(SYMBOL_Host, SYMBOL_Nested, pairid, 
                                                                Host_mouse_ortholog, Nested_mouse_ortholog) %>% 
  distinct(pairid, .keep_all = T)
non_ortho_org_pair_mouse <- non_ortho_org_pair_mouse %>% select(SYMBOL_Host, SYMBOL_Nested, pairid, 
                                                                Host_human_GeneName, Nested_human_GeneName) %>% na.omit() %>% 
  distinct(pairid, .keep_all = T)

# Export tables for use in the expression data analysis present in genocode_human_clean.R and gencode_mouse_clean.R
write.csv(non_ortho_org_pair_human, "conserved_genes_not_host_nested_HUMAN_TO_MOUSE.csv")
write.csv(non_ortho_org_pair_mouse, "conserved_genes_not_host_nested_MOUSE_TO_HUMAN.csv")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### Correlations between non-pairs and pairs, conserved in organisation and non-conserved ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Ensure correct libraries are loaded
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# NON_CONSERVED GENE
h2m <- read.csv("conserved_genes_not_host_nested_HUMAN_TO_MOUSE.csv")
m2h <- read.csv("conserved_genes_not_host_nested_MOUSE_TO_HUMAN.csv")

### ### ### ### ###
  #### 1. hg19 ###
### ### ### ### ###

# Read in data - PAIR IN MOUSE NOT HUMAN
  # NOTE pair id is the respective pair_id of expression
# Files generated within gencode_mouse_clean.R and gencode_human_clean.R
hg19_pseudopair <- read.csv("pseudopair_in_human_not_mouse_cor_hg19.csv")
mm10_pair <- read.csv("pair_in_mouse_not_human_cor_mm10.csv")
p <- rbind(hg19_pseudopair,mm10_pair) %>% ggplot(aes(cor, ortho)) +
  geom_violin()+
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.1) + 
  theme_bw() + 
  scale_y_discrete(labels=c("Matching Pair mm10","Orthologos 'non-pair' hg19")) + 
  ylab("") + 
  xlab("Spearman correlation coefficient")
p


# Export plot comparing pairs and non-pairs in respective species (mouse) in hg19 expression data 
pdf("pseudopairs_in_human.pdf", width = 5, height = 3)
p
dev.off()

# pheatmap 
heatmap <- merge(hg19_pseudopair, mm10_pair, by = "mouse_pair_id")
  # in this instance, x = non_pair and y = real_pair
heatmap <-heatmap %>% dplyr::select(cor.x, cor.y, mouse_pair_id)
rownames(heatmap) <- heatmap$mouse_pair_id
# Remove NAs for heatmaps
heatmap <- heatmap %>% na.omit()
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
pdf("pseudopairs_in_human_heatmap.pdf", width = 2, height = 8)
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, 
         labels_col=c("Orthologos 'non-pair' hg19","Matching Pair mm10"))
dev.off()

# Generation of difference between matched spearman correlation values 
heatmap <- heatmap %>% na.omit()
a <- heatmap %>% mutate(human_minus_mouse= cor.x-cor.y) %>% 
  select(-mouse_pair_id) %>% mutate(comparison = "mouse_not_human")




### ### ### ### ###
#### 2. mm10 ###
### ### ### ### ###

# Read in data - PAIR IN HUMAN NOT MOUSE
# NOTE pair id is the respective pair_id of expression
mm10_pseudopair <- read.csv("pseudopair_in_human_not_mouse_cor_mm10.csv") %>% rename(pair_id =human_pair_id)
hg19_pair <- read.csv("pairs_in_human_not_mouse_cor_hg19.csv")%>% rename(pair_id =mouse_pair_id)
p <- rbind(mm10_pseudopair,hg19_pair) %>% ggplot(aes(cor, ortho)) +
  geom_violin()+
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.1) + 
  theme_bw() + 
  scale_y_discrete(labels=c("Matching Pair hg19","Orthologos 'non-pair' mm10")) + 
  ylab("") + 
  xlab("Spearman correlation coefficient")
p

pdf("pseudopairs_in_mouse.pdf", width = 5, height = 3)
p
dev.off()


# pheatmap
heatmap <- merge(mm10_pseudopair, hg19_pair, by = "pair_id")
# in this instance, x = non_pair and y = real pair
heatmap <-heatmap %>% dplyr::select(cor.x, cor.y, pair_id)
rownames(heatmap) <- heatmap$mouse_pair_id
# Remove NAs for heatmaps
heatmap <- heatmap %>% na.omit()
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)

pdf("pseudopairs_in_mouse_heatmap.pdf", width = 2, height = 8)
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, 
         labels_col=c("Orthologos 'non-pair' mm10","Matching Pair hg19"))
dev.off()

b <- heatmap %>% mutate(human_minus_mouse= cor.x-cor.y) %>% 
  select(-pair_id) %>% mutate(comparison = "human_not_mouse")



### ### ### ### ###
### 3. Read in pairs which are conserved ### 
### ### ### ### ###

# Read in correlation for heatmap 
cor_pairs <- read.csv("cor_human_mouse_conserved.csv")
p <- cor_pairs %>% tidyr::pivot_longer(3:4,names_to = "ortho", values_to = "cor") %>% 
  ggplot(aes(cor, ortho)) +
  geom_violin()+
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.1) + 
  theme_bw() + 
  scale_y_discrete(labels=c("Orthologous Pair hg19","Orthologous Pair mm10")) + 
  ylab("") + 
  xlab("Spearman correlation coefficient")
p

pdf("conserved_pairs_cor_boxplot.pdf", width = 5, height = 3)
p
dev.off()

# pheatmap 
heatmap <-cor_pairs %>% dplyr::select(cor_human, cor_mouse, conserved_pair_id)
rownames(heatmap) <- heatmap$conserved_pair_id
heatmap <- heatmap %>% na.omit()
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)

pdf("conserved_pairs_heatmap.pdf", width = 2, height = 8)
pheatmap(heatmap[,-3], cluster_cols = F, show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, 
         labels_col=c("Orthologous Pair mm10","Orthologous Pair hg19"))
dev.off()
c <- heatmap %>% mutate(human_minus_mouse= cor_human-cor_mouse) %>% 
  select(-conserved_pair_id) %>% mutate(comparison = "conserved")







### ### ### ### ### ### ### ### 
# Correlation difference analysis # 
### ### ### ### ### ### ### ### 

# rbind change daya 
rownames(a) <- NULL
rownames(b) <- NULL
rownames(c) <- NULL
pdf("correlation_differences_conservation.pdf", width = 5, height = 3)
rbind(
  a[,c(3:4)],
  b[,c(3:4)],
  c[,c(3:4)]
) %>% ggplot(aes(human_minus_mouse, comparison)) + 
  geom_violin()+
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.1)
dev.off()

correlation_pair_difference <- rbind(
  a[,c(3:4)],
  b[,c(3:4)],
  c[,c(3:4)]
)

write.csv(correlation_pair_difference, "correlation_pair_differences.csv")


# Levene Test for statistical comparison of distributions # 
library(car)
# Test between conserved AND human_not_mouse 
leveneData <- correlation_pair_difference %>% filter(!comparison == "mouse_not_human")
leveneTest(human_minus_mouse ~ comparison, data = leveneData)
# Levene Test result conserved AND human_not_mouse p = 2.978e-05
# Test between conserved AND mouse_not_human 
leveneData <- correlation_pair_difference %>% filter(!comparison == "human_not_mouse")
leveneTest(human_minus_mouse ~ comparison, data = leveneData)
# Levene Test result conserved AND mouse_not_human p =0.00237
