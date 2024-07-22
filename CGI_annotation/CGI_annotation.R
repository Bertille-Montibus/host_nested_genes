# CGI Annotation # 
# Author: James Cain 
# Paper: Global identification of mammalian host and nested gene pairs reveal tissue-specific transcriptional interplay

### Set working Directory ###

setwd("~/Local_Documents/NestedGene_Project/analysis_July2023/promoter_annotation")

### Load Libraries ### 

library(dplyr)
library(ggplot2)
library(ggsankey)


### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### Hg19 analysis ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### 


# Loading in host/nested gene dataset, filtered removed complex loci and TEC transcripts 
hg19_all <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")

# biotypes to retain for the human genome 
keep <- c("protein_coding","processed_transcript","nonsense_mediated_decay",
          "retained_intron","lncRNA","lincRNA","sense_intronic","snRNA","miRNA","misc_RNA",
          "snoRNA","non_stop_decay", "polymorphic_pseudogene","antisense", "sense_overlapping","rRNA","scRNA","vault_RNA")
load("../total_of_genes /common_biotypes.rda")

# Generation of TSS positions for all transcripts 
sense <- hg19_all %>% filter(V6 == "+") %>% 
  mutate(TSS_end = V2 + 2,
         TSS_start= V2 - 2) %>% 
  select(V1, TSS_start, TSS_end, V5, V7, V4) %>% distinct(V4, .keep_all = T)
antisense <- hg19_all %>% filter(V6 == "-") %>% 
  mutate(TSS_end = V3 + 2,
         TSS_start= V3 - 2) %>% 
  select(V1, TSS_start, TSS_end, V5, V7, V4) %>% distinct(V4, .keep_all = T)

hg19TSS <- rbind(sense,antisense)

# Export bed file
write.table(hg19TSS, "hg19TSS.bed", col.names = F, row.names = F, sep = '\t', quote=FALSE) 

# Within terminal: 
  # $   bedtools intersect -a hg19_cgis.bed -b hg19TSS.bed -F 1.0 -wo > hg19TSS_overlap.bed

# Import TSS overlapping regions with Illingworth defined CGI positions 
hg19_CGIs <- read.table("hg19TSS_overlap.bed")

# Filter to keep specific biotypes 
hg19_CGIs <- hg19_CGIs %>% filter(V9 %in% keep) %>% filter(V9 %in% common_types ) %>% 
  mutate(Promoter = "CGI")
hg19_CGIs_tomerge <-hg19_CGIs %>% select(5,6,7,10,8,9,12)

# Generation of CGI and non-CGI column 
hg19_non_CGIs <- hg19_all %>% filter(!V4 %in% hg19_CGIs$V10) %>% 
  mutate(Promoter = "Non-CGI") %>% select(-6)
hg19_non_CGIs

# Match colnames specifically 
colnames(hg19_non_CGIs) <- colnames(hg19_CGIs_tomerge)

# Bind into master table with CGI annotations
hg19_promoter_annos <- rbind(hg19_CGIs_tomerge, hg19_non_CGIs)

# Distribution of CGIs across all genes # 
library(viridis)
hg19_promoter_annos %>% group_by(V9, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(V9, freq, fill = Promoter)) + 
  geom_col(position = position_stack(reverse = T), color = "white") + 
  coord_flip() + theme_classic() + 
  scale_fill_manual(values = c("forestgreen", "ivory3")) + 
  ylab("Proportion containing promoter type (%)") + xlab("")

# Addition of Host & Nested Gene Information to the table #
hg19_hn <- read.csv("hn_gencode_hg19_IDs.csv")
a <- hg19_promoter_annos %>% filter(V8 %in% hg19_hn$SYMBOL_Host) %>% mutate(hn = "host")
b <- hg19_promoter_annos %>% filter(V8 %in% hg19_hn$SYMBOL_Nested) %>% mutate(hn = "nested")
c <- hg19_promoter_annos %>% mutate(hn = "all_genes")

# Generation of distribution plots per biotype
count_a <- a %>% group_by(hn,V9, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
count_b <-b %>% group_by(hn,V9, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
count_c <-c %>% group_by(hn,V9, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

# Plotting 
plot <- rbind(count_a, count_b, count_c) %>% 
  ggplot(aes(hn, freq, fill = Promoter)) + 
  geom_col(position = position_stack(reverse = T)) + 
  theme_classic() +   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("forestgreen", "ivory3")) + 
  ylab("Proportion containing promoter type (%)") + xlab("")+
  facet_wrap(.~V9) + 
  scale_y_continuous(labels = scales::percent)

# Export data of CGIs per biotype 
write.csv(rbind(count_a, count_b, count_c), "plots/CGIanno_human_data_biotypes.csv")

# Export plot 
pdf("plots/biotypes_CGIanno_human.pdf", width = 9, height = 4)
plot
dev.off()


### Total proportion of Nested or Host genes associated with a CGI ###
count <- c(
  # all genes 
  c %>% filter(Promoter == "CGI") %>% pull(V8) %>% unique() %>% length(),
  c %>% filter(Promoter == "Non-CGI") %>% pull(V8) %>% unique()%>% length(),
  # host genes 
  a %>% filter(Promoter == "CGI") %>% pull(V8) %>% unique()%>% length(),
  a %>% filter(Promoter == "Non-CGI") %>% pull(V8) %>% unique()%>% length(),
  # nested genes
  b %>% filter(Promoter == "CGI") %>% pull(V8) %>% unique()%>% length(),
  b %>% filter(Promoter == "Non-CGI") %>% pull(V8) %>% unique()%>% length()
)

Promoter <- rep(c("CGI", "Non-CGI"),3)

hn <- c("All", "All", "Host", "Host", "Nested", "Nested")

plot <- data.frame(count, Promoter, hn) %>% group_by(hn) %>% 
  mutate(freq = count / sum(count)) %>% 
  ggplot(aes(hn, freq, fill = Promoter)) + 
  geom_bar(stat="identity", position = "stack") + 
  theme_classic() + 
  scale_fill_manual(values = c("forestgreen", "ivory3")) + 
  scale_y_continuous(labels = scales::percent)+ 
  ylab("Proportion containing promoter type") + 
  xlab("")

data <- data.frame(count, Promoter, hn) %>% group_by(hn) %>% 
  mutate(freq = count / sum(count))

data

write.csv(data, "plots/CGIanno_human_data.csv")

pdf("plots/CGIanno_human.pdf", width = 4, height = 3)
plot
dev.off()










### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### mm10 analysis ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### 

### ### ### ### mm10 analysis ### ### ### ###

# Import Data 
mm10_all <- read.table("wgEncodeGencodeCompVM25_fileforhost_filtered.bed")

# Keep 
keep <- c("protein_coding", "lincRNA","sense_overlapping", "antisense","processed_transcript","snRNA",
          "sense_intronic","miRNA","snoRNA","misc_RNA","rRNA","ribozyme", "scaRNA","polymorphic_pseudogene","bidirectional_promoter_lncRNA","macro_lncRNA"
          ,"3prime_overlapping_ncRNA","sRNA", "scRNA","Mt_tRNA","Mt_rRNA")

### ### ### All Genes / Transcripts  ### ### ###

# Generation of TSS Regions 
sense <- mm10_all %>% filter(V3 == "+") %>% filter(V2 != "chrM") %>% 
  mutate(TSS_end = V4 + 2, 
         TSS_start= V4 - 2) %>% 
  select(V2, TSS_start, TSS_end,V3 ,V1, V6,V7) %>% distinct(TSS_start, .keep_all = T) %>% 
  mutate(TSS_start = as.integer(TSS_start)) %>% 
  mutate(TSS_end = as.integer(TSS_end))
antisense <- mm10_all %>% filter(V3 == "-")%>% filter(V2 != "chrM") %>%
  mutate(TSS_end = V5 + 2, 
         TSS_start= V5 - 2) %>% 
  select(V2, TSS_start, TSS_end,V3 ,V1, V6,V7) %>% distinct(TSS_start, .keep_all = T)%>% 
  mutate(TSS_start = as.integer(TSS_start)) %>% 
  mutate(TSS_end = as.integer(TSS_end))
mm10TSS <- rbind(sense, antisense)

# Export bed file
write.table(mm10TSS, "mm10TSS.bed", col.names = F, row.names = F, sep = '\t', quote=FALSE)

# In Terminal: 
  # $ bedtools intersect -a mm10_cgis.bed -b mm10TSS.bed -F 1.0 -wo > mm10TSS_overlap.bed

# See bedtools file for overlap, validated using IGV at HM13
mm10_CGIs <- read.table("mm10TSS_overlap.bed")

# Keep only the select biotypes etc 
mm10_CGIs <- mm10_CGIs %>% filter(V11 %in% keep) %>% filter(V11 %in% common_types ) %>% 
  mutate(Promoter = "CGI")

# Generation of CGI and non-CGI column 
mm10_non_CGIs <- mm10_all %>% filter(!V1 %in% mm10_CGIs$V9) %>%  filter(V7 %in% keep) %>% 
  mutate(Promoter = "Non-CGI") %>% select(-6)

# ID Conversion table 
mm10_IDs <- read.delim("ID_conversion_table_mm10") %>% select(1,5)
colnames(mm10_IDs) <- c("V1", "symbol")
mm10_non_CGIs <- merge(mm10_IDs, mm10_non_CGIs, by = "V1")

# Addition of Gene IDs
a <- mm10_CGIs %>% select(5,6,7,9,10,11,13)
b <- mm10_non_CGIs %>% select(3,5,6,1,2,7,8)
colnames(a) <- colnames(b)
mm10_promoter_annos <- rbind(a, b)


### Addition of Host & Nested Gene Information to the table ### 
mm10_hn <- read.csv("hn_mouse_master.csv")
a <- mm10_promoter_annos %>% filter(symbol %in% mm10_hn$SYMBOL_Host) %>% mutate(hn = "host")
b <- mm10_promoter_annos %>% filter(symbol %in% mm10_hn$SYMBOL_Nested) %>% mutate(hn = "nested")
c <- mm10_promoter_annos %>% mutate(hn = "all_genes")

# Generation of distribution plots per biotype
count_a <- a %>% group_by(hn,V7, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
count_b <-b %>% group_by(hn,V7, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
count_c <-c %>% group_by(hn,V7, Promoter) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

plot <- rbind(count_a, count_b, count_c) %>% 
  ggplot(aes(hn, freq, fill = Promoter)) + 
  geom_col(position = position_stack(reverse = T)) + 
  theme_classic() +   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("forestgreen", "ivory3")) + 
  ylab("Proportion containing promoter type (%)") + xlab("")+
  facet_wrap(.~V7) + 
  scale_y_continuous(labels = scales::percent)
plot
write.csv(rbind(count_a, count_b, count_c), "plots/CGIanno_mouse_data_biotypes.csv")



pdf("plots/biotypes_CGIanno_mouse.pdf", width = 9, height = 4)
plot
dev.off()

### Total proportion of Nested or Host genes associated with a CGI ###
count <- c(
  # all genes 
  c %>% filter(Promoter == "CGI") %>% pull(symbol) %>% unique() %>% length(),
  c %>% filter(Promoter == "Non-CGI") %>% pull(symbol) %>% unique()%>% length(),
  # host genes 
  a %>% filter(Promoter == "CGI") %>% pull(symbol) %>% unique()%>% length(),
  a %>% filter(Promoter == "Non-CGI") %>% pull(symbol) %>% unique()%>% length(),
  # nested genes
  b %>% filter(Promoter == "CGI") %>% pull(symbol) %>% unique()%>% length(),
  b %>% filter(Promoter == "Non-CGI") %>% pull(symbol) %>% unique()%>% length()
)

Promoter <- rep(c("CGI", "Non-CGI"),3)

hn <- c("All", "All", "Host", "Host", "Nested", "Nested")

plot <- data.frame(count, Promoter, hn) %>% group_by(hn) %>% 
  mutate(freq = count / sum(count)) %>% 
  ggplot(aes(hn, freq, fill = Promoter)) + 
  geom_bar(stat="identity", position = "stack") + 
  theme_classic() + 
  scale_fill_manual(values = c("forestgreen", "ivory3")) + 
  scale_y_continuous(labels = scales::percent)+ 
  ylab("Proportion containing promoter type") + 
  xlab("")
plot
data <- data.frame(count, Promoter, hn) %>% group_by(hn) %>% 
  mutate(freq = count / sum(count))

data

write.csv(data, "plots/CGIanno_mouse_data.csv")

pdf("plots/CGIanno_mouse.pdf", width = 4, height = 3)
plot
dev.off()





