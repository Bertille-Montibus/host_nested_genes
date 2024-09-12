# Author: James Cain 
# Paper: Global identification of mammalian host and nested gene pairs reveal tissue-specific transcriptional interplay


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## -- ## -- ## -- ## -- Setting up workspace  -- ## -- ## -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Set your working directory 
setwd("~/Local_Documents/ENCODE_file/ENCODE_2/new_lists/hg19")

# Load required packages
library(tximport)
library(rhdf5)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(EnvStats)

# Define ggplot themes
bartheme <- theme(
  plot.title = element_text(family = "Helvetica", face = "italic", size = 20, hjust = 0.5), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"), 
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 12, face = "bold"), 
  legend.title = element_text(colour = "black", size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, colour = "black"), 
  axis.text.x = element_text(size = 12, colour = "black", angle = 90, vjust = 0.5, hjust = 1)
)

densitytheme <- theme(
  plot.title = element_text(family = "Helvetica", face = "italic", size = 20, hjust = 0.5), 
  panel.grid.major = element_line(colour = "grey"), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"), 
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14, face = "bold"), 
  legend.title = element_text(colour = "black", size = 12, face = "bold"),
  axis.text.y = element_text(face = "bold", size = 14, colour = "black"), 
  axis.text.x = element_text(face = "bold", size = 14, colour = "black")
)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## -- ## -- ## -- ## -- Creating Normalised counts -- ## -- ## -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Importing files that were aligned using kallisto 
sample_id <- dir(file.path("../../kalliso_output_gencode/"))
sample_id
kal_dirs <- file.path("../../kalliso_output_gencode", sample_id)
files <- file.path(kal_dirs, "abundance.h5")
txi_kallisto  <- tximport(files, type = "kallisto", txOut = TRUE)
colnames(txi_kallisto$counts) <- sample_id
counts <- as.data.frame(txi_kallisto$counts)
colnames(counts)
str(counts)

# Convert counts to an integer for use in DESeq2
head(counts)
int <- as.data.frame(sapply(counts,as.integer))
rownames(int) <- rownames(counts)
counts <- int



### ### ### ### ### ### ###
### Metadata Preparation ### 
### ### ### ### ### ### ###

# Generation of coldata for DESeq2  
metadata <- read.csv("../../metadata.csv")
colnames(metadata)
sortedmeta <- metadata %>% filter(Paired.end == "1") %>% select(File.accession, Biosample.term.name)
coldata <- data.frame(row.names=sortedmeta[,1], tissue =as.factor(sortedmeta$Biosample.term.name))
coldata$tissue <- gsub(" ", "_", coldata$tissue)


### ### ### ### ### ### ###
  ### Running DEseq2 ### 
### ### ### ### ### ### ###

# Running DESeq2 
counts <- counts[,rownames(coldata)]
DESeq2Experiment <- as.matrix(counts)

# Generation of DESeq2 object from the count data 
DESeq2Experiment <- DESeqDataSetFromMatrix(countData=DESeq2Experiment, colData=coldata, design=~tissue)
DESeq2Experiment$tissue

# Normalised Count Matrix Generation and export 
DESeq2Experiment <- estimateSizeFactors(DESeq2Experiment)
normalised_counts <- counts(DESeq2Experiment, normalized=TRUE)
write.csv(normalised_counts, "normalisedcounts_gencode.csv")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## -- ## -- ## -- ## -- Loading in Host Nested Gene Data -- ## -- ## -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Loading in the host and nested gene list, protocadherins removed, VDJ, same start and end and same transcript name. 
hn_gencode <- read_delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt")

# Remove version name from transcript ID
hn_gencode$ENST_Host = gsub("\\..*","",hn_gencode$ENST_Host)
hn_gencode$ENST_Nested = gsub("\\..*","",hn_gencode$ENST_Nested)
# Generation of unique ID list for host nested gene pairs 
ids_host <- hn_gencode %>% select(pairid, ENST_Host, Strandness)
ids_nested <- hn_gencode %>% select(pairid, ENST_Nested, Strandness)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## -- ## -- ## -- ## -- Extracting Expression of Host / Nested Genes  -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Removal of version name from transcript ID
rownames(normalised_counts) = gsub("\\..*","",rownames(normalised_counts))

# Generation of host gene expression data 
host_expression <- as.data.frame(normalised_counts)
host_expression$ENST_Host <- rownames(host_expression)
host_expression <- left_join(host_expression, ids_host, by = "ENST_Host")
host_expression <- host_expression %>% filter(ENST_Host %in% ids_host$ENST_Host)

# Generation of nested gene expression data 
nested_expression <- as.data.frame(normalised_counts)
nested_expression$ENST_Nested <- rownames(nested_expression)
nested_expression <- left_join(nested_expression, ids_nested, by = "ENST_Nested")
nested_expression <- nested_expression %>% filter(ENST_Nested %in% ids_nested$ENST_Nested)

# Extraction of expression of shared host and nested gene PAIRS #


# Sorting to match pair ID expression in all data sets
nest <- nested_expression %>% filter(pairid %in% host_expression$pairid ) %>% arrange(pairid) %>% select(-ENST_Nested) %>% as.data.frame()
host <- host_expression %>% filter(pairid %in% nested_expression$pairid ) %>% arrange(pairid) %>% select(-ENST_Host) %>% as.data.frame()

# Transposing the dataframe swapping rows with columns in nest2 and host2 
rownames(nest) <- nest$pairid
nest <- nest[,-29]
nest2 <- data.frame(t(nest[-1]))
colnames(nest2) <- rownames(nest)
rownames(host) <- host$pairid
host <- host[,-29]
host2 <- data.frame(t(host[-1]))
colnames(host2) <- rownames(host)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## -- ## -- ## -- ## -- Analysis of expression data -- ## -- ## -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

## Generation of 'all' genes expression dataset ##

# Import file with all transcripts annotations; 'filtered' file is removal of complex loci and TEC transcripts
gene_annotation <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")
names <- c("chr", "start", "end", "transcript_ID", "symbol", "strand", "biotype")
colnames(gene_annotation) <- names

# Removal of transcript verison ID
gene_annotation$transcript_ID <- gsub("\\..*","",gene_annotation$transcript_ID)

# Annotate Gene IDs with symbols
x <- gene_annotation %>% select(transcript_ID, strand, biotype)
all_expression <- as.data.frame(normalised_counts)
all_expression$transcript_ID <- rownames(all_expression)
all_expression %>% left_join(x) %>% nrow() 
all_expression %>% left_join(x) %>% remove_missing() %>% nrow() 
all_expression <- all_expression %>% left_join(x) %>% remove_missing() %>% pivot_longer(cols = ENCFF576OBS:ENCFF464TEM, names_to = "sample", values_to = "norm_count")

# Merging with tissue-meta, to add metadata information
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
all_expression <- all_expression %>% left_join(tissue_meta)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Large Table of Expression across: All, Host and Nested Genes  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Addition of column identifying group expression data is from 
all_expression2 <-all_expression %>% mutate(hn = "all") %>% select(transcript_ID, biotype, sample, norm_count, tissue, hn)

# Addition of: Host data with biotypes 
to_merge <- gene_annotation %>% select(ENST_Host=transcript_ID, biotype)
host_expression2 <- host_expression %>% pivot_longer(cols = ENCFF576OBS:ENCFF464TEM, names_to = "sample", values_to = "norm_count") %>% mutate(hn = "host") %>% 
  left_join(tissue_meta) %>% left_join(to_merge) %>% 
  select(transcript_ID=ENST_Host, biotype, sample, norm_count, tissue, hn)

# Addition of: Nested data with biotypes
to_merge <- gene_annotation %>% select(ENST_Nested=transcript_ID, biotype)
nested_expression2 <- nested_expression %>% pivot_longer(cols = ENCFF576OBS:ENCFF464TEM, names_to = "sample", values_to = "norm_count")%>% mutate(hn = "nested") %>% 
  left_join(tissue_meta) %>% left_join(to_merge) %>% 
  select(transcript_ID=ENST_Nested, biotype, sample, norm_count, tissue, hn)

# Merging all data and assessing host nested gene expression, only keeping common biotypes 
expression_bigtable <- rbind(all_expression2, host_expression2, nested_expression2)

# Finding total numbers to determine order of bars (total 5859868), proportion of biotypes
expression_bigtable %>% filter(hn == "all")%>% count(biotype) %>% arrange(desc(n)) %>% mutate(proportion=(n/5859868)*100)
expression_bigtable %>% filter(hn == "host") %>% count(biotype) %>% arrange(desc(n)) %>% mutate(proportion=(n/363692)*100)
expression_bigtable %>% filter(hn == "nested") %>% count(biotype) %>% arrange(desc(n)) %>% mutate(proportion=(n/334656)*100)


# Attributing common and unique biotypes across the dataset 
common <- c("protein_coding","lncRNA","retained_intron","processed_transcript","nonsense_mediated_decay","polymorphic_pseudogene", "lincRNA", "antisense")
nested_only <- c("snRNA", "misc_RNA", "miRNA", "snoRNA", "rRNA", "sense_intronic")

# Plotting, ordering bars by the most highly common biotype in ALL genes 
expression_bigtable %>% distinct(transcript_ID, .keep_all = T) %>% count(biotype) %>% arrange(n)
a <- expression_bigtable %>% filter(biotype %in% common) %>% ggplot(aes(factor(biotype, levels =common), log10(norm_count), fill = hn)) + 
  geom_boxplot(width = 0.5) + xlab("biotype") + bartheme + 
  scale_fill_manual(values = c("#fbe4afd6","#1b124870", "#cb6082d6")) + ylab("Normalised Counts (Log10)") + coord_flip()
b <- expression_bigtable %>% filter(biotype %in% nested_only) %>% ggplot(aes(factor(biotype, levels =nested_only), log10(norm_count), fill = hn)) + 
  geom_boxplot(width = 0.5) + xlab("biotype") + bartheme + 
  scale_fill_manual(values = c("#fbe4afd6","#1b124870", "#cb6082d6")) + ylab("Normalised Counts (Log10)")+ coord_flip()

# Exporting plots (CHECK)
pdf("human_plots/hg_expression_bias_common.pdf", height = 5, width = 7)
a
dev.off()
pdf("plots/hg_expression_bias_nestedonly.pdf", height = 3.5, width = 7)
b
dev.off()

# Specifically looking at protein_coding biotype
a <- expression_bigtable %>% filter(biotype %in% common) %>% filter(biotype == "protein_coding") %>% ggplot(aes(tissue, log10(norm_count), fill = hn)) + 
  geom_boxplot(width = 0.5) + xlab("biotype") + bartheme + 
  scale_fill_manual(values = c("#fbe4afd6","#1b124870", "#cb6082d6")) + ylab("Normalised Counts (Log10)") + coord_flip() 

# Printing hg19 correletion plot for host/nested expression 
pdf("plots/hg_expression_bias_common_alltissues_proteincoding.pdf", height = 8, width = 9)
a
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### Export of expression data for RShiny Application  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

tissue_ids <- colnames(nest)
datalist = list()
for (i in tissue_ids) {
  group1 <- nest[[i]]
  group2 <- host[[i]]
  data <- data.frame(nest=group1, host=group2, pair_id = rownames(nest), sample = i)
  datalist[[i]] <- data
}
expression_data <-  do.call(rbind, datalist)
meta <- data.frame(sample = metadata$File.accession, tissue = metadata$Biosample.term.name)
expression_data <- merge(expression_data, meta, by = "sample")
write.csv(expression_data,"expression_data.csv")


### ### ### ### ### ### ### ### ### ### ### ### ### 
  ### Strand Specific Expression Analysis  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### 

# Same orientation, using genes with matched pair ids only. 
same_strand_host <- host_expression %>% filter(Strandness == "same")
same_strand_nested <- nested_expression %>% filter(Strandness == "same") 
same_strand_host <- same_strand_host %>% filter(pairid %in% same_strand_nested$pairid )
same_strand_nested <- same_strand_nested %>% filter(pairid %in% same_strand_host$pairid )
pair_ids <- same_strand_host$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <- i
  group1 <- same_strand_nested %>% 
    filter(Strandness == "same") %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Nested, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group1 <- log(group1+1)
  group2 <- same_strand_host %>% 
    filter(Strandness == "same") %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Host, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(nest2), pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
cor_data_same <-  do.call(rbind, datalist)
correlation_plot_same <- cor_data_same %>% ggplot(aes(cor)) + geom_density(size=1.5)+geom_vline(aes(xintercept=0.7), color="red", linetype="dashed", size=1)+densitytheme+ggtitle("0.7 Red Line")
cor_sort_same <- cor_data_same %>% filter(cor >= 0.7)

# Opposite orientation, using genes with matched pair ids only. 
opposite_strand_host <- host_expression %>% filter(Strandness == "opposite")
opposite_strand_nested <- nested_expression %>% filter(Strandness == "opposite") 
opposite_strand_host <- opposite_strand_host %>% filter(pairid %in% opposite_strand_nested$pairid )
opposite_strand_nested <- opposite_strand_nested %>% filter(pairid %in% opposite_strand_host$pairid )
pair_ids <- opposite_strand_host$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <- i
  group1 <- opposite_strand_nested %>% 
    filter(Strandness == "opposite") %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Nested, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group1 <- log(group1+1)
  group2 <- opposite_strand_host %>% 
    filter(Strandness == "opposite") %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Host, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(nest2), pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
cor_data_opposite <-  do.call(rbind, datalist)
correlation_plot_opposite <- cor_data_opposite %>% ggplot(aes(cor)) + geom_density(size=1.5)+geom_vline(aes(xintercept=0.7), color="red", linetype="dashed", size=1)+densitytheme+ggtitle("0.7 Red Line")
cor_sort_opposite <- cor_data_opposite %>% filter(cor >= 0.7)

# Merging 'same' and 'opposite' strand data into one dataset and performing correlation analysis on them 
a <- cor_data_opposite %>% mutate(Strandness = "opposite")
b <- cor_data_same %>% mutate(Strandness = "same")
cor <- rbind(a,b)

# Split data based on anti, none, cor 
cor_low <- cor %>% filter(cor <= 0)
cor_mid <- cor %>% filter(cor > 0 & cor < 0.5)
cor_high <- cor %>% filter(cor >= 0.5)

# Plotting data 
a <- cor_low %>% 
  ggplot(aes(cor, x = Strandness, fill = Strandness)) + 
  geom_boxplot(size=1) + 
  geom_jitter(alpha = 0.01, width = 0.2, color = "grey30") + 
  scale_color_manual(values = c("firebrick3", "deepskyblue3")) + 
  scale_fill_manual(values = c("firebrick3", "deepskyblue3")) + bartheme+ 
  theme(legend.position = "none")

b <-cor_mid %>% 
  ggplot(aes(cor, x = Strandness, fill = Strandness)) + 
  geom_boxplot(size=1) + 
  geom_jitter(alpha = 0.01, width = 0.2, color = "grey30") + 
  scale_color_manual(values = c("firebrick3", "deepskyblue3")) + 
  scale_fill_manual(values = c("firebrick3", "deepskyblue3")) + bartheme+ 
  theme(legend.position = "none")

c <-cor_high %>% 
  ggplot(aes(cor, x = Strandness, fill = Strandness)) + 
  geom_boxplot(size=1) + 
  geom_jitter(alpha = 0.01, width = 0.2, color = "grey30") + 
  scale_color_manual(values = c("firebrick3", "deepskyblue3")) + 
  scale_fill_manual(values = c("firebrick3", "deepskyblue3")) + bartheme + 
  theme(legend.position = "none")

ggarrange(a,b,c, nrow = 1, ncol = 3)

# Total correlation 
correlation_plot <- cor %>% 
  ggplot(aes(cor, color = Strandness)) + 
  geom_density(size=1.5)+geom_vline(aes(xintercept=0.5), color="grey", linetype="dashed", size=1)+
  geom_density(size=1.5)+geom_vline(aes(xintercept=0), color="grey", linetype="dashed", size=1) + 
  scale_color_manual(values = c("firebrick3", "deepskyblue3")) + theme_classic() + coord_flip()
correlation_plot

cor_sort <- rbind(cor_sort_opposite, cor_sort_same)
cor_sort <- merge(cor_sort, pair_meta, by ="pair_id")

# Exporting the plot 
pdf("plots/hg_correlation_strand_spearman.pdf", height = 3, width = 4)
correlation_plot + coord_flip()
dev.off()

pdf("plots/hg_correlation_strand_spearman_split.pdf", height = 3, width = 4)
ggarrange(a,b,c, nrow = 1, ncol = 3)
dev.off()


# Obtaining stdev for each transcript across all datasets 
a <- expression_bigtable %>% group_by(hn, transcript_ID) %>% summarise(stdev = log10(sd(norm_count))) %>% 
  ggplot(aes(x = hn, y = stdev, fill = hn))+
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  bartheme + 
  scale_fill_manual(values = c("#fbe4afd6","#1b124870", "#cb6082d6")) + ylab("Standard deviation of normalised counts (log10)") + 
  coord_flip()

# Performing separately for match group sizes # 
stdev_table <- expression_bigtable %>% group_by(hn, transcript_ID) %>% summarise(stdev = log10(sd(norm_count))) 
stdev_table %>% filter(hn == "nested") %>% nrow()
stdev_table %>% filter(hn == "host") %>% nrow()

# For host genes 
a <- stdev_table %>% filter(hn == "host") %>% 
  distinct(transcript_ID, .keep_all = T) 
b <- stdev_table %>% filter(hn == "all") %>% 
  filter(!transcript_ID %in% a$transcript_ID) %>% 
  distinct(transcript_ID, .keep_all = T) %>% 
  mutate(hn = "non_host_nested")

# Plotting 
random_stdev_table_host <- rbind(a,b)
p <- random_stdev_table_host %>% ggplot(aes(x = hn,y = stdev, fill = hn)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  coord_flip() + bartheme + xlab("") + theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("#1b124870","grey80")) 

# performing this 1000 times recording the p-value of each KS.test # 
num_iterations <- 1000
ks_test_results <- data.frame(iteration = numeric(num_iterations), 
                              ks_statistic = numeric(num_iterations), 
                              p_value = numeric(num_iterations))
for (i in 1:num_iterations) {
  set.seed(i)
  sample_b <- b %>% sample_n(7920)
  ks_test <- ks.test(a$stdev, sample_b$stdev)
  ks_test_results$iteration[i] <- i
  ks_test_results$ks_statistic[i] <- ks_test$statistic
  ks_test_results$p_value[i] <- ks_test$p.value
}

# Plotting 
ksplot <-ks_test_results %>% arrange(p_value) %>%
  ggplot(aes(x=p_value)) + 
  geom_boxplot() + 
  geom_point(aes(y = 0), position = position_jitter(height = 0.05, width = 0), alpha = 0.5, color = "black") +
  xlim(-0.01,0.1) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = 0.05, linetype = "dashed")+
  xlab("Kolmogorov–Smirnov test p-value")


# Export plots
pdf("plots/human_stdev_host_group.pdf", width =5, height = 3)
p
dev.off()

pdf("plots/human_stdev_host_group_KSTEST.pdf", width =5, height = 1.5)
ksplot
dev.off()

# For nested genes 
a <- stdev_table %>% filter(hn == "nested") %>% 
  distinct(transcript_ID, .keep_all = T) 
b <- stdev_table %>% filter(hn == "all") %>% 
  filter(!transcript_ID %in% a$transcript_ID) %>% 
  distinct(transcript_ID, .keep_all = T) %>% 
  mutate(hn = "non_host_nested")
random_stdev_table_nested <- rbind(a,b)
p <- random_stdev_table_nested %>% ggplot(aes(x = hn,y = stdev, fill = hn)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  coord_flip() + bartheme + xlab("") + theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("#cb6082d6","grey80")) 

# performing this 1000 times recording the p-value of each KS.test # 
num_iterations <- 1000
ks_test_results <- data.frame(iteration = numeric(num_iterations), 
                              ks_statistic = numeric(num_iterations), 
                              p_value = numeric(num_iterations))
for (i in 1:num_iterations) {
  set.seed(i)
  sample_b <- b %>% sample_n(10719)
  ks_test <- ks.test(a$stdev, sample_b$stdev)
  ks_test_results$iteration[i] <- i
  ks_test_results$ks_statistic[i] <- ks_test$statistic
  ks_test_results$p_value[i] <- ks_test$p.value
}

# Plotting 
ksplot <-ks_test_results %>% arrange(p_value) %>%
  ggplot(aes(x=p_value)) + 
  geom_boxplot() + 
  geom_point(aes(y = 0), position = position_jitter(height = 0.05, width = 0), alpha = 0.5, color = "black") +
  xlim(-0.01,0.1) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = 0.05, linetype = "dashed")+
  xlab("Kolmogorov–Smirnov test p-value")
ksplot

# Export plot 
pdf("plots/human_stdev_nested_group.pdf", width =5, height = 3)
p
dev.off()

pdf("plots/human_stdev_nested_group_KSTEST.pdf", width =5, height = 1.5)
ksplot
dev.off()

### ### ### ### ### ### ### ### 
## Tau Analysis ## 
### ### ### ### ### ### ### ### 

# Generation of Tissue specificty Table
# Updating the code for Tau not zero 
x <- as.data.frame(read_tsv("tspex_output/hg19/tau.tsv"))
colnames(x) <- c("ENST", "tau")
x$ENST <- gsub("\\..*","",x$ENST)
y <- hn_gencode %>% filter(pairid %in% rownames(nest)) %>% select(ENST=ENST_Nested, pairid, SYMBOL = SYMBOL_Nested)
z <- hn_gencode %>% filter(pairid %in% rownames(nest)) %>% select(ENST=ENST_Host,pairid, SYMBOL = SYMBOL_Host)
y1 <- x %>% filter(ENST %in% y$ENST)
z2 <- x %>% filter(ENST %in% z$ENST)
test <- left_join(y,y1)
test2 <- left_join(z,z2)
test$hostornested <- "nested"
test2$hostornested <- "host"
tau_table <- rbind(test,test2)
write.csv(tau_table, "hg_tautable.csv")
x$hostornested <- "all_genes"
x$pairid <- "blank"
x$SYMBOL <- "blank"
head(x)
tau_table <- rbind(tau_table,x)
# remove pair ID and unique()
# Checking number of unique colums 
tau_table %>% filter(hostornested == "nested") %>% select(-pairid) %>% unique() %>% nrow()
tau_table %>% filter(hostornested == "host") %>% select(-pairid) %>% unique() %>% nrow()
tau_table %>% filter(hostornested == "all_genes") %>% nrow()
filtered_tau_table <- tau_table %>% select(-pairid) %>% unique()
# Tau == 0 means the gene is expressed in NO TISSUES at all 
# Whereas the stdev calculation removes these valuses
# Attempt plot and remove 0 value 
b <- filtered_tau_table %>% filter(tau > 0) %>% ggplot(aes(x = hostornested,y = tau, fill = hostornested)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  coord_flip() + bartheme + xlab("") + theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("#fbe4afd6","#1b124870", "#cb6082d6")) 


# Comparing tau with random matched sample numbers #

# Matching nested numbers, independent samples met 
a <- tau_table %>% filter(hostornested == "nested") %>% 
  distinct(ENST, .keep_all = T) 
b <- tau_table %>% filter(hostornested == "all_genes") %>% 
  filter(!SYMBOL %in% a$SYMBOL) %>% 
  distinct(ENST, .keep_all = T) %>% 
  mutate(hostornested = "non_host_nested")
random_tau_table_nested <- rbind(a,b)
# Plotting 
plot <- random_tau_table_nested %>% filter(tau > 0) %>% ggplot(aes(x = hostornested,y = tau, fill = hostornested)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  coord_flip() + bartheme + xlab("") + theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("#cb6082d6","grey80")) 
plot
# performing this 1000 times recording the p-value of each KS.test # 
num_iterations <- 1000
ks_test_results <- data.frame(iteration = numeric(num_iterations), 
                              ks_statistic = numeric(num_iterations), 
                              p_value = numeric(num_iterations))
for (i in 1:num_iterations) {
  set.seed(i)
  sample_b <- b %>% sample_n(10679)
  ks_test <- ks.test(a$tau, sample_b$tau)
  ks_test_results$iteration[i] <- i
  ks_test_results$ks_statistic[i] <- ks_test$statistic
  ks_test_results$p_value[i] <- ks_test$p.value
}

# Plotting 
ksplot <-ks_test_results %>% arrange(p_value) %>%
  ggplot(aes(x=p_value)) + 
  geom_boxplot() + 
  geom_point(aes(y = 0), position = position_jitter(height = 0.05, width = 0), alpha = 0.5, color = "black") +
  xlim(-0.01,0.1) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = 0.05, linetype = "dashed")+
  xlab("Kolmogorov–Smirnov test p-value")
ksplot

# EXPORT plot
pdf("plots/human_tau_nested_group.pdf", width =5, height = 3)
plot
dev.off()
pdf("plots/human_tau_nested_groupKSTEST.pdf", width =5, height = 1.5)
ksplot
dev.off()

# Repeat analysis for host genes 
a <- tau_table %>% filter(hostornested == "host") %>% 
  distinct(ENST, .keep_all = T) 
b <- tau_table %>% filter(hostornested == "all_genes") %>% 
  distinct(ENST, .keep_all = T) %>% 
  mutate(hostornested = "non_host_nested")
random_tau_table_host <- rbind(a,b)
# Plotting 
plot <- random_tau_table_host %>% filter(tau > 0) %>% ggplot(aes(x = hostornested,y = tau, fill = hostornested)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  coord_flip() + bartheme + xlab("") + theme(legend.title = element_blank()) + 
  scale_fill_manual(values = c("#1b124870","grey80")) + ggtitle("7397 genes, ks.test = < 2.2e-16")

# performing this 1000 times recording the p-value of each KS.test # 
num_iterations <- 1000
ks_test_results <- data.frame(iteration = numeric(num_iterations), 
                              ks_statistic = numeric(num_iterations), 
                              p_value = numeric(num_iterations))
for (i in 1:num_iterations) {
  set.seed(i)
  sample_b <- b %>% sample_n(7397)
  ks_test <- ks.test(a$tau, sample_b$tau)
  ks_test_results$iteration[i] <- i
  ks_test_results$ks_statistic[i] <- ks_test$statistic
  ks_test_results$p_value[i] <- ks_test$p.value
}

# Plotting 
ksplot <-ks_test_results %>% arrange(p_value) %>%
  ggplot(aes(x=p_value)) + 
  geom_boxplot() + 
  geom_point(aes(y = 0), position = position_jitter(height = 0.05, width = 0), alpha = 0.5, color = "black") +
  xlim(-0.01,0.1) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_vline(xintercept = 0.05, linetype = "dashed")+
  xlab("Kolmogorov–Smirnov test p-value")
ksplot

pdf("plots/human_tau_host_group.pdf", width =5, height = 3)
plot
dev.off()
pdf("plots/human_tau_host_groupKSTEST.pdf", width =5, height = 1.5)
ksplot
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            ## -- Exploring expression data - Orthologs -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Importing Data # 
ortho_pair <- read.table("../../ortholog_info/hn_gencode_human_host_and_nestedwithortholog.txt")
ortho_pair
ortho_pair_org <- read.table("../../ortholog_info/Conserved_pairs_mouse_human.txt")
ortho_pair_org

# find pairs which ARE and ARE NOT conserved
a <- ortho_pair %>% filter(pairid %in% ortho_pair_org$pairid_human) %>% mutate(conserved_pair = "yes")
b <- ortho_pair %>% filter(!pairid %in% ortho_pair_org$pairid_human)%>% mutate(conserved_pair = "no")
ortho_pair_master <- rbind(a,b)


### ### ### ### ### ### ###
# Correlation of value of orthopair #
### ### ### ### ### ### ###

# Modification of correlation code specifically for orthologous pairs
# Pull out expression values for pairs which are conserved 
a <- ortho_pair %>% filter(pairid %in% ortho_pair_org$pairid_human) %>% mutate(conserved_pair = "yes")
h <- host_expression %>% filter(pairid %in% a$pairid)
n <- nested_expression %>% filter(pairid %in% a$pairid)
h <- h %>% filter(pairid %in% n$pairid )
n <- n %>% filter(pairid %in% h$pairid )
pair_ids <- h$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <- i
  group1 <- n %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Nested, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group1 <- log(group1+1)
  group2 <- h %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Host, -pairid,-Strandness)  %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(nest2), pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
cor_data_a <-  do.call(rbind, datalist) %>% mutate(ortho = "conserved_pair")
cor_data_a %>% ggplot(aes(cor)) + geom_density(size=1.5)+densitytheme


### ### ### ### ### ### ### ### ### ###
# Correlation of value pairs which aren't conserved #
### ### ### ### ### ### ### ### ### ###

b <- ortho_pair %>% filter(!pairid %in% ortho_pair_org$pairid_human) %>% mutate(conserved_pair = "no")
h <- host_expression %>% filter(pairid %in% b$pairid)
n <- nested_expression %>% filter(pairid %in% b$pairid)
h <- h %>% filter(pairid %in% n$pairid )
n <- n %>% filter(pairid %in% h$pairid )
pair_ids <- h$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <- i
  group1 <- n %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Nested, -pairid, -Strandness) %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group1 <- log(group1+1)
  group2 <- h %>% 
    filter(pairid == pair_id) %>% 
    select(-ENST_Host, -pairid,-Strandness)  %>% 
    t() %>% as.data.frame() %>% pull(V1)
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(nest2), pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
cor_data_b <-  do.call(rbind, datalist) %>% mutate(ortho = "non_conserved_pair")
cor_data_b %>% ggplot(aes(cor)) + geom_density(size=1.5)+densitytheme

# Plotting 
plot <- rbind(cor_data_a%>% drop_na(cor),
              cor_data_b%>% drop_na(cor)) %>% 
  ggplot(aes(ortho,cor)) + 
  geom_violin(width=.5)+
  geom_boxplot(size=1, alpha=0.1, width = 0.2) +
  theme_bw() + coord_flip()

pdf("plots/human_conserved_correlation_ALLPAIRS.pdf", width = 6, height = 3)
plot
dev.off()

# Exporting conserved pair id 
conserved_pair_id <- ortho_pair_org %>% filter(pairid_human %in% cor_data_a$pair_id) %>% 
  mutate(conserved_pair_id = paste0("conserved_pairid_",row_number()))

# Exporting tables 
write.csv(conserved_pair_id,"../../ortholog_info/human_analysis_conserved_pairids.csv")
write.csv(cor_data_a,"../../ortholog_info/conserved_pair_human_spearman_correlation.csv")
write.csv(cor_data_b,"../../ortholog_info/non_conserved_pair_human_spearman_correlation.csv")



### ### ### ### ### ### ### ### ### ###
# Orthologous Expression Analysis  #
### ### ### ### ### ### ### ### ### ###

# Use of files generated in ortholog_convervation_analysis.R script also present in this folder 
# Extraction of expression values to conserved pairs which are either: host/nested or not, between species 

## USAGE OF Untitled.R location in ortholog_info folder to generate files with expression data ortho and non-orthos ## 
non_ortho_org_pair_human <- read.csv("../../ortholog_info/conserved_genes_not_host_nested_HUMAN_TO_MOUSE.csv")
non_ortho_org_pair_mouse <- read.csv("../../ortholog_info/conserved_genes_not_host_nested_MOUSE_TO_HUMAN.csv")

### ### ### ### ### ### ### ### ### ### ### ###
# Set 1. Genes which are host nested in human but not in mouse 
### ### ### ### ### ### ### ### ### ### ### ###

## Extraction of expression of genes ##

# Import conversion table # 
hg19_ids <- read_tsv("../../ortholog_info/ID_conversion_table_hg19") %>% select(1,5)
colnames(hg19_ids) <- c("ENST", "GeneName")
hg19_ids$ENST <- gsub("\\..*", "", hg19_ids$ENST)

# Import longest transcript data - filter gene names exclusive to these 
longest_transcipt <- read.table("../../ortholog_info/GENCODEV36lift37_longest_transcript.txt")
longest_transcipt$V4 <- gsub("\\..*", "", longest_transcipt$V4 )
hg19_ids <- hg19_ids %>% filter(ENST %in% longest_transcipt$V4 )


# Generation of psuedo-host gene expression data 
# a_expression = pseudo_host 
a_expression <- as.data.frame(normalised_counts)
a_expression$ENST <- rownames(a_expression)
a_expression <- left_join(a_expression, hg19_ids, by = "ENST")
a_expression <- a_expression %>% filter(GeneName %in% non_ortho_org_pair_mouse$Host_human_GeneName)
a_expression <- merge(a_expression,
                      non_ortho_org_pair_mouse %>% select(GeneName=Host_human_GeneName,pairid),
                      by = "GeneName")


# Generation of b gene expression data 
# b_expression = pseudo_nested
b_expression <- as.data.frame(normalised_counts)
b_expression$ENST <- rownames(b_expression)
b_expression <- left_join(b_expression, hg19_ids, by = "ENST")
b_expression <- b_expression %>% filter(GeneName %in% non_ortho_org_pair_mouse$Nested_human_GeneName)
b_expression <- merge(b_expression,
                      non_ortho_org_pair_mouse %>% select(GeneName=Nested_human_GeneName,pairid),
                      by = "GeneName")

# Matched pair ID expression in all data sets is 2341, sort data 
# Matched orthologous pairs is 
a_expression <- a_expression %>% unique()%>% filter(pairid %in% b_expression$pairid ) %>% arrange(pairid) %>% select(-ENST) %>% as.data.frame()
b_expression <- b_expression %>% unique()%>% filter(pairid %in% a_expression$pairid ) %>% arrange(pairid) %>% select(-ENST) %>% as.data.frame() 

# Transposing the dataframe swapping rows with columns in nest2 and host2 
rownames(a_expression) <- a_expression$pairid
a_expression <- a_expression[,-29]
a_expression2 <- data.frame(t(a_expression[-1]))
colnames(a_expression2) <- rownames(a_expression)

rownames(b_expression) <- b_expression$pairid
b_expression <- b_expression[,-29]
b_expression2 <- data.frame(t(b_expression[-1]))
colnames(b_expression2) <- rownames(b_expression)

## Extraction of expression of genes ##
# correlation between 'pseudo-pairs' in the mouse genome
pair_ids <- b_expression$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <-i
  group1 <- b_expression %>% 
    filter(pairid == pair_id) %>% 
    select( -pairid) %>% 
    t() %>% as.data.frame() %>% pull(pair_id)
  group1 <- as.numeric(group1[-1])
  group1 <- log(group1+1)
  group2 <- a_expression %>% 
    filter(pairid == pair_id) %>% 
    select( -pairid) %>% 
    t() %>% as.data.frame() %>% pull(pair_id)
  group2 <- as.numeric(group2[-1])
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(b_expression2)[1:27], pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
pseudopairs_in_mouse_not_human_cor <-  do.call(rbind, datalist) %>% mutate(ortho = "pseudopairs_in_mouse_not_human")
colnames(pseudopairs_in_mouse_not_human_cor)[1] <- "mouse_pair_id"

# Quick histogram shows relative positive correlation
hist(pseudopairs_in_mouse_not_human_cor$cor)

# Export table # 
write.csv(pseudopairs_in_mouse_not_human_cor, "../../ortholog_info/pseudopair_in_human_not_mouse_cor_hg19.csv")



### ### ### ### ### ### ### ### ### ### ### ###
# Set 2. Genes which are host nested in mouse but not in human 
### ### ### ### ### ### ### ### ### ### ### ###

## Extraction of expression of genes ##

# Import conversion table # 
hg19_ids <- read_tsv("../../ortholog_info/ID_conversion_table_hg19") %>% select(1,5)
colnames(hg19_ids) <- c("ENST", "GeneName")
hg19_ids$ENST <- gsub("\\..*", "", hg19_ids$ENST)

# Import longest transcript data - filter gene names exclusive to these 
longest_transcipt <- read.table("../../ortholog_info/GENCODEV36lift37_longest_transcript.txt")
longest_transcipt$V4 <- gsub("\\..*", "", longest_transcipt$V4 )
hg19_ids <- hg19_ids %>% filter(ENST %in% longest_transcipt$V4 )

# Generation of psuedo-host gene expression data 
# a_expression = pseudo_host 
a_expression <- as.data.frame(normalised_counts)
a_expression$ENST <- rownames(a_expression)
a_expression <- left_join(a_expression, hg19_ids, by = "ENST")
a_expression <- a_expression %>% filter(GeneName %in% non_ortho_org_pair_human$SYMBOL_Host)
a_expression <- merge(a_expression,
                      non_ortho_org_pair_human %>% select(GeneName=SYMBOL_Host,pairid),
                      by = "GeneName")

# Generation of b gene expression data 
# b_expression = pseudo_nested
b_expression <- as.data.frame(normalised_counts)
b_expression$ENST <- rownames(b_expression)
b_expression <- left_join(b_expression, hg19_ids, by = "ENST")
b_expression <- b_expression %>% filter(GeneName %in% non_ortho_org_pair_human$SYMBOL_Nested)
b_expression <- merge(b_expression,
                      non_ortho_org_pair_human %>% select(GeneName=SYMBOL_Nested,pairid),
                      by = "GeneName")

# Matched pair ID expression in all data sets is 2341, sort data 
# Matched orthologous pairs is 
a_expression <- a_expression %>% unique()%>% filter(pairid %in% b_expression$pairid ) %>% arrange(pairid) %>% select(-ENST) %>% as.data.frame()
b_expression <- b_expression %>% unique()%>% filter(pairid %in% a_expression$pairid ) %>% arrange(pairid) %>% select(-ENST) %>% as.data.frame() 

# Transposing the dataframe swapping rows with columns in nest2 and host2 
rownames(a_expression) <- a_expression$pairid
a_expression <- a_expression[,-29]
a_expression2 <- data.frame(t(a_expression[-1]))
colnames(a_expression2) <- rownames(a_expression)

rownames(b_expression) <- b_expression$pairid
b_expression <- b_expression[,-29]
b_expression2 <- data.frame(t(b_expression[-1]))
colnames(b_expression2) <- rownames(b_expression)

## Extraction of expression of genes ##
# correlation between 'pseudo-pairs' in the mouse genome
pair_ids <- b_expression$pairid
datalist = list()
plot_list = list()
pair_meta <- hn_gencode %>% select(pairid, SYMBOL_Nested, SYMBOL_Host)
tissue_meta <- metadata %>% select(sample =File.accession, tissue = Biosample.term.name)
colnames(pair_meta) <- c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")
for (i in pair_ids) {
  pair_id <-i
  group1 <- b_expression %>% 
    filter(pairid == pair_id) %>% 
    select( -pairid) %>% 
    t() %>% as.data.frame() %>% pull(pair_id)
  group1 <- as.numeric(group1[-1])
  group1 <- log(group1+1)
  group2 <- a_expression %>% 
    filter(pairid == pair_id) %>% 
    select( -pairid) %>% 
    t() %>% as.data.frame() %>% pull(pair_id)
  group2 <- as.numeric(group2[-1])
  group2 <- log(group2+1)
  data <- data.frame(nest=group1, host=group2, sample = rownames(b_expression2)[1:27], pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
    stat_cor(method = "pearson", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
}
pairs_in_human_not_mouse <-  do.call(rbind, datalist) %>% mutate(ortho = "pairs_in_human_not_mouse")
colnames(pairs_in_human_not_mouse)[1] <- "mouse_pair_id"

# Quick histogram shows relative positive correlation
hist(pairs_in_human_not_mouse$cor)

# Export table # 
write.csv(pairs_in_human_not_mouse, "../../ortholog_info/pairs_in_human_not_mouse_cor_hg19.csv")




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ## -- ## -- Exploring expression data - Heatmaps -- ## -- ## 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Generation of data for use in pheatmap package 
all_genes_expression <- as.data.frame(normalised_counts)
all_genes_expression$ENST <- rownames(all_genes_expression)
all_genes_expression <- all_genes_expression %>% filter(!ENST %in% ids_host$ENST_Host ) %>% 
  filter(!ENST %in% ids_nested$ENST_Nested ) 
all_genes_expression <- all_genes_expression[,-29]

# Generation of pheatmap metadata 
meta
meta_heatmap <- as.data.frame(colnames(all_genes_expression))
colnames(meta_heatmap) <- "sample"
meta_heatmap <- merge(meta_heatmap, meta, by = "sample")
rownames(meta_heatmap) <- meta_heatmap$sample
meta_heatmap <- data.frame(row.names =  meta_heatmap$sample, 
                           tissue = meta_heatmap$tissue)
# Load in the library 
library(pheatmap)

# taking expression from stat.cor pair ids (correlation >0.7)
cor_nested_exp <- nested_expression %>% filter(pairid %in% cor_sort$pair_id) %>% arrange(pairid)
cor_nested_exp_heatmap <- cor_nested_exp[,-c(29,30)]
rownames(cor_nested_exp_heatmap) <- cor_nested_exp$pairid
list_of_genes <- cor_nested_exp$pairid
cor_host_exp <- host_expression %>% filter(pairid %in% list_of_genes) %>% arrange(pairid)
cor_host_exp_heatmap <- cor_host_exp[,-c(29,30)]
rownames(cor_host_exp_heatmap) <- cor_host_exp$pairid

# Customising pheatmap
  # Arrange based on tissue
# Generation of pheatmap with the same order of rows (genes)
meta_heatmap <- meta_heatmap %>% arrange(tissue)
newCols <- colorRampPalette(grDevices::rainbow(length(unique(meta_heatmap$tissue))))
mycolors <- newCols(length(unique(meta_heatmap$tissue)))
names(mycolors) <- unique(meta_heatmap$tissue)
mycolors <- list(tissue = mycolors)
cor_host_exp_heatmap <- cor_host_exp_heatmap[rownames(meta_heatmap)]
brks <- seq(0,14,length.out=100)  
one <- pheatmap(log(cor_host_exp_heatmap+1),cluster_cols = F,annotation_col =meta_heatmap, clustering_method = "average",show_rownames=F,border_color = FALSE, annotation_colors = mycolors, breaks=brks)
ph <- pheatmap(log(cor_host_exp_heatmap+1),cluster_cols = F,annotation_col =meta_heatmap, clustering_method = "average")
order <- rownames(cor_host_exp_heatmap[ph$tree_row[["order"]],])
cor_nested_exp_heatmap$pairid <- rownames(cor_nested_exp_heatmap)
cor_nested_exp_heatmap <- cor_nested_exp_heatmap[match(order, cor_nested_exp_heatmap$pairid),]
cor_nested_exp_heatmap <- cor_nested_exp_heatmap[,-29]
cor_nested_exp_heatmap <- cor_nested_exp_heatmap[rownames(meta_heatmap)]
two <- pheatmap(log(cor_nested_exp_heatmap+1),cluster_cols = F,cluster_rows = F,annotation_col =meta_heatmap, clustering_method = "average",show_rownames=F,border_color = FALSE, annotation_colors = mycolors, breaks=brks)
cor_nested_exp_heatmap
cor_host_exp_heatmap

# Exporting Heatmaps 
pdf("plots/hg_host_heatmap.pdf", height = 10, width = 10)
one
dev.off()

pdf("plots/hg_nested_heatmap.pdf", height = 10, width = 10)
two
dev.off()


# Idenficiation of clusters containing testis specific genes 
out <- pheatmap(log(cor_host_exp_heatmap+1), cluster_cols = F,annotation_col =meta_heatmap, clustering_method = "average", cutree_rows = 100)
out
clusters <- as.data.frame(sort(cutree(out$tree_row, k=100))) 
colnames(clusters) <- "cluster"
clusters$id <-rownames(clusters)
pdf("plots/cluster_identification.pdf", height = 100, width = 10)
out
dev.off()

# Extraction of testis specific genes 
testis_pairid <- clusters %>% filter(cluster == 1)
testis_pairid <- rownames(testis_pairid)
testis_hn_genes <- hn_gencode %>% filter(pairid %in% testis_pairid )
write.csv(testis_hn_genes, file = "testis_hn_genes_logged.csv")

# Correlation with the tau 
tau_testis <- tau_table %>% filter(pairid %in% testis_hn_genes$pairid) %>% select(pair_id=pairid, tau, hostornested)
cor_value <-cor_sort %>% filter(pair_id %in% testis_hn_genes$pairid)
a <- merge(tau_testis, cor_value, by = "pair_id") %>% ggplot(aes(tau,cor, color = hostornested)) + 
  geom_point(size =2,alpha = 0.5) +
  geom_smooth(size = 1,method=lm, aes(fill=hostornested)) + 
  bartheme + 
  scale_color_manual(values = c("#1b124870", "#cb6082d6")) + 
  scale_fill_manual(values = c("#1b124870", "#cb6082d6")) + 
  theme(legend.position="none") + 
  scale_y_continuous(limits = c(0.65,1)) +
  scale_x_continuous(limits = c(0.4,1))

# For all host nested genes that are highly correlated 
tau_hn <- tau_table %>% select(pair_id=pairid, tau, hostornested)
cor_value <-cor_sort 
b_data <- left_join(tau_hn, cor_value) %>% distinct(pair_id, hostornested, tau, cor,.keep_all = T)
b <- b_data %>% filter(hostornested != "all_genes")%>% ggplot(aes(tau,cor, color = hostornested)) + 
  geom_point(size =2, alpha = 0.5) +
  geom_smooth(size = 1,method=lm, aes(fill=hostornested), se= F) + 
  bartheme + 
  scale_color_manual(values = c("#1b124870", "#cb6082d6")) + 
  scale_fill_manual(values = c("#1b124870", "#cb6082d6")) + 
  theme(legend.position="none")  + 
  scale_y_continuous(limits = c(0.65,1)) +
  scale_x_continuous(limits = c(0.4,1)) 
b

# For all host nested genes that are not highly correlated
tau_hn <- tau_table %>% select(pair_id=pairid, tau, hostornested)
cor_value <- merge(cor, pair_meta, by ="pair_id")
c_data <- left_join(tau_hn, cor_value) %>% distinct(pair_id, hostornested, tau, cor,.keep_all = T)
c <- c_data %>% filter(hostornested != "all_genes")%>%ggplot(aes(tau,cor, color = hostornested)) + 
  geom_point(size =2, alpha = 0.1) +
  geom_smooth(size = 1,method=lm, aes(fill=hostornested),se= F)+ 
  bartheme + 
  scale_color_manual(values = c("#1b124870", "#cb6082d6")) + 
  scale_fill_manual(values = c("#1b124870", "#cb6082d6")) + 
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0.4,1)) 
c

##

pdf(file = "plots/hg_cor_tau_relationship_2_highcor.pdf", width = 4, height = 4)
b
dev.off()


pdf(file = "plots/hg_cor_tau_relationship_3_allgenes.pdf", width = 4, height = 4)
c
dev.off()


# SAVING TABLE FOR THE MASTER TABLE COLLATION 
write.csv(cor_value, file = "../../../../NestedGene_Project/matser_table_generation/cor_value_human.csv")
write.csv(tau_table, file = "../../../../NestedGene_Project/matser_table_generation/tau_value_human.csv")
write.csv(variance, file = "../../../../NestedGene_Project/matser_table_generation/stdev_value_human.csv")

