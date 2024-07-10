# Analysis of the biotypes of the host and nested genes #


library("ggplot2")
library("dplyr")
library(RColorBrewer)
library(pals)
library("ComplexHeatmap")
library("viridis")

##################
#### HUMAN #######
##################

#Open the table generated after filtering of the pairs and orienatation was added
hn_gencode <- read.table("Host_nested_genes_GENCODEV36lift37_strand.txt",sep = '\t', header = T)

#Open the reference table with all the transcripts we used at the beginning
all <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")


#Add the strand information to the biotypes
#Nested genes
nested_strand = hn_gencode[,c(5,16)]
#Remove duplications when the nested gene is always in the same orientation when nested in multiple hosts
nested_strand = nested_strand[!duplicated(nested_strand),]
#Find the nested genes which can be in both orientation when nested in multiple hosts
both  = nested_strand[duplicated(nested_strand[,1]),]
#Remove them from the main table to add them back with the right annotation
nested_strand = nested_strand[!(nested_strand$SYMBOL_Nested %in% both$SYMBOL_Nested),]
both$Strandness = "Both"
nested_strand = rbind(nested_strand, both)
#Add the orientation information to the table with the biotypes of all selected transcripts for the analysis
#Recover the biotypes of the transcripts
biotype_nested_count = merge.data.frame(nested_strand, all[,c(5,7)], by.x = 1, by.y = 1)
biotype_nested_count_2 = biotype_nested_count %>% dplyr::count(V7, Strandness) %>% as.data.frame()

#Host genes
host_strand = hn_gencode[,c(12,16)]
#Remove duplications when the host gene is always in the same orientation when it contains multiple nested genes
host_strand = host_strand[!duplicated(host_strand),]
#Find the host genes which can be in both orientation when it contains multiple nested genes
both  = host_strand[duplicated(host_strand[,1]),]
#Remove them from the main table to add them back with the right annotation
host_strand = host_strand[!(host_strand$SYMBOL_Host %in% both$SYMBOL_Host),]
both$Strandness = "Both"
host_strand = rbind(host_strand, both)
#Add the orientation information to the table with the biotypes of all selected transcripts for the analysis
#Recover the biotypes of the transcripts
biotype_host_count = merge.data.frame(host_strand, all[,c(5,7)], by.x = 1, by.y = 1)
biotype_host_count_2 = biotype_host_count %>% dplyr::count(V7, Strandness) %>% as.data.frame()


#Make the bar plot from the data with the counts
ggplot(biotype_nested_count_2, aes(fill=Strandness, y=n, x= V7)) +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))


ggplot(biotype_host_count_2, aes(fill=Strandness, y=n, x= V7)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

#Comparison of the different distributions of biotypes (host, nested and all transcripts) to see if there is a difference
#Start with the biotypes found both in mouse and human

#List of the common biotypes between mouse and human
common_types = c("antisense","lincRNA","miRNA","misc_RNA","polymorphic_pseudogene","processed_transcript","protein_coding","rRNA",
                  "sense_intronic","sense_overlapping","snRNA","snoRNA","lncRNA")

#Calculate the percentage for all transcripts for these biotypes
a <- all %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "all") %>%
  dplyr::rename(biotype=V7)
sum(a$n)

#Calculate the percentage for host transcripts for these biotypes
b <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Host) %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "host") %>%
  dplyr::rename(biotype=V7)
sum(b$n)

#Calculate the percentage for nested transcripts for these biotypes
c <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Nested) %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "nested") %>%
  dplyr::rename(biotype=V7)
sum(c$n)

#Pool the datasets
data <- rbind(a,b,c)
remove <- data %>% filter(category == "all" & n < 5) %>% pull(biotype)
# Remove biotypes which are present in less than five transcripts
data <-data %>% filter(!biotype %in%  remove)
data$category <- factor(data$category, levels=c("nested", "host", "all"))

#Plot the data
data %>%  filter(biotype %in% common_types ) %>%
  ggplot(aes(category, prop, fill = biotype)) +
  geom_col(position = position_stack(reverse = T), color = "white") +
  coord_flip() + theme_classic() +
  scale_fill_viridis(option="magma", discrete = T) +
  ylab("Proportion total transcripts (%)") + xlab("")

#Biotypes specific to human
specific_types = setdiff(unique(all$V7), common_types)

#Calculate the percentage for all transcripts for these biotypes
a <- all %>% dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "all") %>%
    dplyr::rename(biotype=V7)
sum(a$n)

#Calculate the percentage for host transcripts for these biotypes
b <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Host) %>%dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "host") %>%
    dplyr::rename(biotype=V7)
sum(b$n)

#Calculate the percentage for nested transcripts for these biotypes
c <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Nested) %>% dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "nested") %>%
    dplyr::rename(biotype=V7)
sum(c$n)

#Pool the datasets
data <- rbind(a,b,c)
data$category <- factor(data$category, levels=c("nested", "host", "all"))

#Plot the data
data %>%  ggplot(aes(category, prop, fill = biotype)) +
    geom_col(position = position_stack(reverse = T), color = "white") +
    coord_flip() + theme_classic() +
    scale_fill_viridis(option="cividis", discrete = T) +
    ylab("Proportion total transcripts (%)") + xlab("")

#Heatmap to see the biotype correspondence between Host/nested genes

#Add an ID to the host/nested gene pairs
hn_gencode$pairid <- paste(seq.int(nrow(hn_gencode)),"_ID",sep = "")

#Save the table
write.table(hn_gencode, file = "Host_nested_genes_GENCODEV36lift37_strand_ID.txt", sep = '\t')

#Use the orientation already defined as an annotation
hn_gencode$Strandness = as.factor(hn_gencode$Strandness)
Anno = hn_gencode[,c(17,16)]
row.names(Anno) = Anno$pairid

#Make a table with the biotypes for the heatmap
Gencode_biotype = hn_gencode[,c(14,7,17)]
row.names(Gencode_biotype)=Gencode_biotype$pairid
Gencode_biotype = Gencode_biotype[,-3]

#Make the biotypes as factor to have one color per biotype
Gencode_biotype$biotype_Host = as.factor(Gencode_biotype$biotype_Host)
Gencode_biotype$biotype_Nested = as.factor(Gencode_biotype$biotype_Nested)

#Transform the table as a matrix for the heatmap function
Gencode_biotype_h = as.matrix(Gencode_biotype)

#Create the objects needed for the annotation of the heatmap
row_ha = rowAnnotation(orientation = Anno$Strandness, col = list(orientation= c("same"="grey50", "opposite"="coral1")))
col =polychrome(n = nlevels(Gencode_biotype$biotype_Host))
names(col) = levels(Gencode_biotype$biotype_Host)
biotype_anno = rowAnnotation(biotype = Gencode_biotype$biotype_Host, col = list(biotype = col))

#Plot the heatmap
Heatmap(Gencode_biotype_h, name = "Biotype", col = inferno(length(test)),
        column_title = "Heatmap_Biotype",
        row_split = Gencode_biotype_h[,1:2],
        row_title = NULL,
        show_row_names = FALSE,
        row_gap = unit(0, "mm"),
        right_annotation = row_ha,
        left_annotation = biotype_anno)

##################
#### MOUSE #######
##################

#Open the table generated after filtering of the pairs and orienatation was added
hn_gencode <- read.table("Host_nested_genes_wgEncodeGencodeCompVM25_strand.txt",sep = '\t', header = T)

#Open the reference table with all the transcripts we used at the beginning
all <- read.table("wgEncodeGencodeCompVM25_fileforhost_filtered.bed")


#Add the strand information to the biotypes
#Nested genes
nested_strand = hn_gencode[,c(5,16)]
#Remove duplications when the nested gene is always in the same orientation when nested in multiple hosts
nested_strand = nested_strand[!duplicated(nested_strand),]
#Find the nested genes which can be in both orientation when nested in multiple hosts
both  = nested_strand[duplicated(nested_strand[,1]),]
#Remove them from the main table to add them back with the right annotation
nested_strand = nested_strand[!(nested_strand$SYMBOL_Nested %in% both$SYMBOL_Nested),]
both$Strandness = "Both"
nested_strand = rbind(nested_strand, both)
#Add the orientation information to the table with the biotypes of all selected transcripts for the analysis
#Recover the biotypes of the transcripts
biotype_nested_count = merge.data.frame(nested_strand, all[,c(5,7)], by.x = 1, by.y = 1)
biotype_nested_count_2 = biotype_nested_count %>% dplyr::count(V7, Strandness) %>% as.data.frame()

#Host genes
host_strand = hn_gencode[,c(12,16)]
#Remove duplications when the host gene is always in the same orientation when it contains multiple nested genes
host_strand = host_strand[!duplicated(host_strand),]
#Find the host genes which can be in both orientation when it contains multiple nested genes
both  = host_strand[duplicated(host_strand[,1]),]
#Remove them from the main table to add them back with the right annotation
host_strand = host_strand[!(host_strand$SYMBOL_Host %in% both$SYMBOL_Host),]
both$Strandness = "Both"
host_strand = rbind(host_strand, both)
#Add the orientation information to the table with the biotypes of all selected transcripts for the analysis
#Recover the biotypes of the transcripts
biotype_host_count = merge.data.frame(host_strand, all[,c(5,7)], by.x = 1, by.y = 1)
biotype_host_count_2 = biotype_host_count %>% dplyr::count(V7, Strandness) %>% as.data.frame()

#Make the bar plot from the data with the counts
ggplot(biotype_nested_count_2, aes(fill=Strandness, y=n, x= V7)) +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

ggplot(biotype_host_count_2, aes(fill=Strandness, y=n, x= V7)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

#Comparison of the different distributions of biotypes (host, nested and all transcripts) to see if there is a difference
#Start with the biotypes found both in mouse and human

#List of the common biotypes between mouse and human
common_types = c("antisense","lincRNA","miRNA","misc_RNA","polymorphic_pseudogene","processed_transcript","protein_coding","rRNA",
                    "sense_intronic","sense_overlapping","snRNA","snoRNA","lncRNA")

#Calculate the percentage for all transcripts for these biotypes
a <- all %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "all") %>%
  dplyr::rename(biotype=V7)
sum(a$n)

#Calculate the percentage for host transcripts for these biotypes
b <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Host) %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "host") %>%
  dplyr::rename(biotype=V7)
sum(b$n)

#Calculate the percentage for nested transcripts for these biotypes
c <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Nested) %>% dplyr::count(V7) %>% filter(V7 %in% common_types ) %>%
  mutate(prop=n/sum(n)*100) %>%
  mutate(category = "nested") %>%
  dplyr::rename(biotype=V7)
sum(c$n)

#Pool the datasets
data <- rbind(a,b,c)
remove <- data %>% filter(category == "all" & n < 5) %>% pull(biotype)
# Remove biotypes which are present in less than five transcripts
data <-data %>% filter(!biotype %in%  remove)
data$category <- factor(data$category, levels=c("nested", "host", "all"))

#Plot the data
data %>%  filter(biotype %in% common_types ) %>%
  ggplot(aes(category, prop, fill = biotype)) +
  geom_col(position = position_stack(reverse = T), color = "white") +
  coord_flip() + theme_classic() +
  scale_fill_viridis(option="magma", discrete = T) +
  ylab("Proportion total transcripts (%)") + xlab("")

#Biotypes specific to mouse
specific_types = setdiff(unique(all$V7), common_types)

#Calculate the percentage for all transcripts for these biotypes
a <- all %>% dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "all") %>%
    dplyr::rename(biotype=V7)
sum(a$n)

#Calculate the percentage for host transcripts for these biotypes
b <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Host) %>%dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "host") %>%
    dplyr::rename(biotype=V7)
sum(b$n)

#Calculate the percentage for nested transcripts for these biotypes
c <- all %>% filter(V5 %in% hn_gencode$SYMBOL_Nested) %>% dplyr::count(V7) %>% filter(V7 %in% specific_types) %>%
    mutate(prop=n/sum(n)*100) %>%
    mutate(category = "nested") %>%
    dplyr::rename(biotype=V7)
sum(c$n)

#Pool the datasets
data <- rbind(a,b,c)
data$category <- factor(data$category, levels=c("nested", "host", "all"))

#Plot the data
data %>%  ggplot(aes(category, prop, fill = biotype)) +
    geom_col(position = position_stack(reverse = T), color = "white") +
    coord_flip() + theme_classic() +
    scale_fill_viridis(option="cividis", discrete = T) +
    ylab("Proportion total transcripts (%)") + xlab("")

#Heatmap to see the biotype correspondence between Host/nested genes

#Add an ID to the host/nested gene pairs
hn_gencode$pairid <- paste(seq.int(nrow(hn_gencode)),"_ID",sep = "")

#Save the table
write.table(hn_gencode, file = "Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", sep = '\t')

    #Use the orientation already defined as an annotation
    hn_gencode$Strandness = as.factor(hn_gencode$Strandness)
    Anno = hn_gencode[,c(17,16)]
    row.names(Anno) = Anno$pairid

    #Make a table with the biotypes for the heatmap
    Gencode_biotype = hn_gencode[,c(14,7,17)]
    row.names(Gencode_biotype)=Gencode_biotype$pairid
    Gencode_biotype = Gencode_biotype[,-3]

    #Make the biotypes as factor to have one color per biotype
    Gencode_biotype$biotype_Host = as.factor(Gencode_biotype$biotype_Host)
    Gencode_biotype$biotype_Nested = as.factor(Gencode_biotype$biotype_Nested)

    #Transform the table as a matrix for the heatmap function
    Gencode_biotype_h = as.matrix(Gencode_biotype)

    #Create the objects needed for the annotation of the heatmap
    row_ha = rowAnnotation(orientation = Anno$Strandness, col = list(orientation= c("same"="grey50", "opposite"="coral1")))
    col =polychrome(n = nlevels(Gencode_biotype$biotype_Host))
    names(col) = levels(Gencode_biotype$biotype_Host)
    biotype_anno = rowAnnotation(biotype = Gencode_biotype$biotype_Host, col = list(biotype = col))

    #Plot the heatmap
    Heatmap(Gencode_biotype_h, name = "Biotype", col = inferno(length(test)),
            column_title = "Heatmap_Biotype",
            row_split = Gencode_biotype_h[,1:2],
            row_title = NULL,
            show_row_names = FALSE,
            row_gap = unit(0, "mm"),
            right_annotation = row_ha,
            left_annotation = biotype_anno)
