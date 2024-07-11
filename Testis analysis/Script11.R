#Analysis of the cluster of pairs with co-expression in testis using single cell dataset #

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(tidyverse)
library(cluster)
library(factoextra)

#Dataset used was from https://www.sciencedirect.com/science/article/pii/S2666379121002536
#Accession number is https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153947
# Downloaded tables contains the integrated, normalized counts per cell for 3 normal samples "GSE153947_Normal_integrated_data.tsv.gz" and the metadata information ("GSE153947_Cell_metadata.tsv.gz")

Meta = read.delim("GSE153947_Cell_metadata.tsv.gz")
data = read.delim("GSE153947_Normal_integrated_data.tsv.gz")

#Keep only the cells for which we have metadata information and expression data
keep = rownames(Meta) %in% colnames(data)
data = data[,keep]
Meta = Meta[keep,]

#Open the table with the pairs of genes with co-expression in testis
hn_testis_pairs = read.csv("testis_hn_genes_logged.csv")


#Make a list of the names of the host and nested genes from these pairs
all_testis = c(hn_testis_pairs$SYMBOL_Host,hn_testis_pairs$SYMBOL_Nested)

#Keep only the data for these genes
keep_2 = row.names(data) %in% all_testis
Testis_hn_data = data[keep_2,]
write.table(Testis_hn_data, file = "Testis_hn_data.txt", sep = "\t", col.names = T, row.names =T)

#Recover the data for the host and the nested genes separalty to find the pairs with data for both genes
Testis_host_data = Testis_hn_data[row.names(Testis_hn_data) %in% hn_testis_pairs$SYMBOL_Host,]
Testis_nested_data = Testis_hn_data[row.names(Testis_hn_data)%in% hn_testis_pairs$SYMBOL_Nested,]

#Add the pair ID to the table and use it as a row names
#Host
Testis_host_data = merge.data.frame(Testis_host_data,hn_testis_pairs[,c(12,17)], by.x =0, by.y = 1)
row.names(Testis_host_data) = Testis_host_data$pairid
#Remove the fisrt and last colum which are now gene names and ID respectively
Testis_host_data= Testis_host_data[,c(-1,-15534)]

#Nested
Testis_nested_data = merge.data.frame(Testis_nested_data,hn_testis_pairs[,c(5,17)], by.x =0, by.y = 1)
row.names(Testis_nested_data) = Testis_nested_data$pairid
#Remove the fisrt and last colum which are now gene names and ID respectively
Testis_nested_data= Testis_nested_data[,c(-1,-15534)]

#Keep the pairs where both the nested gene and the host have expression data
#Host with data in nested table
Keep_3 = rownames(Testis_host_data) %in% rownames(Testis_nested_data)
Testis_host_data = Testis_host_data[Keep_3,]
#Nested genes with data in the host table
Keep_4 = rownames(Testis_nested_data) %in% rownames(Testis_host_data)
Testis_nested_data = Testis_nested_data[Keep_4,]
#64 pairs with data for both the host and the nested gene

#Put the sample ID as rownames of the Metadata table
Meta$sample  = row.names(Meta)
#Keep only the tissue in a smaller table
tissue_meta <- Meta %>% select(sample = "sample", tissue = "Cluster_identity")

#Correlation analysis between the expression in the single cells of the two genes of the 64 pairs
#Select the IDs of the 64 pairs
pair_ids = rownames(Testis_nested_data)
#Make empty objects to collect the data
datalist = list()
plot_list = list()
#Make a metadata table for the pairs with the pair ID and the gene symbol of the Host and the Nested gene
pair_meta = hn_testis_pairs[,c(17,5,12)]
colnames(pair_meta) = c("pair_id", "SYMBOL_Nested", "SYMBOL_Host")

#Conduct the Spearman correlation analysis and make a graph for each pairs which are then stored as lists
for (i in pair_ids) {
  group1 <- Testis_nested_data[i,]
  group2 <- Testis_host_data[i,]
  data <- data.frame(nest=as.numeric(group1), host=as.numeric(group2), sample = colnames(Testis_nested_data), pair_id = i)
  data <- merge(data, tissue_meta, by = "sample")
  data <-merge(data, pair_meta, by = "pair_id")
  p <- ggscatter(data, x = "nest", y = "host", fill = "tissue", alpha = 0.5 ,color = "tissue", add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) +
    stat_cor(method = "spearman", label.y.npc="top", label.x.npc = "left") + ggtitle(paste0(i,"\n","Nested Gene = ",data$SYMBOL_Nested,"\n","Host Gene = ",data$SYMBOL_Host))+font("legend.text", size =4)
  plot_list[[i]] <-  p
  data <- data %>% group_by(pair_id) %>% summarize(cor=cor(nest, host, method = "spearman"))
  datalist[[i]] <- data
  print(i)
}

#Make a table with the correlation coefficient
cor_data <-  do.call(rbind, datalist)

#Plot with the distribution of the Correlation
cor_data %>% ggplot(aes(cor)) + geom_histogram(color="black", fill="grey", binwidth = 0.05) + xlim(-0.5, 0.7)

#Print all the scatter plot of the pearson corralation into a unique PDF document
pdf("hg_scatter_plots_testis_sc.pdf")
for (i in pair_ids) {
  print(plot_list[[i]])
}
dev.off()

#To plot individual plots use with i replaced by a pair ID
print(plot_list[[i]])

#Segregation of the data to be able to determine is genes are co-expressed or mutually exclusive
#For this we use as a threshold of expression a quarter of the maximum of the expression for each genes
#Data are transpose to be able to select column
transposed_host <- as.data.frame(t(Testis_host_data))
transposed_nested <- as.data.frame(t(Testis_nested_data))
#Make sure that the order of the column is the same in both tables
transposed_host<-transposed_host[names(transposed_nested)]
transposed_nested<-transposed_nested[names(transposed_host)]

#Create an empty vector to store the co-expression data
v <- c()

# Arrange meta data for barcode merging
meta_filtered <- Meta %>% select(Cluster_identity)
meta_filtered$barcode <- rownames(meta_filtered)

#Create two empty lists to store the results
datalist = list()
datalist2 = list()

#Loop to repeat the analysis on all pairs
for (i in colnames(transposed_host)) {
  data_host <- transposed_host %>% select(e=all_of(i))
  #Calculate the quarter of the maximum expression and select the cells above this treshold
  quartermax <- max(as.numeric(data_host$e), na.rm = T)/4
  quarter_host <-  data_host$e > quartermax
  data_host$above_host <- quarter_host
  data_host$barcode <- rownames(data_host)
  #Count the cells with expression
  count_over_host <- data_host %>% count(above_host) %>% pull(n)
  count_over_host <- count_over_host[2]

  data_nested <- transposed_nested %>% select(e=all_of(i))
  #Calculate the quarter of the maximum expression and select the cells above this treshold
  quartermax <- max(as.numeric(data_nested$e), na.rm = T)/4
  quarter_nested <-  data_nested$e > quartermax
  data_nested$above_nested <- quarter_nested
  data_nested$barcode <- rownames(data_nested)
  #Count the cells with expression
  count_over_nested <- data_nested %>% count(above_nested) %>% pull(n)
  count_over_nested <- count_over_nested[2]

  #Recover the data for both the host and the nested gene
  x = merge(data_nested, data_host, by = "barcode")
  #Count the cells with co-expression
  table <- as.data.frame(table(x$above_nested+x$above_host))
  prop <- table[3,2]
  percentage <- prop/length(quarter_host+quarter_nested)
  v <- c(v, percentage)

  #Count the cells per categories et per cell type (co-expression, expression of the host only or the nested only and no expression for both)
  x <- x %>% mutate(co_expressed = above_host+above_nested) %>% select(barcode, co_expressed, above_host, above_nested) %>% mutate(pairid = i) %>% mutate(co_expressed=factor(co_expressed, levels =c("0", "1","2")))
  x2 <- merge(x,meta_filtered, by="barcode")
  x2 <- x2 %>% mutate(above_host = factor(above_host,levels = c("TRUE","FALSE"))) %>%
    mutate(above_nested = factor(above_nested,levels = c("TRUE","FALSE")))
  x <- merge(x,meta_filtered, by="barcode")
  x <- x %>% count(pairid, Cluster_identity,co_expressed, .drop = FALSE)
  x2 <- x2 %>% count(pairid,Cluster_identity,above_host, above_nested, .drop = FALSE)

  datalist[[i]] <- x
  datalist2[[i]] <- x2
  print(i)

}

#Make a dataframe with the percentage of co-expression per pair globally
proportion <- data.frame(pairid = colnames(transposed_host),prop_above_quarter = v)

#Table with number of cells with co-expression, mono-expression and no-expression per pair and cell type
big_data = do.call(rbind, datalist)

#Table with number of cells with co-expression, expression of host only, expression of nested only and no-expression per pair and cell type
big_data2 = do.call(rbind, datalist2)

#Save the tables
write.csv(big_data, file = "big_data_new.csv")
write.csv(big_data2, file = "big_data2_new.csv")

#Similar loop and strategy use to calculate the percentages of cells with host expression and the percentage of cells with nested expression
datalist = list()

for (i in colnames(transposed_host)) {
  data_host <- transposed_host %>% select(e=i)
  #Calculate the quarter of the maximum expression and select the cells above this threshold
  quartermax <- max(as.numeric(data_host$e), na.rm = T)/4
  quarter_host <-  data_host$e > quartermax
  data_host$above_host <- quarter_host
  data_host$barcode <- rownames(data_host)
  #Count the cells with expression
  count_over_host <- data_host %>% count(above_host) %>% pull(n)
  count_over_host <- count_over_host[2]

  data_nested <- transposed_nested %>% select(e=i)
  #Calculate the quarter of the maximum expression and select the cells above this threshold
  quartermax <- max(as.numeric(data_nested$e), na.rm = T)/4
  quarter_nested <-  data_nested$e > quartermax
  data_nested$above_nested <- quarter_nested
  data_nested$barcode <- rownames(data_nested)
  #Count the cells with expression
  count_over_nested <- data_nested %>% count(above_nested) %>% pull(n)
  count_over_nested <- count_over_nested[2]

  #Count the cells per categories et per cell type (co-expression, expression of the host only or the nested only and no expression for both)
  x <- merge(data_nested, data_host, by = "barcode")
  x <- x %>% mutate(co_expressed = above_host+above_nested) %>% select(barcode, co_expressed, above_host, above_nested) %>% mutate(pairid = i) %>% mutate(co_expressed=factor(co_expressed, levels =c("0", "1","2")))
  #Select the cells with expression of the host
  x1 <- x %>% filter(above_host == TRUE)%>% mutate(only = "host_expression")
  #Select the cells with expression of the nested gene
  x2 <-  x %>% filter(above_nested == TRUE)%>% mutate(only = "nested_expression")
  #Count the number of cells per category
  x5 = rbind(x1,x2)
  x6 <- x5 %>% count(pairid,only,.drop = FALSE)
  #Calculate the percentage
  x6$n = (x6$n/15532)*100
  #Make a table out of it
  x6 = as.data.frame(t(x6))
  colnames(x6) = x6[2,]
  x6 = x6[-1:-2,]
  rownames(x6)[1]=i

  datalist[[i]] <- x6
  print(i)

}

#Recover the table with the percentage of expression in all cells
Host_nested_expression = do.call(rbind, datalist)
Host_nested_expression$host_expression = as.numeric(Host_nested_expression$host_expression)
Host_nested_expression$nested_expression = as.numeric(Host_nested_expression$nested_expression)

#Save the table
write.table(Host_nested_expression, file = "Percentage_cells_Host_nested_expression.txt", col.names = T, row.names = T, sep = "\t")


#Similar loop and strategy use to calculate the percentages of cells with co-expression, expression of the host only, expression of the nested gene only and of none of them per pair in all cells
datalist = list()
for (i in colnames(transposed_host)) {
  data_host <- transposed_host %>% select(e=i)
  #Calculate the quarter of the maximum expression and select the cells above this threshold
  quartermax <- max(as.numeric(data_host$e), na.rm = T)/4
  quarter_host <-  data_host$e > quartermax
  data_host$above_host <- quarter_host
  data_host$barcode <- rownames(data_host)
  #Count the cells with expression
  count_over_host <- data_host %>% count(above_host) %>% pull(n)
  count_over_host <- count_over_host[2]

  data_nested <- transposed_nested %>% select(e=i)
  #Calculate the quarter of the maximum expression and select the cells above this treshold
  quartermax <- max(as.numeric(data_nested$e), na.rm = T)/4
  quarter_nested <-  data_nested$e > quartermax
  data_nested$above_nested <- quarter_nested
  data_nested$barcode <- rownames(data_nested)
  #Count the cells with expression
  count_over_nested <- data_nested %>% count(above_nested) %>% pull(n)
  count_over_nested <- count_over_nested[2]

  #Count the cells per categories et per cell type (co-expression, expression of the host only or the nested only and no expression for both)
  x <- merge(data_nested, data_host, by = "barcode")
  x <- x %>% mutate(co_expressed = above_host+above_nested) %>% select(barcode, co_expressed, above_host, above_nested) %>% mutate(pairid = i) %>% mutate(co_expressed=factor(co_expressed, levels =c("0", "1","2")))
  #Select the cells with expression of both and add annotation
  x1 <- x %>% filter(above_host == TRUE & above_nested == TRUE)%>% mutate(only = "co_expression")
  #Select the cells with expression of the host only and add annotation
  x2 <-  x %>% filter(above_host == TRUE & above_nested == FALSE)%>% mutate(only = "host_only")
  #Select the cells with expression of the nested gene only and add annotation
  x3 <- x %>% filter(above_host == FALSE & above_nested == TRUE)%>%  mutate(only = "nested_only")
  #Select the cells with no expression of the nested and the host genes and add annotation
  x4 <- x %>% filter(above_host == FALSE & above_nested == FALSE) %>% mutate(only = "none")
  #Merge all the table together
  x5 = rbind(x1,x2,x3,x4) %>% mutate(only = factor(only, levels =c("co_expression", "host_only","nested_only","none")))
  #Count the number of cells per categories
  x6 <- x5 %>% count(pairid,only,.drop = FALSE)
  #Calculate the percentage
  x6$n = (x6$n/15532)*100
  x6 = as.data.frame(t(x6))
  colnames(x6) = x6[2,]
  x6 = x6[-1:-2,]
  rownames(x6)[1]=i

  datalist[[i]] <- x6
  print(i)

}

global_proportion = do.call(rbind, datalist)

#Save the table
write.table(global_proportion, file = "Percentage_all_cells.txt", col.names = T, row.names = T, sep = "\t")

#Make a plot for individual pair separated per cell type or not to show the strategy of threshold of expression
#Here the pair ID can be changed to make a plot foe a different pair - Example 2151_ID
group1 <- Testis_nested_data["2151_ID",]
group2 <- Testis_host_data["2151_ID",]
#Recover the data and organise them for ggplot
data <- data.frame(nest=as.numeric(group1), host=as.numeric(group2), sample = colnames(Testis_nested_data), pair_id = "2151_ID")
data <- merge(data, tissue_meta, by = "sample")
data <-merge(data, pair_meta, by = "pair_id")
quartermax_host <- max(as.numeric(data$host), na.rm = T)/4
quartermax_nested <- max(as.numeric(data$nest), na.rm = T)/4

p = ggplot(data, aes(nest,host)) +
  geom_rect(aes(xmin = quartermax_nested, xmax = Inf, ymin = quartermax_host , ymax = Inf),fill = "pink", alpha = 0.02)+
  geom_rect(aes(xmin = -Inf, xmax = quartermax_nested, ymin = quartermax_host , ymax = Inf),fill = "orange", alpha = 0.005)+
  geom_rect(aes(xmin = quartermax_nested, xmax = Inf, ymin = -Inf , ymax = quartermax_host),fill = "lightblue", alpha = 0.05)+
  geom_point(aes(colour = factor(tissue)),size = 1.8, alpha = 0.7) +
  geom_hline(yintercept=quartermax_host, linetype="dashed", color = "grey",  linewidth =1.5) +
  geom_vline(xintercept = quartermax_nested, linetype="dashed", color = "grey", linewidth=1.5) +
  theme(legend.position="none")+
  coord_fixed(ratio = 1) +
  xlab("Nested") + ylab("Host") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#Plot for all cells
p
#Plot scattered per cell type
p + facet_wrap(~ tissue, ncol=3)

#Because we want to study co-expresion we need to use the pairs were both genes are at least expressed 
#We decided to use 1% of expression in all the cells as a threshold
#Select the pair with at least 1% of cells with expression of each genes
Host_nested_expression_filtered <- Host_nested_expression %>% filter(host_expression > 1 & nested_expression > 1)
#Recover the name of the selected pairs
pairs_filtered = row.names(Host_nested_expression_filtered)

#Keep only the data for the pairs which were filtered
#Keep only the pairs which were filtered
big_data = big_data[big_data$pairid %in% pairs_filtered,]
big_data2 = big_data2[big_data2$pairid %in% pairs_filtered,]

#Calculate the percentage of cells per category
big_data$co_expressed <- as.factor(big_data$co_expressed)
#Count the number of cells per cell type
cell_type_count <- meta_filtered %>% count(Cluster_identity) %>% select(Cluster_identity,total = n)
#Add this to the previous table
new_big_data <- merge(big_data, cell_type_count, by = "Cluster_identity")
#Calculate percentage
new_big_data <- new_big_data %>% mutate(prop_total=n/total)

# Order the cell type to match the trajectory during spermatogenesis and then add somatic cells at the beginning
unique(new_big_data$Cluster_identity)
order1 <- as.character(unique(new_big_data$Cluster_identity))[c(14,1,7,15,11,2,10,3,6)]
order2 <- setdiff(as.character(unique(new_big_data$Cluster_identity)),order1)
order1 <- factor(order1, levels=unique(order1))
order2<- factor(order2, levels=unique(order2))
order <- lvls_union(list(order2, order1))

#Plot the profile of co-expression for all the pairs in somatic cells and during spermatogenesis
new_big_data %>% filter(co_expressed == 2) %>% filter(Cluster_identity %in% order) %>%
  ggplot(aes(x = factor(Cluster_identity, level = order), y=prop_total, group=pairid)) +
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + theme_classic() + ylab("Cells co-expressing host/nested gene pair (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none") + xlab("") +
  geom_vline(xintercept = 6.5, color = "darkgrey", size = 2, alpha = 0.5, linetype = "dashed")

#Group the pairs by expression profile to try to find patterns
#Make a new table with only the co-expression and few selected columns
explore <- new_big_data %>% filter(co_expressed == 2) %>% filter(Cluster_identity %in% order)
wide <- explore %>% select(Cluster_identity, pairid, prop_total)
#kmean clustering to group profiles
matrix <- reshape2::acast(wide, Cluster_identity~pairid, value.var="prop_total")
matrix <- t(matrix)

#Plot the results to find the optimal number of cluster
fviz_nbclust(matrix, kmeans, method = "wss", k.max = 20)

# Bend in the knee appears to be 3 clusters - Select this
matrix <- matrix[ , which(apply(matrix, 2, var) != 0)]
k2 <- kmeans(matrix, centers = 3, nstart = 25)

#Recover the results of kmean clustering
clusters <- as.data.frame(k2$cluster)
clusters <- clusters %>% select(cluster=`k2$cluster`)
clusters$pairid <- rownames(clusters)

# Merge back with explore data to be able to plot it
explore_clusters <- merge(explore,clusters,by="pairid")
explore_clusters <- explore_clusters %>% mutate(cluster = as.factor(cluster))

#Make the graph with the profile of co-expression in the 3 different groups
explore_clusters %>%
  ggplot(aes(x = factor(Cluster_identity, level = order), y=prop_total,color = cluster, group=pairid)) +
  geom_point() + geom_line() + facet_grid(.~cluster) + theme_bw()+
  scale_colour_brewer(palette = "Dark2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") +
  ylab("Cells co-expressing host/nested gene pair (%)")

#Count the number of pairs per group
count(explore_clusters, cluster, pairid) %>% count(cluster)
##cluster  n
#1       1 30
#2       2 1
#3       3 3

#Use a similar approach for the profile of exclusive expression of the host or the nested gene
#This time need to use big_data2
#Add the cell counts to calculate percentages

new_big_data2 <- merge(big_data2, cell_type_count, by = "Cluster_identity")
w <- new_big_data2 %>% filter(above_host == TRUE & above_nested == TRUE) %>% mutate(proportion_only = n/total) %>% mutate(only = "co_expression")
x <- new_big_data2 %>% filter(above_host == TRUE & above_nested == FALSE) %>% mutate(proportion_only = n/total) %>% mutate(only = "host_only")
y <- new_big_data2 %>% filter(above_host == FALSE & above_nested == TRUE) %>% mutate(proportion_only = n/total) %>% mutate(only = "nested_only")
z <- new_big_data2 %>% filter(above_host == FALSE & above_nested == FALSE) %>% mutate(proportion_only = n/total) %>% mutate(only = "none")
All_expression_data <- rbind(w,x,y,z)

# Spread table and count proportion per cell identity and pairid, replace NA values with 0 
percent = All_expression_data %>% select(Cluster_identity, pairid, proportion_only, only) %>%
  mutate(only = as.factor(only)) %>%
  pivot_wider(names_from = only, values_from = proportion_only) %>% replace(is.na(.), 0)

percent$Cluster_identity <- factor(percent$Cluster_identity , levels = c("Endothelial cells","Fibrotic peritubular myoid cells","Leydig cells",
                                                                         "Macrophages","PMCs","Sertoli cells",
                                                                         "Undifferentiated spermatogonia","Diff.spermatogonia/Preleptotene",
                                                                         "Leptotene","Zygotene","Pachytene",
                                                                         "Diplotene","Meiotic divisions",
                                                                         "Early spermatids","Late spermatids"))

#Do the kmean clustering analysis for the host only and nested only
#Host
#Plot the profile of host only expression for all the pairs in somatic cells and during spermatogenesis
percent %>%
  ggplot(aes(x = Cluster_identity , y=host_only, group=pairid)) +
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + theme_classic() + ylab("Cells expressing exclusively host gene (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none") + xlab("") +
  geom_vline(xintercept = 6.5, color = "darkgrey", size = 2, alpha = 0.5, linetype = "dashed")+ ylim(c(0,1))

#kmean clustering to group profiles
matrix <- reshape2::acast(percent, Cluster_identity~pairid, value.var="host_only")
matrix <- t(matrix)

#Plot the results to find the optimal number of cluster
fviz_nbclust(matrix, kmeans, method = "wss", k.max = 20)

# Bend in the knee appears to be 4 clusters - Select this
k2 <- kmeans(matrix, centers = 4, nstart = 25)

#Recover the results of kmean clustering
clusters <- as.data.frame(k2$cluster)
clusters <- clusters %>% select(cluster=`k2$cluster`)
clusters$pairid <- rownames(clusters)

# Merge back with percent data to be able to plot it
explore_clusters <- merge(percent,clusters,by="pairid")
explore_clusters <- explore_clusters %>% mutate(cluster = as.factor(cluster))

#Make the graph with the profile of co-expression in the 3 different groups
explore_clusters %>%
  ggplot(aes(x = Cluster_identity, y=host_only,color = cluster, group=pairid)) +
  geom_point() + geom_line() + facet_grid(.~cluster) + theme_bw()+
  scale_colour_brewer(palette = "Dark2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") +
  ylab("Cells expressing exclusively the host (%)")+ ylim(c(0,1))

#Count the number of pairs per group
count(explore_clusters, cluster, pairid) %>% count(cluster)
#cluster  n
#1         1 2
#2         2 9
#3         3 1
#4         4 22

#Nested

#Plot the profile of nested only expression for all the pairs in somatic cells and during spermatogenesis
percent %>%
  ggplot(aes(x = Cluster_identity , y=nested_only, group=pairid)) +
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + theme_classic() + ylab("Cells expressing exclusively nested gene (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none") + xlab("") +
  geom_vline(xintercept = 6.5, color = "darkgrey", size = 2, alpha = 0.5, linetype = "dashed")+ ylim(c(0,1))

#kmean clustering to group profiles
matrix <- reshape2::acast(percent, Cluster_identity~pairid, value.var="nested_only")
matrix <- t(matrix)

#Plot the results to find the optimal number of cluster
fviz_nbclust(matrix, kmeans, method = "wss", k.max = 20)

# Bend in the knee appears to be 4 clusters - Select this
k2 <- kmeans(matrix, centers = 4, nstart = 25)

#Recover the results of kmean clustering
clusters <- as.data.frame(k2$cluster)
clusters <- clusters %>% select(cluster=`k2$cluster`)
clusters$pairid <- rownames(clusters)

# Merge back with percent data to be able to plot it
explore_clusters <- merge(percent,clusters,by="pairid")
explore_clusters <- explore_clusters %>% mutate(cluster = as.factor(cluster))

#Make the graph with the profile of co-expression in the 3 different groups
explore_clusters %>%
  ggplot(aes(x = Cluster_identity, y=nested_only,color = cluster, group=pairid)) +
  geom_point() + geom_line() + facet_grid(.~cluster) + theme_bw()+
  scale_colour_brewer(palette = "Dark2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") +
  ylab("Cells expressing exclusively the host (%)")+ ylim(c(0,1))

#Count the number of pairs per group
count(explore_clusters, cluster, pairid) %>% count(cluster)
#cluster  n
#1         1 2
#2         2 1
#3         3 3
#4         4 28

#Individual plots with the percentage of cells with co-expression, host only expression and nested only expression
#Example here for the pair with mutually exclusive expression but can be done for any of the pairs with expression in testis

#Prepare the table for the graph
graph = percent %>% select(-none) %>%
  pivot_longer(cols=c("host_only","nested_only","co_expression"),names_to = "type", values_to = "value")

#Select the pair with the pairid
ID_2510_table_host_nested = graph[graph$pairid == "2510_ID",]

#Make the graph
ggplot(data=ID_2510_table_host_nested, aes(x=Cluster_identity, y=value, group=type)) +
  geom_line(aes(color=type))+
  geom_point(aes(color=type))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


