# Comparison of the size of host and nested genes with genes in the genome #

library("dplyr")
library(ggplot2)
library(viridis)

##################
#### HUMAN #######
##################

#Open the host/nested genes table and the reference with all the transcripts used to define the pairs
hn_gencode <- read.delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt", row.names = 1)
all <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")

#To calculate the length of each gene in the reference, take the longest transcript as that what we selected for the pairs too
all$Transcript_length = all$V3 - all$V2
all = arrange(all, -all$Transcript_length)
all = all[!duplicated(all$V5),]
write.table(all,file = "GENCODEV36lift37_longest_transcript.txt", sep = "\t", col.names = T, row.names = F)

#Calculate the gene size for the host and nested genes
hn_gencode$Host_size = hn_gencode$geneEnd_Host - hn_gencode$geneStart_Host
hn_gencode$Nested_size = hn_gencode$geneEnd_Nested - hn_gencode$geneStart_Nested

#Pool the data to make a table for the graph
#Reference
all$annotation = "All"
all = all[,c(8,9)]
colnames(all)[1] = "Size"
#Host
Host_size = hn_gencode[,c(12,18)]
Host_size$annotation = "Host"
Host_size = Host_size[,c(2,3)]
colnames(Host_size)[1] = "Size"
#Nested genes
Nested_size = hn_gencode[,c(5,19)]
Nested_size$annotation = "Nested"
Nested_size = Nested_size[,c(2,3)]
colnames(Nested_size)[1] = "Size"

#Pool all of them
data = rbind(all,Nested_size,Host_size)

#Make the plot with the distribution of the size
ggplot(data, aes(x=Size, fill=annotation)) +
  geom_density(alpha=.5)+
  scale_x_continuous(trans='log10')

#Statistics to test for the difference between the size distributions
#Comparison of the distribution using a KS test as the variance is different

ks.test(log10(all$Size),log10(Host_size$Size))

#Results
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  log10(all$Size) and log10(Host_size$Size)
#D = 0.51141, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(log10(all$Size),log10(Nested_size$Size))

#Results
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  all$Size and Nested_size$Size
#D = 0.337, p-value < 2.2e-16
#alternative hypothesis: two-sided

#Make a plot to test the correlation between the size of the host and the nested genes
ggplot(hn_gencode, aes(x=Nested_size, y=Host_size)) +
  geom_point(size=0.7)+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_viridis(alpha = 0.3,option = "C")+
  labs(title="Correlation_host_nested_size_Human",x="Nested_gene_size", y = "Host_gene_size")

#Correlation analysis - Test of the data normality
shapiro.test(sample(log10(hn_gencode$Nested_size),size=5000))
#Results
#Shapiro-Wilk normality test
#data:  sample(log10(hn_gencode$Nested_size), size = 5000)
#W = 0.94877, p-value < 2.2e-16 --> Data are not normally distributed

shapiro.test(sample(log10(hn_gencode$Host_size),size=5000))
#Results
#Shapiro-Wilk normality test
#data:  sample(log10(hn_gencode$Host_size), size = 5000)
#W = 0.99318, p-value = 9.853e-15 --> Data are not normally distributed

#Test for correlation using the non-parametric spearman correlation
corr <- cor.test(x=log10(hn_gencode$Nested_size), y=log10(hn_gencode$Host_size), method = 'spearman')
#Spearman's rank correlation rho
#data:  hn_gencode$Nested_size and hn_gencode$Host_size
#S = 3.0942e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.1719101

##################
#### MOUSE #######
##################

#Open the host/nested genes table and the reference with all the transcripts used to define the pairs
hn_gencode = read.delim("Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", row.names = 1)
all <- read.delim("wgEncodeGencodeCompVM25_fileforhost_filtered.bed", header = F)

#To calculate the length of each gene in the reference, take the longest transcript as that what we selected for the pairs too
all$Transcript_length = all$V3 - all$V2
all = arrange(all, -all$Transcript_length)
all = all[!duplicated(all$V5),]
write.table(all,file = "wgEncodeGencodeCompVM25_longest_transcript.txt", sep = "\t", col.names = T, row.names = F)

#Calculate the gene size for the host and nested genes
hn_gencode$Host_size = hn_gencode$geneEnd_Host - hn_gencode$geneStart_Host
hn_gencode$Nested_size = hn_gencode$geneEnd_Nested - hn_gencode$geneStart_Nested

#Pool the data to make a table for the graph
#Reference
all$annotation = "All"
all = all[,c(8,9)]
colnames(all)[1] = "Size"
#Host
Host_size = hn_gencode[,c(12,18)]
Host_size$annotation = "Host"
Host_size = Host_size[,c(2,3)]
colnames(Host_size)[1] = "Size"
#Nested genes
Nested_size = hn_gencode[,c(5,19)]
Nested_size$annotation = "Nested"
Nested_size = Nested_size[,c(2,3)]
colnames(Nested_size)[1] = "Size"

#Pool all of them
data = rbind(all,Nested_size,Host_size)

#Make the plot with the distribution of the size
ggplot(data, aes(x=Size, fill=annotation)) +
  geom_density(alpha=.5)+
  scale_x_continuous(trans='log10')

#Statistics to test for the difference between the size distributions
#Comparison of the distribution using a KS test as the variance is different

ks.test(log10(all$Size),log10(Host_size$Size))

#Results
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  log10(all$Size) and log10(Host_size$Size)
#D = 0.49214, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(log10(all$Size),log10(Nested_size$Size))

#Results
#Asymptotic two-sample Kolmogorov-Smirnov test
#data:  log10(all$Size) and log10(Nested_size$Size)
#D = 0.39145, p-value < 2.2e-16
#alternative hypothesis: two-sided

#Make a plot to test the correlation between the size of the host and the nested genes
ggplot(hn_gencode, aes(x=Nested_size, y=Host_size)) +
  geom_point(size=0.7)+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_viridis(alpha = 0.3,option = "C")+
  labs(title="Correlation_host_nested_size_Mouse",x="Nested_gene_size", y = "Host_gene_size")

#Correlation analysis - Test of the data normality
shapiro.test(sample(log10(hn_gencode$Nested_size),size=5000))
#Results
#Shapiro-Wilk normality test
#data:  sample(log10(hn_gencode$Nested_size), size = 5000)
#W = 0.88579, p-value < 2.2e-16 --> Data are not normally distributed

shapiro.test(sample(log10(hn_gencode$Host_size),size=5000))
#Results
#Shapiro-Wilk normality test
#data:  sample(log10(hn_gencode$Host_size), size = 5000)
#W = 0.98535, p-value < 2.2e-16 --> Data are not normally distributed

#Test for correlation using the non-parametric spearman correlation
corr <- cor.test(x=log10(hn_gencode$Nested_size), y=log10(hn_gencode$Host_size), method = 'spearman')
#Spearman's rank correlation rho
#data:  log10(hn_gencode$Nested_size) and log10(hn_gencode$Host_size)
#S = 5.3888e+10, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho
#0.2516947
