# Analysis of the number of isoforms per genes to see if there is a difference for the host and nested genes #

library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
library(dplyr)
library(twosamples)

##################
#### HUMAN #######
##################

#Open the table with the host/nested gene pairs data
hn_gencode = read.delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt", row.names = 1)
#Open the table with information about all the transcripts used in the analysis
all <- read.table("wgEncodeGencodeComprehensiveV36lift37_fileforhost_filtered.bed")

#Count the number of isoforms per gene
Isoforms = all %>% count(all[,5])
#Add the annotation column for later
Isoforms$list = "All"

#Save the table
write.table(all, file = "All_HUMAN_isoform_number.txt", col.names = T, sep = "\t")

#Select only the genes which are hosts
Host = merge.data.frame(hn_gencode[,c(17,12,14)],Isoforms,by.x = 2, by.y = 1, all.x = T)
#Remove the duplication for host genes involved in multiple pairs
Host = Host[!duplicated(Host$SYMBOL_Host),]
#Remove the pair id and the biotype
Host = Host[,-2:-3]
#Change the annotation for this list
Host$list = "Host"


#select only the genes which are nested 
Nested = merge.data.frame(hn_gencode[,c(17,5,7)],Isoforms,by.x = 2, by.y = 1, all.x = T)
#Remove the duplication for host genes involved in multiple pairs
Nested = Nested[!duplicated(Nested$SYMBOL_Nested),]
#Remove the pair id and the biotype
Nested = Nested[,-2:-3]
#Change the annotation for this list
Nested$list = "Nested"

#Make a big table for all of them
df = Host
colnames(Nested) = colnames(df)
colnames(Isoforms) = colnames(df)
df = rbind(df,Isoforms)
df = rbind(df,Nested)

#Statistical test to test if the cumulative distribution of isoforms number are different
#A two-sample test based on the DTS test statistic (dts_stat). This is the recommended two-sample test in this package because of its power.
#The DTS statistic is the reweighted integral of the distance between the two ECDFs.
# Between All and Host
t_all_host = two_sample(Isoforms$n,Host$n)
#Test Stat    P-Value
#4251.10443    0.00025

t_all_nested = two_sample(Isoforms$n,Nested$n)
#Test Stat   P-Value
#2860.0933    0.0065

#Make the graph
ggplot(df, aes(x = n , color=list)) +
  scale_x_continuous(trans = "log10")+
  stat_ecdf(geom = "step")+
  annotate(geom="text", x=30, y=0.5, label=paste("p_value",t_all_host[2],sep = " : ") ,
           color="darkgreen")+
  annotate(geom="text", x=30, y=0.45, label=paste("p_value",t_all_nested[2],sep = " : ") ,
           color="darkblue")

#Test of this difference in isoform number is only due to the longer size of the host genes as it was observed before that longer genes tend to have more isoforms
#Regression analysis on the relationship between gene size (taking the longest isoform) and number of isoform if you are a host or not

#Determine gene size using the longest isoform
#Calculate the size
all$size = all[,3] - all[,2]
#Organise with the longest first
all = arrange(all, -all$size)
#Remove the duplication
all = all[!duplicated(all[,5]),]
#Keep only the size and the gene name
all= all[,c(5,8)]

#Add the size to the number of isoform table
Isoform_size = merge.data.frame(Isoforms,all, by.x = 1, by.y = 1)

#Add the information is the gene is host or not
Host = hn_gencode$SYMBOL_Host

for (i in 1:nrow(Isoform_size)){
  if (Isoform_size[i,1] %in% Host) {
    Isoform_size[i,"Host"] = "TRUE"
  } else {
    Isoform_size[i,"Host"] = "FALSE"
  }
  print(i)
}

Isoform_size$Host = as.logical(Isoform_size$Host)
#Make individual tables
Isoform_size_Host = Isoform_size[Isoform_size$Host == TRUE,]
Isoform_size_other = Isoform_size[Isoform_size$Host == FALSE,]

#Plot and regression analysis
ggplot() +
  geom_point(data = Isoform_size_other, aes(x = size, y = n), color = "#fbe4afff") +
  geom_smooth(data = Isoform_size_other, aes(x = size, y = n), method = "lm", se = TRUE, color = "#fbe4afff", level = 0.95) +  geom_point(data = Isoform_size_Host, aes(x = size, y = n),colour = "#1b1248c7") +
  geom_smooth(data = Isoform_size_Host, aes(x = size, y = n), method = "lm", se = TRUE, color = "#1b1248c7", level = 0.95) +
  labs(title = "Regression analysis Human",
       x = "Gene Size",
       y = "Isoforms number") +
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  theme_minimal()

#Calculate the regression model 
#Use the model to predict the number of isoforms for the median size
median(Isoform_size$size)
#10059
lm_other = lm(log(n) ~ log(size),data= Isoform_size_other)
summary(lm_other)
log_isoform_number = log(10059)*0.267480 - 1.465221
exp(log_isoform_number)
#[1] 2.718103
lm_host = lm(log(n) ~ log(size), data = Isoform_size_Host)
summary(lm_host)
log_isoform_number = log(10059)*0.219088 -0.655560
exp(log_isoform_number)
#[1] 3.91026


##################
#### MOUSE #######
##################

#Open the table with the host/nested gene pairs data
hn_gencode = read.delim("Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", row.names = 1)
#Open the table with information about all the transcripts used in the analysis
all <- read.table("wgEncodeGencodeCompVM25_fileforhost_filtered.bed")

#Count the number of isoforms per gene
Isoforms = all %>% count(all[,5])
#Add the annotation column for later
Isoforms$list = "All"

#Save the table
write.table(all, file = "All_mouse_isoform_number.txt", col.names = T, sep = "\t")

#Select only the genes which are hosts
Host = merge.data.frame(hn_gencode[,c(17,12,14)],Isoforms,by.x = 2, by.y = 1, all.x = T)
#Remove the duplication for host genes involved in multiple pairs
Host = Host[!duplicated(Host$SYMBOL_Host),]
#Remove the pair id and the biotype
Host = Host[,-2:-3]
#Change the annotation for this list
Host$list = "Host"


#select only the genes which are nested 
Nested = merge.data.frame(hn_gencode[,c(17,5,7)],Isoforms,by.x = 2, by.y = 1, all.x = T)
#Remove the duplication for host genes involved in multiple pairs
Nested = Nested[!duplicated(Nested$SYMBOL_Nested),]
#Remove the pair id and the biotype
Nested = Nested[,-2:-3]
#Change the annotation for this list
Nested$list = "Nested"

#Make a big table for all of them
df = Host
colnames(Nested) = colnames(df)
colnames(Isoforms) = colnames(df)
df = rbind(df,Isoforms)
df = rbind(df,Nested)

#Statistical test to test if the cumulative distribution of isoforms number are different
#A two-sample test based on the DTS test statistic (dts_stat). This is the recommended two-sample test in this package because of its power.
#The DTS statistic is the reweighted integral of the distance between the two ECDFs.
# Between All and Host
t_all_host = two_sample(Isoforms$n,Host$n)
#Test Stat    P-Value
#1676.45941    0.00025

t_all_nested = two_sample(Isoforms$n,Nested$n)
#Test Stat   P-Value
#972.68656   0.00025

#Make the graph
ggplot(df, aes(x = n , color=list)) +
  scale_x_continuous(trans = "log10")+
  stat_ecdf(geom = "step")+
  annotate(geom="text", x=30, y=0.5, label=paste("p_value",t_all_host[2],sep = " : ") ,
           color="darkgreen")+
  annotate(geom="text", x=30, y=0.45, label=paste("p_value",t_all_nested[2],sep = " : ") ,
           color="darkblue")

#Test of this difference in isoform number is only due to the longer size of the host genes as it was observed before that longer genes tend to have more isoforms
#Regression analysis on the relationship between gene size (taking the longest isoform) and number of isoform if you are a host or not

#Determine gene size using the longest isoform
#Calculate the size
all$size = all[,3] - all[,2]
#Organise with the longest first
all = arrange(all, -all$size)
#Remove the duplication
all = all[!duplicated(all[,5]),]
#Keep only the size and the gene name
all= all[,c(5,8)]

#Add the size to the number of isoform table
Isoform_size = merge.data.frame(Isoforms,all, by.x = 1, by.y = 1)

#Add the information is the gene is host or not
Host = hn_gencode$SYMBOL_Host

for (i in 1:nrow(Isoform_size)){
  if (Isoform_size[i,1] %in% Host) {
    Isoform_size[i,"Host"] = "TRUE"
  } else {
    Isoform_size[i,"Host"] = "FALSE"
  }
  print(i)
}

Isoform_size$Host = as.logical(Isoform_size$Host)
#Make individual tables
Isoform_size_Host = Isoform_size[Isoform_size$Host == TRUE,]
Isoform_size_other 

#Plot and regression analysis
ggplot() +
  geom_point(data = Isoform_size_other, aes(x = size, y = n), color = "#fbe4afff") +
  geom_smooth(data = Isoform_size_other, aes(x = size, y = n), method = "lm", se = TRUE, color = "#fbe4afff", level = 0.95) +  geom_point(data = Isoform_size_Host, aes(x = size, y = n),colour = "#1b1248c7") +
  geom_smooth(data = Isoform_size_Host, aes(x = size, y = n), method = "lm", se = TRUE, color = "#1b1248c7", level = 0.95) +
  labs(title = "Regression analysis Mouse",
       x = "Gene Size",
       y = "Isoforms number") +
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  theme_minimal()

#Calculate the regression model 
#Use the model to predict the number of isoforms for the median size
median(Isoform_size$size)
#Mouse
#8604
lm_other = lm(log(n) ~ log(size),data= Isoform_size_other)
summary(lm_other)
log_isoform_number = log(8604)*0.216858 - 1.139432
exp(log_isoform_number)
#[1] 2.28256
lm_host = lm(log(n) ~ log(size), data = Isoform_size_Host)
summary(lm_host)
log_isoform_number = log(8604)*0.185376 -0.564797
exp(log_isoform_number)
#[1] 3.048673