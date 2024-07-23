# Analysis of the number of isoforms per genes to see if there is a difference for the host and nested genes #

library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
library(dplyr)
library(twosamples)
library(car)
library(funModeling)
library(emmeans)

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

#ANCOVA to evaluate if the relationship between size and isoform number is influenced by the host status

#Linear model testing
lm_full <- lm(log(n)~log(size) + 
                Host + 
                log(size) +
                Host:log(size),
              data = Isoform_size)
Anova(lm_full)
#Anova Table (Type II tests)

#Response: log(n)
#Sum Sq    Df   F value    Pr(>F)    
#log(size)      12732.6     1 18647.898 < 2.2e-16 ***
#  Host             395.0     1   578.556 < 2.2e-16 ***
#  log(size):Host    26.6     1    38.957 4.372e-10 ***
#  Residuals      29419.4 43087                        
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The interaction log(size):Host is significant(p=4.372e-10), so the evolution of isoform numbers according to gene size is
#different in host and non-host. 

#Then test if the isoform numbers of host and non-host are different when adjusted to take into account the effect of age.
#Estimation of the adjusted predicted means 
emmeans(lm_full, specs=pairwise ~ Host)

#Host  emmean      SE    df lower.CL upper.CL
#FALSE   1.39 0.00636 43087     1.37      1.4
#TRUE   1.68 0.01001 43087     1.66      1.7

#Results are given on the log (not the response) scale. 
#Confidence level used: 0.95 

#$contrasts
#contrast     estimate     SE    df t.ratio p.value
#FALSE - TRUE   -0.294 0.0119 43087 -24.768  <.0001

#The marginal mean for isoform number for host genes is exp(1.68) = 5.365556 and for the other genes is exp(1.39)= 4.01485 (when size effect is removed)
#The difference between the two is significantly different p<0.0001

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
Isoform_size_other = Isoform_size[Isoform_size$Host == FALSE,]

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

#ANCOVA to evaluate if the relationship between size and isoform number is influenced by the host status

#Linear model testing
lm_full <- lm(log(n)~log(size) + 
                Host + 
                log(size) +
                Host:log(size),
              data = Isoform_size)
Anova(lm_full)
#Anova Table (Type II tests)

#Response: log(n)
#Sum Sq    Df   F value    Pr(>F)    
#log(size)       7519.8     1 16183.242 < 2.2e-16 ***
#  Host             208.3     1   448.340 < 2.2e-16 ***
#  log(size):Host     9.3     1    20.101 7.368e-06 ***
#  Residuals      17616.0 37911                        
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The interaction log(size):Host is significant(p=7.368e-06), so the evolution of isoform numbers according to gene size is
#different in host and non-host. 

#Then test if the isoform numbers of host and non-host are different when adjusted to take into account the effect of gene size.
#Estimation of the adjusted predicted means 
emmeans(lm_full, specs=pairwise~Host)


#Host  emmean      SE    df lower.CL upper.CL
#FALSE   1.11 0.00513 37911     1.10     1.12
#TRUE   1.36 0.01023 37911     1.34     1.38

#Results are given on the log (not the response) scale. 
#Confidence level used: 0.95 

#$contrasts
#contrast     estimate     SE    df t.ratio p.value
#FALSE - TRUE   -0.247 0.0114 37911 -21.619  <.0001

#The marginal mean for isoform number for host genes is exp(1.36) = 3.896193 and for the other genes is exp(1.11)= 3.034358 (when size effect is removed)
#The difference between the two is significantly different p<0.0001