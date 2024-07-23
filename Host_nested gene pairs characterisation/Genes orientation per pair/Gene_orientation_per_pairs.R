# Analysis of the orientation of the pairs #

##################
#### HUMAN #######
##################

#Open the table generated after filtering of the pairs
hn_gencode <- read.table("Host_nested_genes_GENCODEV36lift37.txt",sep = '\t', header = T)

library(dplyr)
library(ggplot2)

#Count the number of pairs with the different possible orientation using the strand information
strand_count <- hn_gencode %>% group_by(Strand_Nested, Strand_Host) %>% count() %>% as.data.frame()
same_strand <- strand_count$n[1] + strand_count$n[4]
opposing_strand <- strand_count$n[2] + strand_count$n[3]

#Make a pie chart with the information on orientation
one <- data.frame(Strandedness = c("same_strand", "opposing_strand"),
                  n = c(same_strand, opposing_strand)) %>%
  ggplot(aes(x="", y=n, fill=Strandedness))+
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("red", "steelblue")) +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%", "\n","(", n,")")),
            position = position_stack(vjust = 0.5), size = 8)
one

#Add the strand information to the main table
hn_gencode$Strandness = hn_gencode$Strand_Nested == hn_gencode$Strand_Host
summary(hn_gencode$Strandness)
hn_gencode$Strandness = sub(TRUE, 'same', hn_gencode$Strandness)
hn_gencode$Strandness = sub(FALSE, 'opposite', hn_gencode$Strandness)

#Save the new table with orientation information
write.table(hn_gencode, file= "Host_nested_genes_GENCODEV36lift37_strand.txt", col.names = T, row.names = F, sep = '\t')

##################
#### MOUSE #######
##################

#Open the table generated after filtering of the pairs
hn_gencode <- read.table("Host_nested_genes_wgEncodeGencodeCompVM25.txt",sep = '\t', header = T)

library(dplyr)
library(ggplot2)

#Count the number of pairs with the different possible orientation using the strand information
strand_count <- hn_gencode %>% group_by(Strand_Nested, Strand_Host) %>% count() %>% as.data.frame()
same_strand <- strand_count$n[1] + strand_count$n[4]
opposing_strand <- strand_count$n[2] + strand_count$n[3]

#Make a pie chart with the information on orientation
one <- data.frame(Strandedness = c("same_strand", "opposing_strand"),
                  n = c(same_strand, opposing_strand)) %>%
  ggplot(aes(x="", y=n, fill=Strandedness))+
  geom_bar(width = 1, stat = "identity", colour = "black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("red", "steelblue")) +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%", "\n","(", n,")")),
            position = position_stack(vjust = 0.5), size = 8)
one

#Add the strand information to the main table
hn_gencode$Strandness = hn_gencode$Strand_Nested == hn_gencode$Strand_Host
summary(hn_gencode$Strandness)
hn_gencode$Strandness = sub(TRUE, 'same', hn_gencode$Strandness)
hn_gencode$Strandness = sub(FALSE, 'opposite', hn_gencode$Strandness)

#Save the new table with orientation information
write.table(hn_gencode, file= "Host_nested_genes_wgEncodeGencodeCompVM25_strand.txt", col.names = T, row.names = F, sep = '\t')
