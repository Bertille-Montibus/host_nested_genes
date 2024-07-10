# Number of nested genes per host #

library("dplyr")

##################
#### HUMAN #######
##################

#Open the table with the list of host/nested gene pairs
hn_gencode <- read.delim("Host_nested_genes_GENCODEV36lift37_strand_ID.txt", row.names = 1)


#Count first the number of nested genes per host and then the number of host genes with the same number of nested genes
NestedPerHost_human = hn_gencode %>% count(SYMBOL_Host) %>% count(n)
NestedPerHost_human$n = as.character(NestedPerHost_human$n)


##################
#### MOUSE #######
##################

#Open the table with the list of host/nested gene pairs
hn_gencode <- read.delim("Host_nested_genes_wgEncodeGencodeCompVM25_strand_ID.txt", row.names = 1)

NestedPerHost = hn_gencode %>% count(SYMBOL_Host) %>% count(n)
NestedPerHost$n = as.character(NestedPerHost$n)

#Pool the data from mouse and human
df = merge(NestedPerHost, NestedPerHost_human, by.x = 1, by.y = 1, all = T)
#Sort the data by descending number of nested genes per host
df$n = as.numeric(df$n)
#Sort by the number of nested per host
df <- df %>% arrange(desc(n))
df$n = as.character(df$n)

#Caculate the 95th percentile for mouse and human
#Select only the mouse data
Mouse = df[,1:2]
Mouse = Mouse[!is.na(Mouse$nn.x),]
#Caculate the 95th percentile
stat = as.numeric(rep(Mouse$n, times = Mouse$nn.x))
quantile(stat, 0.95)
#Result =3

#Select only the human data
Human = df[,c(1,3)]
Human = Human[!is.na(Human$nn.y),]
stat = as.numeric(rep(Human$n, times = Human$nn.y))
quantile(stat, 0.95)
#Result =4

#Plot the data (in red mouse and in black human)
dotchart(df$nn.y, labels = df$n,
         cex = 0.7, xlab = "Number", ylab = "Number of nested genes/host", log = "x", bg = "green", pt.cex = 1.2, pch = 17)
points(df$nn.x, jitter(1:nrow(df)), col = "red", pch = 19, cex = 1)
legend("bottomright",
       legend = c("Human", "Mouse"),
       col = c("black","red"),
       pch = c(17,19),
       bty = "n",
       pt.cex = 1,
       cex = 1.2,
       text.col = "black",
       horiz = F ,
       inset = c(0.1, 0.1))
abline(h = 26,lty = "dotted",lwd = 2,col ="red")
abline(h = 25,lty = "dotted",lwd = 2)
text(2, 27, "95% percentile", col = 'red')
text(2, 24, "95% percentile")
