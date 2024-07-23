# Gene Ontology Analysis #

#Ontology analysis was perform using gProfiler (https://biit.cs.ut.ee/gprofiler/gost)
#Version used was e109_eg56_p17_1d3191d
#List of ENTREZGENE_ACC was used
#Bonferroni correction and a threshold of significance of 0.05 were applied
#Tables were next downloaded and used to make figure on R

library(dplyr)
library(ggpubr)

##################
#### HUMAN #######
##################

#Open the tables from gProfiler for the host genes
host_ontology <- read.csv("gProfiler_hsapiens_host_genes_human.csv", header = T)

#Calculate the observed/expected ratio using the intersection size and the term size
host_ontology <- host_ontology %>% dplyr::select(-intersections)
host_ontology <- host_ontology %>% mutate(observed = intersection_size/query_size) %>%
  mutate(expected = term_size/effective_domain_size) %>%
  mutate(o_e = observed/expected)

#selection of "Biological Processes" to make the figure including the 20 top terms for the observed/expected ratio, with the number of occurences and the adjusted pvalue
host_ontology %>% filter(source == "GO:BP") %>% slice_max(o_e, n=20) %>% ggplot(aes(x = reorder(term_name, o_e), y = o_e,col = negative_log10_of_adjusted_p_value, size =intersection_size )) +
  geom_point() + coord_flip()+
  labs(size="Count", colour="-log10(adj p-value)") + xlab("Biological Process") + ylab("Observed/Expected")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

#Open the tables from gProfiler for the nested genes
nested_ontology <- read.csv("gProfiler_hsapiens_nested_genes_human.csv", header = T)

#Calculate the observed/expected ratio using the intersection size and the term size
nested_ontology <- nested_ontology %>% dplyr::select(-intersections)
nested_ontology <- nested_ontology %>% mutate(observed = intersection_size/query_size) %>%
  mutate(expected = term_size/effective_domain_size) %>%
  mutate(o_e = observed/expected)

#selection of "Biological Processes" to make the figure including the 20 top terms for the observed/expected ratio, with the number of occurences and the adjusted pvalue
nested_ontology %>% filter(source == "GO:BP") %>% slice_max(o_e, n=20) %>% ggplot(aes(x = reorder(term_name, o_e), y = o_e,col = negative_log10_of_adjusted_p_value, size =intersection_size )) +
  geom_point() + bar_theme + coord_flip()+
  labs(size="Count", colour="-log10(adj p-value)") + xlab("Biological Process") + ylab("Observed/Expected")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

##################
#### MOUSE #######
##################

#Open the tables from gProfiler for the host genes
host_ontology <- read.csv("gProfiler_mmusculus_host_genes_mouse.csv", header = T)

#Calculate the observed/expected ratio using the intersection size and the term size
host_ontology <- host_ontology %>% dplyr::select(-intersections)
host_ontology <- host_ontology %>% mutate(observed = intersection_size/query_size) %>%
  mutate(expected = term_size/effective_domain_size) %>%
  mutate(o_e = observed/expected)

#selection of "Biological Processes" to make the figure including the 20 top terms for the observed/expected ratio, with the number of occurences and the adjusted pvalue
host_ontology %>% filter(source == "GO:BP") %>% slice_max(o_e, n=20) %>% ggplot(aes(x = reorder(term_name, o_e), y = o_e,col = negative_log10_of_adjusted_p_value, size =intersection_size )) +
  geom_point() + coord_flip()+
  labs(size="Count", colour="-log10(adj p-value)") + xlab("Biological Process") + ylab("Observed/Expected")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))

#Open the tables from gProfiler for the nested genes
nested_ontology <- read.csv("gProfiler_mmusculus_nested_genes_mouse.csv", header = T)

#Calculate the observed/expected ratio using the intersection size and the term size
nested_ontology <- nested_ontology %>% dplyr::select(-intersections)
nested_ontology <- nested_ontology %>% mutate(observed = intersection_size/query_size) %>%
  mutate(expected = term_size/effective_domain_size) %>%
  mutate(o_e = observed/expected)

#selection of "Biological Processes" to make the figure including the 20 top terms for the observed/expected ratio, with the number of occurences and the adjusted pvalue
nested_ontology %>% filter(source == "GO:BP") %>% slice_max(o_e, n=20) %>% ggplot(aes(x = reorder(term_name, o_e), y = o_e,col = negative_log10_of_adjusted_p_value, size =intersection_size )) +
  geom_point() + bar_theme + coord_flip()+
  labs(size="Count", colour="-log10(adj p-value)") + xlab("Biological Process") + ylab("Observed/Expected")+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 90, vjust = 0.5, hjust=1))
