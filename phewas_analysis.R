library(tidyverse)
library(PheWAS)
library(readr)
library(VennDiagram)
library(ggVennDiagram)
library(CMplot)
library(ggplot2)
library(RColorBrewer)
library(hudson)
library(cowplot)
library(ggrepel)

dir <- "/Users/goblinking/Projects/MLA_local/BIOVU/"
e_population <-"european/"
a_population <-"african/"

file <- "/PheWAS_MEGA_2orMore_ICD9_10CM_EA.txt_NON_GenderSpecific.Rda_PHEWAS_results_ICD9_10_MinAgeGTE18.txt"
all_snps <-list.files(paste(dir,e_population,sep = ""))
data_list <- vector("list", "length" = length(all_snps))
for (snp_dir in all_snps){
  filename = snp_dir
  df<- read_tsv(paste(dir,e_population,snp_dir,file, sep = "")) %>% 
  filter(beta != "NA" & n_cases > 200)
  data_list[[snp_dir]] <- df
}
done <- do.call(rbind, data_list) 
ea_biovu <-  rownames_to_column(done, var = "snpID")%>% 
  separate(snpID,into=c("SNP","POS","A1","A2"), sep=":")%>% 
  mutate(adj_p_fdr= p.adjust(p, method="fdr")) 
ea_done <- ea_biovu %>%
  filter(adj_p_fdr < 0.05)
  
#AA_Biovu
file <- "/PheWAS_MEGA_2orMore_ICD9_10CM_AA.txt_NON_GenderSpecific.Rda_PHEWAS_results_ICD9_10_MinAgeGTE18.txt"
all_snps <-list.files(paste(dir,a_population,sep = ""))
data_list <- vector("list", "length" = length(all_snps))
for (snp_dir in all_snps){
  filename = snp_dir
  df<- read_tsv(paste(dir,a_population,snp_dir,file, sep = "")) %>% 
  filter(beta != "NA" & n_cases > 200)
  data_list[[snp_dir]] <- df
}
done <- do.call(rbind, data_list)
aa_biovu<-  rownames_to_column(done, var = "snpID") %>% 
separate(snpID,into=c("SNP","POS","A1","A2"), sep=":") %>% 
  mutate(adj_p_fdr= p.adjust(p, method="fdr"))
aa_done <- aa_biovu %>%
  filter(adj_p_fdr < 0.05)


# Manhattan Plots #

df <-  ea_done %>% 
  rename (phenotype = CATEGORY_STRING, unadj_p = p, p = adj_p_fdr,OR= OddsRatio) %>%
  type_convert() %>%
  # left_join(chr_map, by = "POS") %>%
  filter ( phenotype!= "NULL") %>% 
  mutate(OR.direction =case_when( OR >1 ~ "Positive", OR < 1 ~"Negative")) %>% 
  group_by(phenotype) %>% 
  type_convert()


PheWAS<- ggplot(df, aes(x = str_to_title(phenotype), y = -log(p,10)))+
  geom_point(aes(fill = phenotype, shape = OR.direction),size =5,position = "jitter")+
  scale_shape_manual(values = c(21, 21))+
  theme_classic() +
  scale_colour_manual()+
  theme(plot.title = element_text(vjust = 0.3, hjust = .4, size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0, size = 16),
        axis.text.y = element_text(angle = 90, vjust = 0.3, hjust = 0, size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(x = "", y = "-log(P)") +
  ylim(c(0, 13))+
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = .6, alpha = 0.8)

plot(PheWAS)
