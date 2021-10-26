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
  filter(adj_p_fdr < 0.80)
  
#AA_Biovu
file <- "/PheWAS_MEGA_2orMore_ICD9_10CM_AA.txt_NON_GenderSpecific.Rda_PHEWAS_results_ICD9_10_MinAgeGTE18.txt"
all_snps <-list.files(paste(dir,a_population,sep = ""))
data_list <- vector("list", "length" = length(all_snps))
for (snp_dir in all_snps){
  filename = snp_dir
  df<- read_tsv(paste(dir,a_population,snp_dir,file, sep = "")) %>% 
  filter(beta != "NA" & n_cases > 50)
  data_list[[snp_dir]] <- df
}
done <- do.call(rbind, data_list)
aa_biovu<-  rownames_to_column(done, var = "snpID") %>% 
separate(snpID,into=c("SNP","POS","A1","A2"), sep=":") %>% 
  mutate(adj_p_fdr= p.adjust(p, method="fdr"))
aa_done <- aa_biovu %>%
  filter(adj_p_fdr < 0.05)

#Read In & Preprocess####

chr_map <- read_tsv("/Users/goblinking/Dropbox/POS2CHR.txt")

#Mapping of trait to snps

ea_200<- ea_biovu %>% filter(n_cases > 200 & adj_p_fdr < 0.05) %>% 
  distinct(SNP)


#Analysis####

ea <- ea_biovu %>% filter(n_cases > 200 & adj_p_fdr < 0.05) %>% 
  select(-snp) %>% 
  type_convert() %>% 
  separate (A2, into = c("A2","Freq"), sep= "\\.", remove=T, convert = T)


# Manhattan Plots #

df <-  aa_done %>% 
  rename (phenotype = CATEGORY_STRING, unadj_p = p, p = adj_p_fdr,OR= OddsRatio) %>%
  type_convert() %>%
  # left_join(chr_map, by = "POS") %>%
  filter ( phenotype!= "NULL") %>% 
  mutate(OR.direction =case_when( OR >1 ~ "Positive", OR < 1 ~"Negative")) %>% 
  group_by(phenotype) %>% 
  type_convert()


PheWAS_2 <- ggplot(df, aes(x = str_to_title(phenotype), y = -log(p,10)))+
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

plot(PheWAS_2)

str_to_title()
#Chromosome Position Plot 
ggplot(ea_b, aes(x=CHR, y=-log(p))) + 
  # geom_point(aes(col=phenotype)) +
  geom_jitter(aes(color = phenotype, size =5))+
  theme_classic() + 
  theme(panel.grid.minor=element_blank(),legend.position = "none") +
  labs(coor="Phenotype Category", x="Chromosome", y="log(p-value)") +
  geom_hline(yintercept= 1, color="grey", size=.5, alpha=0.8) +
  geom_hline(yintercept= 2, color="red", size=.5, alpha=0.8) +
  geom_text_repel(data=. %>% mutate(label = ifelse(p < 0.05, as.character(SNP), "")), aes(label=label), size=2, box.padding = unit(0.7, "lines"))



#####Mapping of trait to snps

mla_mb <- read_tsv("mla_mb.txt")
ea_final_all_traits <-  ea_final %>%
  type_convert() %>%
  left_join(chr_map, by = "POS") %>% 
  left_join(mla_mb, by = "SNP") %>% 
  distinct(`SNP`, .keep_all = T) %>% 
  select(SNP,mla_trait, Study, OddsRatio) %>% 
  write_csv("trait_snp_beta.csv")

all_traits <-  ea_filtered_200_fdr05 %>%
  type_convert() %>%
  left_join(chr_map, by = "POS") %>% 
  # left_join(mla_mb, by = "SNP") %>% 
  distinct(`SNP`, .keep_all = T) %>% 
  select(SNP,CHR) 