library(tidyverse)
library(readr)
library(RColorBrewer)

dir <- DIRECTORY 
e_population <- EA_COHORT FILES 
a_population <- AA_COHORT FILES 
file <- DEFAULT PHEWAS RESULT FILENAME

#This code chunk is used to combine phewas results for each SNP tested in BioVU. FDR adjustment is applied across all SNPs and Traits and results are filtered to FDR < 0.05. 

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
  
#PheWAS summary statistics resulting from this is at FDR < 0.80 is available on FigShare 

# Manhattan Plotting is performed by first cleaning data and then adding OR direction to determine if OR is increased or decreased. If recreating the total plot as it appears in the figure, use FDR <0.80 to populate the non-significant results. 
df <-  ea_done %>% 
  rename(phenotype = CATEGORY_STRING, unadj_p = p,p = adj_p_fdr,OR=OddsRatio) %>%
  type_convert() %>%
  filter (phenotype!= "NULL") %>% 
  mutate(OR.direction =case_when( OR >1 ~ "Positive", OR < 1 ~"Negative")) %>% 
  group_by(phenotype) 
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
