library(tidyverse)
library(Qtlizer)
library(data.table)
library(RColorBrewer)
setwd("~/R/MAVs")
mla_matched<- fread("matched_snps.txt") 

snp_list <- mla_matched$Input_SNP
snp_list <-paste("hg19:",snp_list, sep = "")

df_mav <- get_qtls(snp_list) %>% 
  filter(str_detect(source,"GTE")) %>% 
  filter(sign_info == "FDR<5%" ) %>%
  rename( SNP = query_term) %>% 
  separate(tissue, into = cbind("tissue", "other"), sep= "-") %>% 
  select(SNP, tissue)  %>% 
  group_by(tissue) %>% 
  distinct(SNP)%>% 
  tally()

df_mav_deepdive <- get_qtls(snp_list) %>% 
  filter(str_detect(source,"GTE")) %>% 
  filter(sign_info == "FDR<5%" ) %>%
  rename( SNP = query_term) %>% 
  separate(tissue, into = cbind("tissue", "other"), sep= "-") 

significant_tissues_qtl <- df_mav_deepdive%>% 
  filter(tissue %in% c("Whole blood","Skin ","Adipose ","Pancreas","Esophagus ","Artery ","Lung","Nerve ","Testis","Spleen","Thyroid","Colon ","Pituitary","Heart ","Muscle skeletal")) 


####Simulation and testing of matched sets 
n <- 200
counts <- NULL
for (i in 1:n){
setchoice <- paste("Set_",i, sep= "")
snp_list <- mla_matched[[setchoice]]
snp_list <-paste("hg19:",snp_list, sep = "")
df <- get_qtls(snp_list) %>% 
  filter(str_detect(source,"GTE")) %>% 
  filter(sign_info == "FDR<5%" ) %>%
  rename( SNP = query_term)   %>% 
  separate(tissue, into = cbind("tissue", "other"), sep= "-") %>% 
  select(SNP, tissue)  %>% 
  group_by(tissue) %>% 
  distinct(SNP)%>% 
  tally()
names(df)[2] <- paste(setchoice, sep="")
counts[[i]] <- df
}
done <- reduce(counts,full_join) 
done[is.na(done)] <- 0
done_background <- done %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(Set_1:Set_200)), sd = sd(c_across(Set_1:Set_200)))%>% 
  select(tissue,mean,sd)
all <- full_join(done_background,df_mav) %>% 
  rename(mav_freq = n)
all[is.na(all)] <- 0

#compute z-score and one-sided probability. FDR correct across tissues. 
all_z<- all%>% 
  mutate(z = (mav_freq-mean)/(sd)) %>% 
  mutate(pval = pnorm(z,lower.tail = F))
fdr <- as.data.frame(p.adjust(all_z$pval, method = "fdr"))
final <- cbind(all_z,fdr) %>% 
  mutate( sig = ifelse( fdr < 0.051,"yes", "no"))

final <- read_tsv("matchedplotqtl.txt")%>% 
  mutate(Tissue = fct_reorder(Tissue, desc(-FDR))) 

  ggplot(final, aes(Tissue, Frequency, fill = SNP)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("indianred3","ivory4")) +
  theme_classic()+
  xlab("Tissue") +
  ylab("# of Distinct eQTL ")+
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(face = "bold", color = "black",size = 10, angle = 45,hjust = 1),
  ) 

