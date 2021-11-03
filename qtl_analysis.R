library(tidyverse)
library(Qtlizer)
library(data.table)
library(RColorBrewer)

#This SNP data is based on the SNPSNAP (https://data.broadinstitute.org/mpg/snpsnap/) data of matched input SNPS (780 out of 908 submitted SNPS and 200 sets of the total matched datasets). It can be found on FigShare as "SNPSNAP_MATCHED_SNPS.csv". 

mla_matched<- fread(SNPSNAP DATA) %>% 
  select(1:201)
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

#Compute z-score and one-sided probability. FDR correct across tissue pvals. 
all_z<- all%>% 
  mutate(z = (mav_freq-mean)/(sd)) %>% 
  mutate(pval = pnorm(z,lower.tail = F))
fdr <- as.data.frame(p.adjust(all_z$pval, method = "fdr"))
final <- cbind(all_z,fdr) %>% 
  mutate( sig = ifelse( fdr < 0.05,"yes", "no"))

#Plot frequency of MAV and simulated tissue frequency among significant tissues in group barplot
ggplot(final, aes(x = reorder(tissue,-mav_freq), y = mav_freq, fill= sig)) +
  geom_bar(stat="identity", width=.91) + 
  scale_fill_manual(values = c("ivory4", "indianred3")) +
  theme_classic()+
  xlab("Tissue") +
  ylab("# of Distinct eQTL ")+
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(face = "bold", color = "black",size = 10, angle = 45,hjust = 1),
  ) 

