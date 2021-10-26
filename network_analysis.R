library(tidyverse)
library (igraph)
library(tidygraph)
library(ggraph)
library(circlize)

#This script was used to generate the chord diagram to interpret shared MAVs between BioVU phewas phenotypes.

df <- read_csv("PheWAS File") %>% 
  mutate(phenotype_simple = case_when(
    str_detect(JD_STRING,"Type")~"T1D",
    str_detect(JD_STRING,"Pso")~"Psoriasis",
    str_detect(JD_STRING,"Rheu")~"RA",
    str_detect(JD_STRING,"embo|Deep|coag") ~"DVT/PE",
    str_detect(JD_STRING,"Gout") ~"Gout",
    str_detect(JD_STRING,"pul") ~"Pulmonary Heart Disease",
    TRUE ~ JD_STRING)) %>%
  select(SNP, JD_STRING) %>% 
  rename(source = SNP, destination = JD_STRING) %>% 
  group_by(destination) %>% 
  distinct(source, .keep_all = T)
netwk <- full_join(df, df, c('source' = 'source')) %>%
  group_by(destination.x, destination.y) %>% summarise(freq = n())  %>% 
  filter( destination.x != destination.y) 
netwk <- unique(as.data.frame(t(apply(netwk, 1, sort )))) %>% 
  select(,2:3)
chordDiagram(netwk,  annotationTrack = c("name", "grid"),transparency = .3, link.lwd = 1, link.lty = 2, link.border = "black",annotationTrackHeight = c(0.03, 0.05))
circos.clear()




