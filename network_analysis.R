library(tidyverse)
library (igraph)
library(tidygraph)
library(ggraph)
library(circlize)

letters <- read_csv("phe10.csv") %>% 
  # mutate(phenotype_simple = case_when(
  #   str_detect(JD_STRING,"Type")~"T1D",
  #   str_detect(JD_STRING,"Pso")~"Psoriasis",
  #   str_detect(JD_STRING,"Rheu")~"RA",
  #   str_detect(JD_STRING,"embo|Deep|coag") ~"DVT/PE",
  #   str_detect(JD_STRING,"Gout") ~"Gout",
  #   str_detect(JD_STRING,"pul") ~"Pulmonary Heart Disease",
  #   TRUE ~ JD_STRING)) %>% 
  select(SNP, JD_STRING) %>% 
  rename(source = SNP, destination = JD_STRING) %>% 
  group_by(destination) %>% 
  distinct(source, .keep_all = T)
 
mat <- full_join(letters, letters, c('source' = 'source')) %>%
  group_by(destination.x, destination.y) %>% summarise(freq = n())  %>% 
  filter( destination.x != destination.y) 
  

mat <- unique(as.data.frame(t(apply(mat, 1, sort )))) %>% 
  select(,2:3)

mat <- read_tsv("connections.txt")


chordDiagram(mat, annotationTrack = "grid",transparency = .1, link.lwd = 2, link.lty = 1, link.border = "black")

chordDiagram(mat, annotationTrack = "grid")
# 
chordDiagram(mat,  annotationTrack = c("name", "grid"),transparency = .3, link.lwd = 1, link.lty = 2, link.border = "black",annotationTrackHeight = c(0.03, 0.05))


# grid.col = c(`Abnormal movement` = "red", `Celiac disease` = "green", DVT = "blue",
#              `Multiple sclerosis`= "grey", `Other demyelinating diseases of central nervous system` = "grey", `Other inflammatory spondylopathies` = "grey", `Primary hypercoagulable state` = "grey", Psoriasis = "grey", `Pulmonary heart disease` = "grey", `Rheumatoid arthritis` = "blue", T1D = "grey")

chordDiagram(mat, grid.col = grid.col, transparency = .3, link.lwd = 2, link.lty = 1, link.border = "black")

circos.clear()


# 
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
# }, bg.border = NA) # here set bg.border to NA is important
# 
# circos.clear()


# chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
# circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim")
#   xplot = get.cell.meta.data("xplot")
#   ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
# 
#     circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
#                 niceFacing = TRUE, adj = c(0, 0.5))
# 
# }, bg.border = NA)


