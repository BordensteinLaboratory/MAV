#devtools::install_github("PheWAS/PheWAS",force=TRUE)
library(PheWAS)
library(sqldf)
library(tidyr)

setwd("./")

PheWASWrapper <- function(file) {
  
  #Read in the genetic risk score.
  
  genotypes <- read.table("/PolyGenenicScore_Final.txt",as.is=T, header=T, sep="\t",quote = "")
  colnames(genotypes)[1] <- "ID" 
  covar <- read.table("/PheWAS_Demographics.txt",as.is=T, header=T, sep="\t",quote = "")
  colnames(covar)[1] <- "ID" 
  PC <- read.table("/eigenvectors_10.txt",as.is=T, header=T, sep="\t",quote = "")
  #Get the related individuals to exclude and get the Individuals over 18.
  relateds <- read.table("/related_MinAgeGTE18_Exclusion_PiHat0.2_Set1.txt",as.is=T, header=T, sep="\t",quote = "")
  AgeRange <- read.table("/Subj_MinAgeGTE18_20191223.txt",as.is=T, header=T, sep="\t",quote = "")
  SubjectLists <- read.table("/ApprovedSubjectLists/BV239_Ferguson_Biomarkers_Ferguson_Biomarkers.txt",as.is=T, header=T, sep="\t",quote = "")
  colnames(SubjectLists)[1] <- "FID"
  
  #Select First 5 PC from PC as covariates in regression
  
  covar_PC <- sqldf('
      select covar.*, PC1,PC2,PC3,PC4,PC5
         from covar, PC
         where covar.ID=PC.FID
	     and covar.ID not in (select FID from relateds)
             and covar.ID in (select FID from SubjectLists)
	     and covar.ID in (select FID from AgeRange)
         ;')
  
  #Force the age to be numeric.

  covar_PC$AgeMedian <- as.numeric(covar_PC$AGE_MEDIAN_ICD)

  #Load the transposed pheWAS data.

  load(paste0("/EuropeanAncestry/",file))

  #Standardize the score to mean 0, std 1.
  
  SCORE_std <- genotypes$SCORE1_AVG
  genotypes <- cbind(genotypes,SCORE_std)
  
  #Run the PheWAS.
  
  results=phewas(phenotypes,genotypes[,c("ID","SCORE_std")],cores=16,min.records=10,covariates=covar_PC[,c("ID","GENDER_EPIC","AgeMedian","PC1","PC2","PC3","PC4","PC5")],additive.genotypes=F)

  PhewasAnnotation <- read.csv("/phecode_definitions1.2.csv")
  
  LCI <- exp(results$beta-(1.96*results$SE))
  UCI <- exp(results$beta+(1.96*results$SE))
  negLogPval=-1*log10(results$p)
  results1 <- cbind(results,LCI,UCI,negLogPval)
  colnames(results1)[colnames(results1)=="OR"] <- "OddsRatio"
  
  
  results_annotate <- sqldf('select snp,JD_CODE, JD_STRING,CATEGORY_STRING, n_cases, n_controls,beta,SE, OddsRatio, LCI,UCI, p,negLogPval
                             from results1, PhewasAnnotation
                             where results1.phenotype=PhewasAnnotation.jd_code and n_cases>10
                             order by p;'  )
  print(nrow(results_annotate))  
  
  #Plot the data. 
  
  jpeg(paste0(file,"_PHEWAS_results_ICD9_10_Plot_MinAgeGTE18.jpg"))

  results_plot <- results_annotate[which(results_annotate$n_cases>100),];
  PheWASData_FDR <- p.adjust(results_plot$p, method="fdr")
  PheWASData_FDRsig <- PheWASData_FDR<0.1
  
  gpl <- ggplot(results_plot, aes(OddsRatio,negLogPval)) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size=25),
          axis.text.x = element_text(colour="grey20",size=20), axis.text.y = element_text(colour="grey20",size=20),
          legend.position="bottom",
          legend.title=element_blank())  +
    ggtitle(file) +
    theme(plot.title = element_text(face="bold", size=25, hjust=0.5)) +
    geom_point(aes(shape=PheWASData_FDRsig,color=PheWASData_FDRsig),size=5) +
    scale_shape_manual(labels=c("FDR>0.1", "FDR<0.1"),values=c(16, 17)) +
    scale_color_manual(labels=c("FDR>0.1", "FDR<0.1"),values=c('black','green4'))+
    xlab("Odds-ratio") +
    ylab("-Log(p-value)") +
    scale_size(guide = "none") 

  print(gpl)
  dev.off()
  
  write.table(results_annotate, paste0(file,"_PHEWAS_results_ICD9_10_MinAgeGTE18.txt") , sep="\t", row.names=FALSE,quote=FALSE)
}

PheWASWrapper("PheWAS_MEGA_2orMore_ICD9_10CM_EA.txt_NON_GenderSpecific.Rda") 
