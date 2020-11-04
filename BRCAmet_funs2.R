#Project: Deep Learning Classifier to estimate the probability of future metastasis from Breast Primary Tumor WES data
#@feBueno fernando.bueno.gutierrez@usal.es August-2020

#This file contains functions for BRCAmet_getInput.R and BRCAmet_funs1.R
#Goal of functions in this file is to select a subset of highly informative mutation_positions 

saveGroupMpInGenes_fun<-function(mut_df){
  
  # head(mut_df)
  # Hugo_Symbol Sample_Barcode mutation_position dMB source
  # 1      ATPAF1   TCGA-BH-A18F     chr1_46665282  No    GDC
  # 2        RORC   TCGA-BH-A18F    chr1_151814660  No    GDC
  
  #head(head(mutRed_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-BH-A18F            ATPAF1  No    GDC
  # 2   TCGA-BH-A18F              RORC  No    GDC
  
  mut_df$mutation_position<-mut_df$Hugo_Symbol
  mut_df$Hugo_Symbol<-NULL
  mutRed_df=unique(mut_df)
  
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByGroupMpInGenes_df.txt",col.names = T, row.names = F, quote = F)
  
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByGroupMpInGenes_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
}#WARNING! will affect also the test set

saveMpWindows_fun<-function(mut_df,nWindowsPerChr_i=1000){
  
  # head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-BH-A18F     chr1_46665282  No    GDC
  # 2   TCGA-BH-A18F    chr1_151814660  No    GDC
  
  #nWindowsPerChr_i:1000, number of windows or ranges, each with equal base-pairs-length, in which each chromosome will be split
    #example_v=c(8,44,3,2,12,3,4,9)
    #split(mv, cut(example_v, 3))
    # $`(1.96,16]`
    # [1]  8  3  2 12  3  4  9
    # 
    # $`(16,30]`
    # numeric(0)
    # 
    # $`(30,44]`
    # [1] 44
  
  # head(mutRed_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-BH-A18F           C1_R188  No    GDC
  # 2   TCGA-BH-A18F           C1_R610  No    GDC
  
  mut_df$startPosition=sub(".*_","",mut_df$mutation_position)
  mut_df$chr=sub("_.*","",mut_df$mutation_position)
  mut_df$chr=gsub("[^0-9.-]", "", mut_df$chr)
  mut_df$chr[mut_df$chr==""]<-"X"
  mut_df$range="toBfilled"
  
  fillRangenForChr_fun<-function(chrNumber_c){
    #updates mut_df$range
    
    #chrNumber_c: one of "1", "2"..  "22", "X"
    
    mpChr_idx=which(mut_df$chr==chrNumber_c)#index of all the mutation positions in chrNumber_c
    mpChr_n=as.numeric(mut_df$startPosition[mut_df$chr==chrNumber_c])
    mpInEachRange_l=split(mpChr_n, cut(mpChr_n, nWindowsPerChr_i))
    getRanges_fun<-function(mpChr_idx){
      
      #mpChr_idx: 1, index of the mutation position in the chromosome, in mut_df order
      
      #range_i: 1 if the mpChr_idx mutation belongs to the leftiest part of the chromosome
      
      listName_c=names(mpInEachRange_l)[sapply(1:length(mpInEachRange_l),function(x){mpChr_n[mpChr_idx] %in% mpInEachRange_l[[x]]})]
      range_i=which(names(mpInEachRange_l) == listName_c)
      
      return(range_i) 
    }
    
    allRangesForChrMp_n=sapply(1:length(mpChr_n),getRanges_fun)#allRangesForChrMp_n: all the ranges (integers from 1 to nWindowsPerChr_i) for each of the mutation_positions n the chromosome
    rangeName_c=paste0("C",chrNumber_c,"_R",allRangesForChrMp_n)
    mut_df$range[mpChr_idx]<<-rangeName_c
  }
  
  sapply(unique(mut_df$chr),fillRangenForChr_fun)
  
  mutRed_df=mut_df[,c("Sample_Barcode","range","dMB","source")]
  colnames(mutRed_df)[2]="mutation_position"
  
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByMpWindows_df.txt",col.names = T, row.names = F, quote = F)
  
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByMpWindows_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
  
}#WARNING! will affect also the test set

saveMutRedConsensus_fun<-function(nTimesMPmustBeReported_i=1, excludeFiles_c=""){
  
  #nTimesMPmustBeReported_i: int (1), number of times that a mutation must be reported in the sum of combination of mutation reduction approaches
  #excludeFiles_c: char ("mutRedByBailey1_df.txt","mutRedByCGI_df.txt"). File names within "./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df" that will be ignored
  
  #head(mutRedCon_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  mutRedConsFileNames_c <- list.files("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df")
  mutRedConsFileNames_c=mutRedConsFileNames_c[!mutRedConsFileNames_c %in% excludeFiles_c]
  filePath_c=paste0("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df","/",mutRedConsFileNames_c[1])
  mutRedCon_df<-read.table(filePath_c,header=T)
  mutRedConsFileNames_c=mutRedConsFileNames_c[-1]
  for(fileName_c in mutRedConsFileNames_c){
    filePath_c=paste0("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df","/",fileName_c)
    newmutRedCon_df=read.table(filePath_c,header=T)
    mutRedCon_df=rbind(mutRedCon_df,newmutRedCon_df)
  }
  mpConsensusSelected_c=mutRedCon_df$mutation_position[table(unique(mutRedCon_df$mutation_position))>0]
  mutRedCon_df=mutRedCon_df[mutRedCon_df$mutation_position %in% mpConsensusSelected_c,]
  
  #save file
  write.table(mutRedCon_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedCon_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedCon_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRedCon_df)
} 

saveMutRedHPAprognosisGenes_fun<-function(mut_df){#HPA: human Protein Atlas. Keeps mutations in genes that are significantly associated with prognosis
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  library(readr)
  
  #download from https://www.proteinatlas.org/humanproteome/pathology and selecting the tissue. 5 tissues were selected: breast, liver and lung (since they are breast common metastasis and were avilable in this dataset), ovarian cervical (since they are related to breast and were avilable in this dataset) 
  #Endometrial cancer was not included. Although it is also related to breast, it points to too many genes
  #For instance, for breast download: https://www.proteinatlas.org/search/prognostic%3Abreast+cancer+AND+sort_by%3Aprognostic+breast+cancer+AND+show_columns%3Aprognostic
  HPAprognosisBreastGenes_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/humanProteinAtlas/prognostic_breast.tsv")#HPAprognosisBreastGenes_multiClass: 4-classes-object, HPAprognosis genes data
  HPAprognosisBreastGenes_df=as.data.frame(HPAprognosisBreastGenes_multiClass[,c(1,60)])
  colnames(HPAprognosisBreastGenes_df)<-c("Gene","BCprognosis")
    dim(HPAprognosisBreastGenes_df)#577   2
  HPAprognosisLiverGenes_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/humanProteinAtlas/prognostic_liver.tsv")#HPAprognosisBreastGenes_multiClass: 4-classes-object, HPAprognosis genes data
  HPAprognosisLiverGenes_df=as.data.frame(HPAprognosisLiverGenes_multiClass[,c(1,60)])
  colnames(HPAprognosisLiverGenes_df)<-c("Gene","BCprognosis")
  HPAprognosisLungGenes_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/humanProteinAtlas/prognostic_lung.tsv")#HPAprognosisBreastGenes_multiClass: 4-classes-object, HPAprognosis genes data
  HPAprognosisLungGenes_df=as.data.frame(HPAprognosisLungGenes_multiClass[,c(1,60)])
  colnames(HPAprognosisLungGenes_df)<-c("Gene","BCprognosis")
  HPAprognosisOvarianGenes_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/humanProteinAtlas/prognostic_ovarian.tsv")#HPAprognosisBreastGenes_multiClass: 4-classes-object, HPAprognosis genes data
  HPAprognosisOvarianGenes_df=as.data.frame(HPAprognosisOvarianGenes_multiClass[,c(1,60)])
  colnames(HPAprognosisOvarianGenes_df)<-c("Gene","BCprognosis")
  HPAprognosisCervicalGenes_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/humanProteinAtlas/prognostic_cervical.tsv")#HPAprognosisBreastGenes_multiClass: 4-classes-object, HPAprognosis genes data
  HPAprognosisCervicalGenes_df=as.data.frame(HPAprognosisCervicalGenes_multiClass[,c(1,60)])
  colnames(HPAprognosisCervicalGenes_df)<-c("Gene","BCprognosis")
  
  full5CancersPrognosisGenes_df=rbind(HPAprognosisBreastGenes_df,HPAprognosisLiverGenes_df)
  full5CancersPrognosisGenes_df=rbind(full5CancersPrognosisGenes_df,HPAprognosisLungGenes_df)
  #full5CancersPrognosisGenes_df=rbind(full5CancersPrognosisGenes_df,HPAprognosisOvarianGenes_df)
  #full5CancersPrognosisGenes_df=rbind(full5CancersPrognosisGenes_df,HPAprognosisCervicalGenes_df)
  dim(full5CancersPrognosisGenes_df)
  length(unique(full5CancersPrognosisGenes_df$Gene))#3838 genes
  
  HPAprognosis_mp_c=mut_df$mutation_position[mut_df$Hugo_Symbol %in% full5CancersPrognosisGenes_df$Gene]
  mutRed_df=mut_df[mut_df$mutation_position %in% HPAprognosis_mp_c,]
  mutRed_df$Hugo_Symbol<-NULL
  
  #save file
  write.table(mutRed_df,"./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByHPAprognosisGenes_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByHPAprognosisGenes_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)

}

saveMutRedBailey3_fun<-function(mut_df){
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  sel3mc3Bailey_12087mut_GRCh37_mafdf=readRDS("./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/sel3mc3Bailey_12087mut_GRCh37_mafdf.rds")
  
  #convert to HG38
  mutation_positions_hg38_c=hg19toHg38_fun(as.character(sel3mc3Bailey_12087mut_GRCh37_mafdf$mutation_position))
  saveRDS(mutation_positions_hg38_c,"./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/mutation_positions_hg38_sel3_c.rds")
  return("")
  
    #length(unique(mutationsFROMlifter_c))#12085
  
  #WITH CONVERSION FROM hg37 TO hg38
  mutRed_df=mut_df[mut_df$mutation_position %in% mutationsFROMlifter_c,]
    dim(mutRed_df)#696    4
    length(unique(mutRed_df$mutation_position))#431
  
  #WITHOUT CONVERSION  
  # mutRed_df=mut_df[mut_df$mutation_position %in% sel3mc3Bailey_12087mut_GRCh37_mafdf$mutation_position,]
  # dim(mutRed_df)#2794   4
  # length(unique(mutRed_df$mutation_position))#292
  
  #save file
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey3_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByBailey3_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
}#Bailey: https://www.cell.com/cell/comments/S0092-8674(18)30237-X  (HG37)

saveMutRedBailey2_fun<-function(mut_df){
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  sel2mc3Bailey_64985mut_GRCh37_mafdf=readRDS("./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/sel2mc3Bailey_64985mut_GRCh37_mafdf.rds")
  
  #convert to HG38
  mutation_positions_hg38_c=hg19toHg38_fun(as.character(sel2mc3Bailey_64985mut_GRCh37_mafdf$mutation_position))
  saveRDS(mutation_positions_hg38_c,"./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/mutation_positions_hg38_sel2_c.rds")
  return("")
  
  
    #length(unique(mutationsFROMlifter_c))#64985
  
  #WITH CONVERSION FROM hg37 TO hg38
  mutRed_df=mut_df[mut_df$mutation_position %in% mutationsFROMlifter_c,]
    dim(mutRed_df)#3254    4
    length(unique(mutRed_df$mutation_position))#3292
  
  #WITHOUT CONVERSION  
  # mutRed_df=mut_df[mut_df$mutation_position %in% sel2mc3Bailey_64985mut_GRCh37_mafdf$mutation_position,]
  #   dim(mutRed_df)#9374    4
  #   length(unique(mutRed_df$mutation_position))#2614
  
  #save file
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey2_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByBailey2_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
}#(HG37)

saveMutRedBailey1_fun<-function(mut_df){
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  sel1mc3Bailey_46561mut_GRCh37_mafdf=readRDS("./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/sel1mc3Bailey_46561mut_GRCh37_mafdf.rds")
  
  #convert to HG38
  mutation_positions_hg38_c=hg19toHg38_fun(as.character(sel1mc3Bailey_46561mut_GRCh37_mafdf$mutation_position))
  saveRDS(mutation_positions_hg38_c,"./BRCAmet/data/mutationsSelectionData/Clare_BaileyPaperMutSel/mutation_positions_hg38_sel1_c.rds")
  return("")
    #length(unique(mutationsFROMlifter_c))#46559
  
  #WITH CONVERSION FROM hg37 TO hg38
  mutRed_df=mut_df[mut_df$mutation_position %in% mutationsFROMlifter_c,]
    dim(mutRed_df)#3982    4
    length(unique(mutRed_df$mutation_position))#3292
  
  #WITHOUT CONVERSION  
  # mutRed_df=mut_df[mut_df$mutation_position %in% sel1mc3Bailey_46561mut_GRCh37_mafdf$mutation_position,]
  #   dim(mutRed_df)#6656    4
  #   length(unique(mutRed_df$mutation_position))#885
    
  #save file
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey1_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByBailey1_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
}#(HG37)

ClareSave3mafdf_rdsBailey_fun<-function(){
  
  #This function is to extract a set of cancer-relevant-mutations based on:
  #Comprehensive Characterization of Cancer Driver Genes and Mutations
  #https://www.cell.com/cell/comments/S0092-8674(18)30237-X
  
  #WARNING!
  #The script was run in Clare as it requires large RAM to open (and work with) the large full-TCGA-maf file (both: "mc3.v0.2.8.PUBLIC.maf" and "mc3_v0_2_8_PUBLIC_maf_data.rds")
  #https://gdc.cancer.gov/about-data/publications/mc3-2017; MC3 Public MAF - mc3.v0.2.8.PUBLIC.maf.gz
  
  #Based on paper, I chose 3 options for mutation selection.
  #All 3 options rely on the follwoing open-access-maf-file used in the aforomentioned paper: https://gdc.cancer.gov/about-data/publications/mc3-2017#MC3 Public MAF - mc3.v0.2.8.PUBLIC.maf.gz 
  #Based on paper, clinical data corresponding to full-TCGA-maf comes from (TCGA): https://gdc.cancer.gov/about-data/publications/pancanatlas  
  #The 3 options are:
  #1) Select from full-TCGA-maf (~1.8M mutation_positions) if count for that mutation_position (in whole TCGA samples) is 3 or more
  #2) Select from full-TCGA-maf (~1.8M mutation_positions) if mutation_position is in the 738 defined genes in S1 of Bailey paper (Supp table excel)
  #3) Select from full-TCGA-maf (~1.8M mutation_positions) if mutation_position is in the 578 amino-acid-code-change in S4 of Bailey paper (Supp table excel)
  
  #Other databases mentioned in paper:
  #Check this other paper they used to get clinical data: An Integrated TCGA Pan-Cancer Clinical Data Resource to drive high quality survival outcome analytics.
  #Database of Evidence for Precision Oncology (DEPO; http://depo-dinglab.ddns.net)
  #An Integrated TCGA Pan-Cancer Clinical Data Resource to drive high quality survival outcome analytics (http://oncokb.org)
  
  setwd("/home/fernando/Escritorio/")
  
  #start in LARGE file directly downloaded from https://gdc.cancer.gov/about-data/publications/mc3-2017#MC3 Public MAF - mc3.v0.2.8.PUBLIC.maf.gz 
  #library(maftools)
  #alldfs=read.maf("mc3.v0.2.8.PUBLIC.maf")#unziped after direct downloaded from https://gdc.cancer.gov/about-data/publications/mc3-2017#MC3 Public MAF - mc3.v0.2.8.PUBLIC.maf.gz 
  #saveRDS(alldfs@data,"mc3_v0_2_8_PUBLIC_maf_data.rds")
  
  #start in LARGE .rds file created in lines above
  fullTCGA_GRCh37_mafdf<-readRDS("mc3_v0_2_8_PUBLIC_maf_data.rds")
  fullTCGA_GRCh37_mafdf$mutation_position=paste0("chr",fullTCGA_GRCh37_mafdf$Chromosome,"_",fullTCGA_GRCh37_mafdf$Start_Position)
  dim(fullTCGA_GRCh37_mafdf)#2162922     115
  length(unique(fullTCGA_GRCh37_mafdf$mutation_position))#1861315
  
  #1st mutation selection approach based on Comprehensive-Characterization-Bailey-paper: mutations that appear >2 in the TCGAmaf file
  mutation_position_ta=table(fullTCGA_GRCh37_mafdf$mutation_position)
  mutation_position_aboveThreshold_ta=mutation_position_ta[mutation_position_ta>2*median(mutation_position_ta)]
  sel1_TCGA_GRCh37_mafdf=fullTCGA_GRCh37_mafdf[fullTCGA_GRCh37_mafdf$mutation_position %in% names(mutation_position_aboveThreshold_ta),]
  dim(sel1_TCGA_GRCh37_mafdf)#177943    115
  length(unique(sel1_TCGA_GRCh37_mafdf$mutation_position))#46561
  saveRDS(sel1_TCGA_GRCh37_mafdf,"sel1mc3Bailey_46561mut_GRCh37_mafdf.rds")
  #subset1 of the mc3 maf file (https://gdc.cancer.gov/about-data/publications/mc3-2017) mentioned in Bailey paper (https://www.cell.com/cell/comments/S0092-8674(18)30237-X); that has 46561 unique mutation_positions in GRCh37 NCBI_Build
  
  #2nd mutation selection approach based on Comprehensive-Characterization-Bailey-paper: 299 genes
  library(openxlsx)
  SuppTable2Bailey_df=read.xlsx("/home/fernando/Escritorio/mmc1.xlsx",sheet=2,colNames =T,startRow = 4)#Supplementary Table 2 from Bailey papar (duircet access from paper as mmc1.xlsx)
  length(unique(SuppTable2Bailey_df$Gene))#299 genes
  sel2_TCGA_GRCh37_mafdf=fullTCGA_GRCh37_mafdf[fullTCGA_GRCh37_mafdf$Hugo_Symbol %in% SuppTable2Bailey_df$Gene,]
  dim(sel2_TCGA_GRCh37_mafdf)#88902    115
  length(unique(sel2_TCGA_GRCh37_mafdf$mutation_position))#64985
  saveRDS(sel2_TCGA_GRCh37_mafdf,"sel2mc3Bailey_64985mut_GRCh37_mafdf.rds")
  #subset2 of the mc3 maf file (https://gdc.cancer.gov/about-data/publications/mc3-2017) mentioned in Bailey paper (https://www.cell.com/cell/comments/S0092-8674(18)30237-X); that has 64985 unique mutation_positions in GRCh37 NCBI_Build
  
  #3rd mutation selection approach based on Comprehensive-Characterization-Bailey-paper: 565 amino-acid changes in Supp.Table.5
  SuppTable4Bailey_df=read.xlsx("/home/fernando/Escritorio/mmc1.xlsx",sheet=5,colNames =T,startRow = 4)
  length(unique(SuppTable4Bailey_df$Mutation))#565 amino-acid-change mutation codes
  sel3_TCGA_GRCh37_mafdf=fullTCGA_GRCh37_mafdf[fullTCGA_GRCh37_mafdf$HGVSp_Short %in% SuppTable4Bailey_df$Mutation,]
  dim(sel3_TCGA_GRCh37_mafdf)#18726   115
  length(unique(sel3_TCGA_GRCh37_mafdf$mutation_position))#12087 mutation_positions
  saveRDS(sel3_TCGA_GRCh37_mafdf,"sel3mc3Bailey_12087mut_GRCh37_mafdf.rds")
  #subset3 of the mc3 maf file (https://gdc.cancer.gov/about-data/publications/mc3-2017) mentioned in Bailey paper (https://www.cell.com/cell/comments/S0092-8674(18)30237-X); that has 12087 unique mutation_positions in GRCh37 NCBI_Build
  
  print("The following 3 files have been saved in current directory (/home/fernando/Escritorio/):")
  print("sel1mc3Bailey_46561mut_GRCh37_mafdf.rds")
  print("sel2mc3Bailey_64985mut_GRCh37_mafdf.rds")
  print("sel3mc3Bailey_12087mut_GRCh37_mafdf.rds")
  print("Details about the source and use of these files can be found in this function's code")
  
}

saveMutRedCosmic_fun<-function(mut_df,timesAboveMedianCount_i=2){
  
  #Saves smaller mut_df in same format but with less mutation_positions: mutRed_df. mutation_positions are selected if they are in (requires login) https://cancer.sanger.ac.uk/cosmic/download [All Mutations in Census Genes/Download Whole File]
  #and they fulfill a specified-in-function-argument-criteria regarding the number of times that this mutation was found in the cancer-cesus database
  #For login, it is enough with having email from a Scintific Insititution (i.e. University)
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  #timesAboveMedianCount_i#number of times of counts that a mutation must have above the median counts (for all mutations) to be kept. 
    #Shall you provide a float (i.e. 2.3), timesAboveMedianCount_i will be effectivelly 3 (the ceiling integer) 
  
  #head(mutRed_df)#mutRed_df: a reduced version of mut_df with fewer mutation_position  
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  
  #read large cancer-census-mutations dataset
  library(readr)
  cosmicMut_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/cosmic/CosmicMutantExportCensus.tsv")#cosmicMut_multiClass: 4-classes-object, cosmic-database mutations
  head(sort(table(cosmicMut_multiClass$`Mutation genome position`),decreasing=T))
  # 9:5073770-5073770 7:140753336-140753336  12:25245350-25245350  12:25245351-25245351  12:25245347-25245347 2:208248388-208248388 
  # 42467                 29145                 28863                  8888                  5969                  5274 
  length(unique(cosmicMut_multiClass$`Mutation genome position`))#922732
  
  #get character with the most frecuent mutations    
  medianMutCounts_i=median(table(cosmicMut_multiClass$`Mutation genome position`))#medianMutCounts_i: integer, mediant of mutation counts #1
  print("median:")
  print(medianMutCounts_i)
  print(class(medianMutCounts_i))
  thresholdMutCounts_i=medianMutCounts_i*timesAboveMedianCount_i#thresholdMutCounts_i: integer, threshold of mutation counts. Mutation-positions above threshold will be selected #2
  selMut_c=as.character(cosmicMut_multiClass$`Mutation genome position`[table(cosmicMut_multiClass$`Mutation genome position`)>thresholdMutCounts_i])#selMut_c: character, selected mutations
    length(unique(selMut_c))#67652
    head(selMut_c)#[1] "4:152503513-152503513" "1:156769051-156769051"
  
  #change to mutation_position format (as in mut_df)   
  selMut_c=gsub(":","_",selMut_c)
  selMut_c<-gsub("-.*","",selMut_c)
  selMut_c<-paste0("chr",selMut_c)
  head(selMut_c)#[1] "chr4_152503513" "chr1_156769051"
  mutRed_df=mut_df[mut_df$mutation_position %in% selMut_c,]
  
  #save file
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByCensus_df.txt",col.names = T, row.names = F, quote = F)
  
  #print main data properties
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByCensus_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
  
}#(HG38)

saveMutRedCGI_fun<-function(mut_df){
  
  #Saves smaller mut_df in same format but with less mutation_positions: mutRed_df. mutation_positions are selected if they are in https://www.cancergenomeinterpreter.org/mutations
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  #head(mutRed_df)#mutRed_df: a reduced version of mut_df with fewer mutation_position  
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC

  #load CGI mutations dataset
  library(readr)
  CGImut_multiClass=read_tsv("./BRCAmet/data/mutationsSelectionData/CancerGenomeInterpreter/catalog_of_validated_oncogenic_mutations_latest/catalog_of_validated_oncogenic_mutations.tsv")
  head(CGImut_multiClass$gdna)
  chr_c=sub(":.*","",CGImut_multiClass$gdna)
  after_c=sub(".*:","",CGImut_multiClass$gdna)
  after_c=sub("_.*","",after_c)
  pos_c=gsub("[^0-9.-]", "", after_c)
  pos_c=as.numeric(sub("\\.","",pos_c))
  mutation_position_suppossedHg37_c=paste0(chr_c,"_",pos_c)#https://www.cancergenomeinterpreter.org/home does not indicate wether it is Hg37 or hg38. We supposed it is hg37
    #Explanation below in this function
  
  mut_df=read.table("/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mut_hg38_df.txt",header=T)
  mut_df$Hugo_Symbol=NULL
  
  #load mutations converted to Hg38 or make conversion
  mutation_position_CGI_hg38_c=try(readRDS("/home/febueno/Documents/BRCAmetProject/BRCAmet/data/mutationsSelectionData/CancerGenomeInterpreter/mutation_position_5601CGI_hg38_c.rds"))
  if(class(mutation_position_CGI_hg38_c)=="try-error"){
    mutation_position_CGI_hg38_c=hg19toHg38_fun(mutation_position_suppossedHg37_c)#~5 minutes
    saveRDS(mutation_position_CGI_hg38_c,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/mutationsSelectionData/CancerGenomeInterpreter/mutation_position_5601CGI_hg38_c.rds")
  }
  
  mutRed_df=mut_df[as.character(mut_df$mutation_position) %in% mutation_position_CGI_hg38_c,]
    length(unique(mutRed_df$mutation_position))#687; and only 23 if you use mutation_position_suppossedHg37_c; thus, we assume that CGI mutations needed to be converted to Hg38
  
  write.table(mutRed_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByCGI_df.txt",col.names = T, row.names = F, quote = F)
  
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByCGI_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mutRed_df)
  
}#(HG37)

saveMutRedByFrec_fun<-function(mut_df,nMutations2select_int=10000){
  #Saves smaller df in same format as input but with less mutation_positions. These being selected by frecuency
  #In case of tie, mutations are selected randomly
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  #nMutations2select_int: number of mutations to be selected
  
  #head(selMut_df) 
  #Sample_Barcode       mutation_position     Tumor_Subtype
  #3    TCGA-3C-AAAU    chr2_144398518           Lum
  
  #Select the nMutations2select_int most frecuent mutations from a population df
  mutationsCount_ta=sort(table(as.character(mut_df$mutation_position)),decreasing=T)
  head(mutationsCount_ta)
  #chr3_179234297  chr3_179218303  chr3_179218294 chr14_104780214   chr17_7675088  chr3_179203765 
  #117              56              38              22              18              16 
  
  # Select number of different mutation-positions that will be included in the GDL input
  selectedMutations_ch=getMutCharSelByFrec_fun(nMutations2select_int,mutationsCount_ta)
  selMut_df=mut_df[mut_df$mutation_position %in% selectedMutations_ch,]
  
  write.table(selMut_df,"/home/febueno/Documents/BRCAmetProject/BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByFrec_df.txt",col.names = T, row.names = F, quote = F)
  
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/reducedMut_df/mutRedByFrec_df.txt was saved. In this:")
  print_mutDfMainProperties_fun(selMut_df)
  
}

getMutCharSelByFrec_fun<-function(nMutations2select_int,mutationsCount_tab){
  #nMutations2select_int: number of different mutation-positions that will be included in the classifier input
  #selectedMutations_ch: Example: "chr3_179234297" "chr3_179218303"   . Mutation-positions selected. Most frecuent mutations in data are selected
  
  #1) Get threshold of mutation-positions count (i.e. number of times that this mutation appears in data) above which, all mutations will be selected. 
  #mutation-positions with threshold just below this, may have to be selected randomly to fill up to nMutations2select_int
  nMutationsThreshold_int=mutationsCount_tab[nMutations2select_int]#For nMutations2select_int=10000, mutationsThreshold_int=4
    #print(paste0("Given nMutations2select_int, the minimum mutation count for a mutation_position to remain in data was: ",nMutationsThreshold_int))
  
  #2) Select all mutation-positions above the threshold
  mutationsAboveThreshold_tab=mutationsCount_tab[mutationsCount_tab>nMutationsThreshold_int]
  
  #3) Select some mutation-positions just in the threshold
  fillMutationsChar<-function(nMutationsThreshold_int){
    #nMutationsThreshold_int: threshold of mutation-positions count (i.e. number of times that this mutation appears in data) above which, all mutations will be selected
    #mutationsInThreshold_tab: table with mutation-positions as names and mutation-count as values. all mutation-positions have nMutationsThreshold_int counts
    
    mutationsInThreshold_tab=mutationsCount_tab[mutationsCount_tab==nMutationsThreshold_int]#All mutation-positions that have counts == mutationsInThreshold_tab
    
    mutationsInThreshold_tab=sample(mutationsInThreshold_tab,nMutations2select_int-length(mutationsAboveThreshold_tab))#randomly sample from mutationsInThreshold_tab till nMutations2select_int
    return(mutationsInThreshold_tab)
  }
  
  mutationsInThreshold_tab=fillMutationsChar(nMutationsThreshold_int)
  selectedMutations_ch=c(names(mutationsAboveThreshold_tab),names(mutationsInThreshold_tab))
  
  return(selectedMutations_ch)
  
}

