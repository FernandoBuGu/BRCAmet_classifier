#Project: Deep Learning Classifier to estimate the probability of future metastasis from Breast Primary Tumor WES data
#@feBueno fernando.bueno.gutierrez@usal.es August-2020

#This file contains the main functions for BRCAmet_getInput.R

#REQUIRES (SCRIPTS)
#BRCAmet_funs2.R


#REQUIRE (PACKAGES): 
    library(maftools)
#Auxiliary functions for these functions or for getMetClassiInput.R are in:
  #Get primary-tumor-breast-cancer samples and mutation data from different data sources. 1 file per data source:
    source("./BRCAmet/scripts/inputSamplesFuns/getGDCinput_funs.R")
    source("./BRCAmet/scripts/inputSamplesFuns/getAKT1input_funs.R")
    source("./BRCAmet/scripts/inputSamplesFuns/getV8input_funs.R")
  #Select mutations subsets based on different approaches
    source("./BRCAmet/scripts/BRCAmet_funs2.R")

createTrainTestSets_fun<-function(mut_df,testSizePosClass_i=100,testSizeNegClass_i=100){
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  #dMBclass_c: character that is either "Yes" or "No"
  
  #testSizeClassX_i: 100, number of samples to move to test set for one of the classes (indicated in the dMBclass_c argument)
  
  #trainTestMut_l: list with two df: $train and $test
  
  randomlySampleTestSamplesFrom1class_fun<-function(mut_df,dMBclass_c,testSizeClassX_i=100){
    
    #head(mut_df)
    # Sample_Barcode mutation_position dMB source
    # 1   TCGA-A8-A093     chr1_32363805  No    GDC
    # 2   TCGA-A8-A093     chr1_52033476  No    GDC
    
    #dMBclass_c: character that is either "Yes" or "No"
    
    #testSizeClassX_i: 100
    
    # head(sampleNames2bTestFromdMBclass_c)
    # [1] "MBC-MBCProject_LvS2IvIy-Tumor-SM-DL3IT" "GENIE-MSK-P-0001149-T01-IM3"
    
    oneRowPerSample_df=mut_df[!duplicated(mut_df$Sample_Barcode),]
    selSamples_idx=which(oneRowPerSample_df$dMB==dMBclass_c)
    randomSelTestSamples_idx=sample(selSamples_idx,testSizeClassX_i)
    allSampleNames_c=as.character(unique(mut_df$Sample_Barcode))
    sampleNames2bTestFromdMBclass_c=allSampleNames_c[randomSelTestSamples_idx]
    
    return(sampleNames2bTestFromdMBclass_c)
    
  }  
  
  sampleNames2bPosTest_c=randomlySampleTestSamplesFrom1class_fun(mut_df,dMBclass_c="Yes",testSizeClassX_i=testSizePosClass_i)
  sampleNames2bNegTest_c=randomlySampleTestSamplesFrom1class_fun(mut_df,dMBclass_c="No",testSizeClassX_i=testSizeNegClass_i)
  sampleNames2bPosNegTest_c=append(sampleNames2bPosTest_c,sampleNames2bNegTest_c)
  
  testMut_df=mut_df[mut_df$Sample_Barcode %in% sampleNames2bPosNegTest_c,]
  trainMut_df=mut_df[!mut_df$Sample_Barcode %in% sampleNames2bPosNegTest_c,]
  
  trainTestMut_l=list("train"=trainMut_df,"test"=testMut_df)
  
  return(trainTestMut_l)
}
      
hg19toHg38_fun<-function(mutation_position_hg37_c){
  #mutation_position_hg37_c: "chr1_115256529","chr17_7577018","ThisOneIsWrong" 
  #mutation_position_hg38_c: "chr1_114713908" "chr17_7673700"  NA 
  
    #will print: "InputNAs: 0;  OutputNAs: 1"
  
  #Converts human genome reference
  
  #WARNING!
  #requires pasting http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz in /home/febueno/R/x86_64-pc-linux-gnu-library/3.6/liftOver/extdata
  
  #As explained in https://www.rdocumentation.org/packages/rtracklayer/versions/1.32.1/topics/liftOver
  library(rtracklayer)
  path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
  chain <- import.chain(path)

  library(doSNOW)
  library(foreach)
  cores <- parallel::detectCores()-2
  cores
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=8000, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  mp_hg19toHg38_fun<-function(mutation_position_c,chain=chain){#convert mutation position 
    
    #mutation_position_c: character, chr1_115256529 that corresponds to hg37 (also called hg19 or GRCh37)
    #mutation_position_c: character, chr1_114713908 that corresponds to hg38 (GRCh38)
    
    #Convertion in the example can be checked in https://genome.ucsc.edu/cgi-bin/hgLiftOver 
    #1)Original Assembly: hg19
    #2)New Assembly : hg38
    #3)Paste: chr1:115256529-115256529
    #4)submit button
    #5)view conversion button
    
    #parse input
    mchr=sub("_.*", "", mutation_position_c)#mchr: my chromosome
    mstart=as.numeric(sub(".*_", "", mutation_position_c))#mstart: my start
    
    
    df <- data.frame(chr=mchr, start=mstart, end=mstart)
    
    create_GRanges<-function(df){
      dfranges=makeGRangesFromDataFrame(df)#https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/makeGRangesFromDataFrame.html  
      converted=unlist(liftOver(dfranges, chain))
      chr=converted@seqnames@values
      pos=converted@ranges@start
      mutation_position_c=paste0(chr,"_",pos)
      return(mutation_position_c)
    }
    
    mutation_position_c=try(create_GRanges(df),silent = T)
    if(class(mutation_position_c)=="try-error"){
      mutation_position_c=NA
    }
    
    return(mutation_position_c)
    
  }  
  
  
  mutation_position_hg38_l=foreach(i=1:length(mutation_position_hg37_c), .options.snow=opts, .packages = c("biomaRt","rtracklayer")) %dopar% mp_hg19toHg38_fun(mutation_position_c=mutation_position_hg37_c[i],chain=chain)
  close(pb)
  stopCluster(cl)
  mutation_position_hg38_c=unlist(mutation_position_hg38_l)
  
  #mutation_position_hg38_c<-unlist(lapply(mutation_position_hg37_c,mp_hg19toHg38_fun,chain=chain))
  
  NAinput=sum(is.na(mutation_position_hg37_c))
  NAoutput=sum(is.na(mutation_position_hg38_c))
  print(paste0("InputNAs: ",NAinput,";  OutputNAs: ",NAoutput))
  
  return(mutation_position_hg38_c)
}

print_mutDfMainProperties_fun<-function(mut_df){
  #Given a (reduced?) mutations_dataframe prints number of samples and features  
  
  # head(mut_df)
  # Hugo_Symbol Sample_Barcode mutation_position dMB source
  # 1       TSSK3   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2       KTI12   TCGA-A8-A093     chr1_52033476  No    GDC
  
  mut_df$dMB_source=paste0(mut_df$dMB,"_",mut_df$source)
  
  print(paste0("number of unique mutation_position: ",length(unique(mut_df$mutation_position))))
  print(paste0("number of unique samples: ",length(unique(mut_df$Sample_Barcode))))
  print("number of samples that developed metastasis:")
  print(table(mut_df$dMB[!duplicated(mut_df$Sample_Barcode)]))
  print("number of samples (positives+negatives) from each data set:")
  print(table(mut_df$source[!duplicated(mut_df$Sample_Barcode)]))
  print("number of positive and negative samples from each data set:")
  print(table(mut_df$dMB_source[!duplicated(mut_df$Sample_Barcode)]))
  print(paste0("number of 1s in the sampleXfeature matrix: ",nrow(mut_df)))
}

saveMutationData_fun<-function(excludeNgeV8_b=F){
  
  #Saves a file with the following format in "./BRCAmet/data/WES/myTransformedData/ready_df/mut_df.txt"
  # head(mut_df)#mut_df is a combination of X mutation_position files of similar format and should contain data from all samples top be considered by metClassi.py
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC

  #excludeNgeV8_b: F, wether the large number of negative V8 should be excluded, since they are not very reliable. For instance, I found there some samples from the AKT1 cohort
  
  #WARNING! uncomment the following lines if mutation data is beeing created for the first time or if it was modified
  # saveReadyGDCdata()
  # saveReadyAKT1data()
  # saveReadyV8data()
  
  GDCmut_df=read.table("./BRCAmet/data/WES/myTransformedData/ready_df/GDCmut_df.txt",header=T)
  AKT1mut_df=read.table("./BRCAmet/data/WES/myTransformedData/ready_df/AKT1mut_hg38_df.txt",header=T)
  V8mut_df=read.table("./BRCAmet/data/WES/myTransformedData/ready_df/V8mut_hg38_df.txt",header=T)
  MBCPmut_df=read.table("./BRCAmet/data/WES/myTransformedData/ready_df/MBCPmut_HG38_df.txt",header=T)
  
  mergedSampleData_df=rbind(GDCmut_df,AKT1mut_df)
  mergedSampleData_df=rbind(mergedSampleData_df,V8mut_df)
  mut_df=unique(rbind(mergedSampleData_df,MBCPmut_df))
  
  mut_df=mut_df[!nchar(as.character(mut_df$mutation_position))<5,]#remove wierd name mutations
  
  if(excludeNgeV8_b){
    mut_df=mut_df[!(mut_df$dMB=="No" & mut_df$source=="V8"),]
  }
  
  write.table(mut_df,"./BRCAmet/data/WES/myTransformedData/ready_df/mut_hg38_df.txt",col.names = T, row.names = F, quote = F)
  
  print("-------------")
  print("./BRCAmet/data/WES/myTransformedData/ready_df/mut_hg38_df.txt was saved. Main properties:")
  print_mutDfMainProperties_fun(mut_df)
  
}

mafPath2sampleMp_fun<-function(path2mafDownload){
  #path2mafDownload: string with path to .maf file. Example: "/home/febueno/Documents/mgdl4/GDCdata/TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf"
  #si38mp_df (si38mp: sample-id and GRCh38-mutation-position) example: 
    #Sample_Barcode  mutation_position
    #TCGA-BH-A0HO      chr1_152760163
  
  MAF=read.maf(path2mafDownload)
  df=as.data.frame(MAF@data)#There is a 1-to-1 correspondance bt Matched_Norm_Sample_Barcode and Tumor_Sample_Barcode. With 979 unique values in each case
  si38mp_df=df[,c("Tumor_Sample_Barcode","Chromosome","Start_Position","Hugo_Symbol")]
  
  #get Sample_Barcode
    #Tumor_Sample_Barcode: Example: "TCGA-BH-A0HO-01A-11W-A050-09"
    #Tumor_Sample_Barcode need to be shortened after the 3th "-" in their name, in order to match the Case.ID in GDCsampleSheet. 
    #This becomes simply "Sample_Barcode" with a unique Barcode for every sample from a case_id (Example: 8332806e-f547-4aae-89af-6d5bec831fd2) 
    #This corresponds, for instance to a Tumor_Sample_Barcode and a Matched_Norm_Sample_Barcode.
  si38mp_df$Sample_Barcode=do.call(paste, c(read.table(text = as.character(si38mp_df$Tumor_Sample_Barcode), sep = "-")[1:3], sep = "-"))
  si38mp_df$Tumor_Sample_Barcode<-NULL
  
  #get mutation_position
  si38mp_df$mutation_position<-paste(si38mp_df$Chromosome,si38mp_df$Start_Position,sep="_")#define mutation_position. Example: chr1_152760163
  
  #get si38mp_df
  si38mp_df=si38mp_df[ , -which(names(si38mp_df) %in% c("Chromosome","Start_Position"))]#Exclude "Chromosome","Start_Position" since that info is already in "mutation_position"
  return(si38mp_df)
  }

saveClassiInput_fun<-function(mutRed_df,percVal_n=0.2,percTest_n=0.2){
  
  # head(mutRed_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr3_13819294  No    GDC
  # 2   TCGA-A8-A093    chr3_179218306  No    GDC
  
  #percVal_n: percentage of samples that will be in the validation set. Used only for MetClass input
  
  #percTest_n: percentage of samples that will be in the test set. Used only for GDL input 
  
  #NOTE: with current version of code, you either select a validation set (for MetClass), or a test set (for GDL)
  
  #1)remove source column, no longer required
  mutRed_noSource_df=mutRed_df[,c(1:3)]
  
  #2)create bsm_df: binary sample mutations matrix
  createPredictorsBinary_fun<-function(mutRed_noSource_df){
    
    # head(mutRed_noSource_df)
    # Sample_Barcode mutation_position Tumor_Subtype
    # 1   TCGA-A8-A093     chr3_13819294            No
    # 2   TCGA-A8-A093    chr3_179218306            No
    
    # bsm_df[1:2,1:3]
    # chr3_13819294 chr3_179218306 chr4_39062858
    # TCGA-A8-A093             1              1             1
    # TCGA-C8-A133             0              0             0
    
    uniqueSample_Barcode_c=as.character(unique(mutRed_noSource_df$Sample_Barcode))
    getBinaryMatrix_fun<-function(sample_idx){
      
      #sample_idx: 1
      
      #head(sample_idx_mutations_bi_n)
      #[1] 0 0 1 0 0 0
      
      sampleMp_c=mutRed_noSource_df$mutation_position[mutRed_noSource_df$Sample_Barcode == uniqueSample_Barcode_c[sample_idx]]
      sample_idx_mutations_bi_n=rep(0,length(unique(mutRed_noSource_df$mutation_position)))#sample_idx_mutations_bi_n: binary numeric of length equal to the number of unique mutations in mutRed_noSource_df, that has 1 if sample_idx had that mutation, 0 otherwise
      columns2b1_idx=which(as.character(unique(mutRed_noSource_df$mutation_position)) %in% as.character(sampleMp_c))#columns2b1_idx: index of the columns that will get value 1 (mutation_position found) for sample_idx 
      sample_idx_mutations_bi_n[columns2b1_idx]<-1
      
      return(sample_idx_mutations_bi_n)
    }
    oneBinaryNumPerSample_ln=lapply(1:length(uniqueSample_Barcode_c),function(x) getBinaryMatrix_fun(x))#oneBinaryNumPerSample_ln, list of numerics of length equal to the number of samples in mutRed_noSource_df
      #head(oneBinaryNumPerSample_ln[[1]])
      # [1] 1 1 1 1 1 1
    bsm_df=as.data.frame(do.call(rbind, oneBinaryNumPerSample_ln))
    colnames(bsm_df)<-unique(mutRed_noSource_df$mutation_position)
    rownames(bsm_df)<-uniqueSample_Barcode_c
    
    return(bsm_df)
  }
  bsm_df=createPredictorsBinary_fun(mutRed_noSource_df)

  #3)save
  saveMetClassInput<-function(bsm_df,mutRed_noSource_df,percVal_n=0.2){
    
    #add label as last column
    addLabel2bsm_fun<-function(bsm_df,mutRed_noSource_df){
      
      # head(mutRed_noSource_df)
      # Sample_Barcode mutation_position dMB
      # 1   TCGA-A8-A093     chr3_13819294  No
      # 2   TCGA-A8-A093    chr3_179218306  No
      
      # bsm_df[1:2,3289:ncol(bsm_df)]
      # chr9_21974697 chr9_8485824 chr9_21974696 chr5_13810194
      # TCGA-A8-A093             0            0             0             0
      # TCGA-C8-A133             0            0             0             0
      
      # bsmLabelAsLastCol_df[1:2,3289:ncol(bsm_df)]
      # chr9_21974697 chr9_8485824 chr9_21974696 chr5_13810194 label
      # TCGA-A8-A093             0            0             0             0   0
      # TCGA-C8-A133             0            0             0             0   0
      
      mutRed_sbts_df=unique(mutRed_noSource_df[,c(1,3)]) 
      positiveSamples_idx=which(mutRed_sbts_df$dMB!="No")
      rm(mutRed_sbts_df)
      label_c=rep(0,nrow(bsm_df))
      label_c[positiveSamples_idx]<-1
      bsmLabelAsLastCol_df=bsm_df
      bsmLabelAsLastCol_df$label<-label_c
      
      return(bsmLabelAsLastCol_df)
    }
    bsmLabelAsLastCol_df=addLabel2bsm_fun(bsm_df,mutRed_noSource_df)
    
    #randomly create training and validation sets
    createRandomValSet_fun<-function(bsmLabelAsLastCol_df,percVal_n=0.2){
      
        # bsmLabelAsLastCol_df[1:2,3289:ncol(bsm_df)]
        # chr9_21974697 chr9_8485824 chr9_21974696 chr5_13810194 label
        # TCGA-A8-A093             0            0             0             0     0
        # TCGA-C8-A133             0            0             0             0     0
      
        #percVal_n=0.2
      
        #trainValDfs_l: list with two df, called train and validation
      
        #randomly shuffle rows
        allRowsShuffled_idx <- sample(nrow(bsmLabelAsLastCol_df))
        bsmLabelAsLastCol_rS_df <- bsmLabelAsLastCol_df[allRowsShuffled_idx, ]#rS row shuffled
        
        #get total number of validation samples 
        nValSamples_i=floor(nrow(bsmLabelAsLastCol_rS_df)*percVal_n)
        
        #separate negatives and positives
        negBsmLabelAsLastCol_rS_df=bsmLabelAsLastCol_rS_df[bsmLabelAsLastCol_rS_df$label==0,]
        posBsmLabelAsLastCol_rS_df=bsmLabelAsLastCol_rS_df[bsmLabelAsLastCol_rS_df$label==1,]
        
        #get number of negative and positive validation
        nNegValSamples_i=floor(nValSamples_i*nrow(negBsmLabelAsLastCol_rS_df)/(nrow(bsmLabelAsLastCol_rS_df)))
        nPosValSamples_i=floor(nValSamples_i*nrow(posBsmLabelAsLastCol_rS_df)/(nrow(bsmLabelAsLastCol_rS_df)))
        
        #divide negatives samples in train and val
        train_negBsmLabelAsLastCol_rS_df=negBsmLabelAsLastCol_rS_df[(nNegValSamples_i+1):nrow(negBsmLabelAsLastCol_rS_df),]
        val_negBsmLabelAsLastCol_rS_df=negBsmLabelAsLastCol_rS_df[1:nNegValSamples_i,]
        
        #divide positive samples in train and val
        train_posBsmLabelAsLastCol_rS_df=posBsmLabelAsLastCol_rS_df[(nPosValSamples_i+1):nrow(posBsmLabelAsLastCol_rS_df),]
        val_posBsmLabelAsLastCol_rS_df=posBsmLabelAsLastCol_rS_df[1:nPosValSamples_i,]
        
        #rbind negative and positive train samples and randomly shuffle rows
        train_bsm_df=rbind(train_negBsmLabelAsLastCol_rS_df,train_posBsmLabelAsLastCol_rS_df)
        all_train_bsm_RowsShuffled_idx <- sample(nrow(train_bsm_df))
        train_bsm_df <- train_bsm_df[all_train_bsm_RowsShuffled_idx, ]
        
        #rbind negative and positive test samples and randomly shuffle rows
        val_bsm_df=rbind(val_negBsmLabelAsLastCol_rS_df,val_posBsmLabelAsLastCol_rS_df)
        all_val_bsm_RowsShuffled_idx <- sample(nrow(val_bsm_df))
        val_bsm_df <- val_bsm_df[all_val_bsm_RowsShuffled_idx, ]
        
        trainValDfs_l=list("train"=train_bsm_df,"validation"=val_bsm_df)
        
        return(trainValDfs_l)
    }
    trainValDfs_l=createRandomValSet_fun(bsmLabelAsLastCol_df,percVal_n)
    
    #save data
    save4metClassiInputFiles_fun<-function(trainValDfs_l){
      #define the 4 data pieces
      Xtr_df=trainValDfs_l$train[,1:(ncol(trainValDfs_l$train)-1)]
      ytr_n=trainValDfs_l$train$label
      Xv_df=trainValDfs_l$validation[,1:(ncol(trainValDfs_l$validation)-1)]
      yv_n=trainValDfs_l$validation$label
      
      #WARNING! IF YOU SHUFFLE COLUMNS, THE TEST SET WILL BE SCREWED
      #randomly shuffle columns
      # allColsShuffled_idx <- sample(ncol(Xtr_df))
      # Xtr_colShuffled_df=Xtr_df[,allColsShuffled_idx]
      # Xv_colShuffled_df=Xv_df[,allColsShuffled_idx]
      Xtr_colShuffled_df=Xtr_df#WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Xv_colShuffled_df=Xv_df
      
      #full experiment-file name
      setwd("/home/febueno/Documents/BRCAmetProject")
      XtrName_c="./metClassiInputByBRCAmet/Xtr.txt"
      XvName_c="./metClassiInputByBRCAmet/Xv.txt"
      ytrName_c="./metClassiInputByBRCAmet/ytr.txt"
      yvName_c="./metClassiInputByBRCAmet/yv.txt"
      
      #save
      write.table(Xtr_colShuffled_df,XtrName_c,col.names = F, row.names = F, quote = F, sep="\t")
      write.table(Xv_colShuffled_df,XvName_c,col.names = F, row.names = F, quote = F, sep="\t")
      write.table(ytr_n,ytrName_c,col.names = F, row.names = F, quote = F, sep="\t")
      write.table(yv_n,yvName_c,col.names = F, row.names = F, quote = F, sep="\t")

      print("4 files have been saved:")
      print(XtrName_c)
      print(XvName_c)
      print(ytrName_c)
      print(yvName_c)
      
      print("dim(Xtr)")
      print(dim(Xtr_colShuffled_df))
      print("dim(Xv)")
      print(dim(Xv_colShuffled_df))
    
      print("These are all training files needed by metClassi.py")
      
      print("-----------------------------")
      print(paste0("input_shape: ",ncol(Xtr_colShuffled_df)))
      
      allLabels_c=append(ytr_n,yv_n)
      print("labels (including train and validation), are:")
      print(table(allLabels_c))
      classWeights0_n=1
      percentageOf0s_f=table(allLabels_c)[1]/(table(allLabels_c)[1]+table(allLabels_c)[2])
      percentageOf1s_f=table(allLabels_c)[2]/(table(allLabels_c)[1]+table(allLabels_c)[2])
      classWeights1_n=round(percentageOf0s_f/percentageOf1s_f,2)
      print("accodingly, class_weight should be:")
      print(paste0("0: ",classWeights0_n,"     1: ",classWeights1_n))
    }
    save4metClassiInputFiles_fun(trainValDfs_l)
  }
  saveMetClassInput(bsm_df,mutRed_noSource_df,percVal_n)
  
}  

saveTestClassiInput_fun<-function(trainMutRed_df,testMut_df,MutReductionMethod_c="test_byFrec"){
  
  # head(trainMutRed_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093    chr3_179218306  No    GDC
  # 2   TCGA-A8-A093     chr4_39062858  No    GDC
  
  # head(testMut_df)
  # Sample_Barcode mutation_position dMB source
  # 400   TCGA-BH-A0W5     chr1_45500381  No    GDC
  # 401   TCGA-BH-A0W5    chr1_237890148  No    GDC
  
  #MutReductionMethod_c="test_byFrec"
  
  #saves 2 .txt files:
  #/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/test_byFrec_Xte.txt
  #/home/febueno/Documents/BRCAmetProject/metClassiInputByBRCAmet/test_byFrec_yte.txt
  
  
  #give to test set the model input format: predictors on which it was trained
  mpInModel_c=as.character(unique(trainMutRed_df$mutation_position))
  
  #reduce test set size to neccesay
  testMutFilteredOnTrainMP_df=testMut_df[testMut_df$mutation_position %in% mpInModel_c,]
  print_mutDfMainProperties_fun(testMutFilteredOnTrainMP_df)
  
  testBmatrix_df=data.frame(matrix(0,nrow=length(unique(testMutFilteredOnTrainMP_df$Sample_Barcode)),ncol=length(mpInModel_c)))#test binary matrix
  colnames(testBmatrix_df)=mpInModel_c
  testSamplesNames_c=as.character(unique(testMutFilteredOnTrainMP_df$Sample_Barcode))
  
  #add 1s to testBmatrix_df for each test sample
  fillMatrixForTestSampleX_fun<-function(testSample_idx){
    
    #testSample_idx:1
    
    # updates testBmatrix_df  
    #testBmatrix_df[1:2,1:5]
    # chr3_179218306 chr4_39062858 chr4_143547438 chr7_95290118 chr9_128691169
    # 1              0             0              0             0              0
    # 2              0             0              0             0              0
    
    mpSampleX_c=as.character(testMutFilteredOnTrainMP_df$mutation_position[testMutFilteredOnTrainMP_df$Sample_Barcode==testSamplesNames_c[testSample_idx]])
    mpSampleX_idx=which(mpInModel_c %in% mpSampleX_c)
    testBmatrix_df[testSample_idx,mpSampleX_idx]<<-rep(1,length(mpSampleX_idx))
    
  }
  for(idx in 1:length(testSamplesNames_c)){fillMatrixForTestSampleX_fun(idx)}
  
  rowPerSampleTestMutFilteredOnTrainMP_df=testMutFilteredOnTrainMP_df[!duplicated(testMutFilteredOnTrainMP_df$Sample_Barcode),]
  yte_n<-rep(0,nrow(rowPerSampleTestMutFilteredOnTrainMP_df))
  posTestSample_idx=which(rowPerSampleTestMutFilteredOnTrainMP_df$dMB=="Yes")
  yte_n[posTestSample_idx]<-rep(1,length(posTestSample_idx))
  
  #full experiment-file name
  XteName_c=paste0("./metClassiInputByBRCAmet/",MutReductionMethod_c,"_Xte.txt")
  yteName_c=paste0("./metClassiInputByBRCAmet/",MutReductionMethod_c,"_yte.txt")
  
  #save
  write.table(testBmatrix_df,XteName_c,col.names = F, row.names = F, quote = F, sep="\t")
  write.table(yte_n,yteName_c,col.names = F, row.names = F, quote = F, sep="\t")
  
  #print
  print("2 files have been saved:")
  print(XteName_c)
  print(yteName_c)
  
  print("dim(testBmatrix_df)")
  print(dim(testBmatrix_df))
  print("dim(yte_n)")
  print(dim(yte_n))
  
  print("These are all testing files needed by metClassi.py")
  
  print("-----------------------------")
  print(paste0("input_shape: ",ncol(testBmatrix_df)))
  
  print("tets labels are:")
  print(table(yte_n))
  
}

createTrainTestSets_fun<-function(mut_df,testSizePosClass_i=100,testSizeNegClass_i=100){
  
  #head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  #testSizePosClass_i: 100, number of samples to move to test set for the positive class
  
  #testSizeNegClass_i: 100, number of samples to move to test set for the negative class
  
  #trainTestMut_l: list with two df: $train and $test
  
  randomlySampleTestSamplesFrom1class_fun<-function(mut_df,dMBclass_c,testSizeClassX_i=100){
    
    #head(mut_df)
    # Sample_Barcode mutation_position dMB source
    # 1   TCGA-A8-A093     chr1_32363805  No    GDC
    # 2   TCGA-A8-A093     chr1_52033476  No    GDC
    
    #dMBclass_c: character that is either "Yes" or "No"
    
    #testSizeClassX_i: 100
    
    # head(sampleNames2bTestFromdMBclass_c)
    # [1] "MBC-MBCProject_LvS2IvIy-Tumor-SM-DL3IT" "GENIE-MSK-P-0001149-T01-IM3"
    
    oneRowPerSample_df=mut_df[!duplicated(mut_df$Sample_Barcode),]
    selSamples_idx=which(oneRowPerSample_df$dMB==dMBclass_c)
    randomSelTestSamples_idx=sample(selSamples_idx,testSizeClassX_i)
    allSampleNames_c=as.character(unique(mut_df$Sample_Barcode))
    sampleNames2bTestFromdMBclass_c=allSampleNames_c[randomSelTestSamples_idx]
    
    return(sampleNames2bTestFromdMBclass_c)
    
  }  
  
  sampleNames2bPosTest_c=randomlySampleTestSamplesFrom1class_fun(mut_df,dMBclass_c="Yes",testSizeClassX_i=testSizePosClass_i)
  sampleNames2bNegTest_c=randomlySampleTestSamplesFrom1class_fun(mut_df,dMBclass_c="No",testSizeClassX_i=testSizeNegClass_i)
  sampleNames2bPosNegTest_c=append(sampleNames2bPosTest_c,sampleNames2bNegTest_c)
  
  testMut_df=mut_df[mut_df$Sample_Barcode %in% sampleNames2bPosNegTest_c,]
  trainMut_df=mut_df[!mut_df$Sample_Barcode %in% sampleNames2bPosNegTest_c,]
  
  trainTestMut_l=list("train"=trainMut_df,"test"=testMut_df)
  
  return(trainTestMut_l)
}

saveTrainTestMetClassiInput<-function(mut_df,remNeg2balance_b=T,nSamplesBalancedTest_i=100,percVal_n=0.2,MutReductionMethod_c="byFrec",extraParameter_i=1000){
  
  # head(mut_df)
  # Hugo_Symbol Sample_Barcode mutation_position dMB source
  # 1      ATPAF1   TCGA-BH-A18F     chr1_46665282  No    GDC
  # 2        RORC   TCGA-BH-A18F    chr1_151814660  No    GDC
  
  #remNeg2balance_b: wether negative samples should be removed in order to have balanced classes
  #nSamplesBalancedTest_i: integer, total number of samples to move to test set, half of them will be positive samples and half will be negative
  #percVal_n: percentage of samples that will be moved from the training set to the validation set. neg/pos proportions will be maintained
  #MutReductionMethod_c: one in {"byFrec", "CGI", "cosmic", "HPA", "Bailey1", "Bailey2", "Bailey3", "consensus", "groupMpInGenes", "createMpWindows", "none"}
  #extraParameter_i: Provide if MutReductionMethod_c in {"byFrec", "createMpWindows"}. Functions called in these cases require an extraparametes. See ./BRCAmet/scripts/mutReduce_funs.R
  
  #general arrangements
  testSizePosClass_i=floor(nSamplesBalancedTest_i/2)#integer of samples to move to test set for the positive class
  testSizeNegClass_i=nSamplesBalancedTest_i-testSizePosClass_i#integer of samples to move to test set for the negative class
  
  #balance?
  if(remNeg2balance_b){mut_df=balanceMut_fun(mut_df)}
  
  #crate train and test set
  trainTestMut_l=createTrainTestSets_fun(mut_df,testSizePosClass_i=testSizePosClass_i,testSizeNegClass_i=testSizeNegClass_i)
  trainMut_df=trainTestMut_l$train
  testMut_df=trainTestMut_l$test
  
  #apply mutation reduction in the training set
  applyMutRed_fun<-function(trainMut_df,MutReductionMethod_c, extraParameter_i){
    
    # head(trainMut_df)
    # Sample_Barcode mutation_position dMB source
    # 3900   TCGA-BH-A42U    chr1_220762956  No    GDC
    # 3901   TCGA-BH-A42U    chr1_228372781  No    GDC
    
    #MutReductionMethod_c: one in {"byFrec", "CGI", "cosmic", "HPA", "Bailey1", "Bailey2", "Bailey3", "consensus", "groupMpInGenes", "createMpWindows", "none"}
    
    #extraParameter_i: Provide if MutReductionMethod_c in {"byFrec", "createMpWindows"}. Functions called in these cases require an extraparameter. See ./BRCAmet/scripts/mutReduce_funs.R
    
    #head(trainMutRed_df)
    # Sample_Barcode mutation_position dMB source
    # 1   TCGA-BH-A0HN    chr3_179203765  No    GDC
    # 2   TCGA-BH-A0HN     chr6_95586980  No    GDC)
    
    switch(MutReductionMethod_c,
      "byFrec"={
        saveMutRedByFrec_fun(trainMut_df,nMutations2select_int=extraParameter_i)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByFrec_df.txt",header=T)))  
      },
      "CGI"={
        saveMutRedCGI_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByCGI_df.txt",header=T)))  
      },
      "cosmic"={
        saveMutRedCosmic_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByCensus_df.txt",header=T)))  
      },
      "HPA"={
        saveMutRedHPAprognosisGenes_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByHPAprognosisGenes_df.txt",header=T)))  
      },
      "Bailey1"={
        saveMutRedBailey1_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey1_df.txt",header=T)))  
      },
      "Bailey2"={
        saveMutRedBailey2_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey2_df.txt",header=T)))  
      },
      "Bailey3"={
        saveMutRedBailey2_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByBailey3_df.txt",header=T)))  
      },
      "consensus"={
        saveMutRedConsensus_fun(trainMut_df)
        return(unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedCon_df.txt",header=T)))  
      },
      "none"={
        return(trainMut_df)
      },
      
      #WARNING! Two of the MutReductionMethod_c, are applied on the whole dataset
      "groupMpInGenes"={
        saveGroupMpInGenes_fun(mut_df)
        mutRed_df=unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByGroupMpInGenes_df.txt",header=T))
        trainTestMut_l=createTrainTestSets_fun(mutRed_df,testSizePosClass_i=testSizePosClass_i,testSizeNegClass_i=testSizeNegClass_i)
        trainMut_df=trainTestMut_l$train
        testMut_df<<-trainTestMut_l$test
        return(trainMut_df)
      },
      "createMpWindows"={
        saveMpWindows_fun(mut_df,nWindowsPerChr_i=extraParameter_i)
        mutRed_df=unique(read.table("./BRCAmet/data/WES/myTransformedData/ready_df/mutReduced_df/mutRedByMpWindows_df.txt",header=T))
        trainTestMut_l=createTrainTestSets_fun(mutRed_df,testSizePosClass_i=testSizePosClass_i,testSizeNegClass_i=testSizeNegClass_i)
        trainMut_df=trainTestMut_l$train
        testMut_df<<-trainTestMut_l$test
        return(trainMut_df)
      }
    )
  }#WARNING! groupMpInGenes and createMpWindows, are applied on the whole dataset.
  trainMutRed_df=applyMutRed_fun(trainMut_df,MutReductionMethod_c,extraParameter_i)
  trainMutRed_df$Hugo_Symbol=NULL
  
  
  #save train set after mutation reduction
  saveClassiInput_fun(trainMutRed_df,percVal_n=percVal_n)
  
  #save test set with mutations accounted for in the model
  saveTestClassiInput_fun(trainMutRed_df,testMut_df)
}

balanceMut_fun<-function(mut_df){
  #will remove negative examples in order to have the same number of positive and negative examples
  
  # head(mut_df)
  # Sample_Barcode mutation_position dMB source
  # 1   TCGA-A8-A093     chr1_32363805  No    GDC
  # 2   TCGA-A8-A093     chr1_52033476  No    GDC
  
  # head(balancedMut_df)
  # Sample_Barcode mutation_position dMB source
  # 863   TCGA-A8-A083    chr1_154947515  No    GDC
  # 864   TCGA-A8-A083    chr1_196902547  No    GDC
  
  rowPerSampleMut_df=mut_df[!duplicated(mut_df$Sample_Barcode),]
  posSamples_idx=which(rowPerSampleMut_df$dMB=="Yes")
  negSamples_idx=which(rowPerSampleMut_df$dMB=="No")
  selNegSamples_idx=sample(negSamples_idx,min(length(posSamples_idx),length(negSamples_idx)))
  allSamples_c=as.character(rowPerSampleMut_df$Sample_Barcode)
  selSamples_c=allSamples_c[c(posSamples_idx,selNegSamples_idx)]
  balancedMut_df=mut_df[as.character(mut_df$Sample_Barcode) %in% selSamples_c,]
  
  return(balancedMut_df)
  
}








