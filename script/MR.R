
RUNMRtest <- function(out_Fin, inputrange,outfile,inputsensi) {

  if ("IEU" %in% inputrange) {
    file0<-list.files("./datatest/",pattern = "^IEU*")
  }else{
    file0<-c()
  }
  
  if ("UKB" %in% inputrange) {
    file1<-list.files("./datatest/",pattern = "^ukb*")
  }else{
    file1<-c()
  }
  
  
  if ("IME" %in% inputrange) {
    file2<-list.files("./datatest/",pattern = "^immune_clump*")
  }else{
    file2<-c()
  }  
  
  if ("MIC" %in% inputrange) {
    file3<-list.files("./datatest/",pattern = "^gut_clump*")
  }else{
    file3<-c()
  } 
  
  if ("MET" %in% inputrange) {
    file4<-list.files("./datatest/",pattern = "^metbolite*")
  }else{
    file4<-c()
  } 
  
  if ("UKBPPP" %in% inputrange) {
    file5<-list.files("./datatest/",pattern = "^ppp*")
  }else{
    file5<-c()
  } 
  
  if ("deCODE" %in% inputrange) {
    file6<-list.files("./datatest/",pattern = "^decode*")
  }else{
    file6<-c()
  } 
  
  if ("GEN" %in% inputrange) {
    file7<-list.files("./datatest/",pattern = "^sciImmiune*")
  }else{
    file7<-c()
  } 
  
  fileall<-c(file0,file1,file2,file3,file4,file5,file6,file7)
  
  sf_all<-data.frame()
  
  print(paste0("total number of exposure:",length(fileall)))
  
  


for (i in 1:length(fileall)) {
  
  met_exp<-fread(paste0("./datatest/",fileall[i])) 
  
  print(paste0("performing MR of ", fileall[i]))
  
  try_out<- try(out <- format_data(
    dat = out_Fin,
    type = "outcome",
    snps = met_exp$SNP,
    header = TRUE,
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos",
    eaf_col = "af_alt",
    ncase_col = "ncase",
    samplesize_col = "samplesize",
    phenotype_col = "Phenotype"
  ),silent=T)
  if('try-error' %in% class(try_out)) {
    print(paste0("error in formatting data",try_out," "))
    next
  }
  
  
  mydata <- harmonise_data(
    exposure_dat=met_exp,
    outcome_dat=out,
    action= 2
  )
  mydata <- subset(mydata, palindromic == FALSE) #去除回文
  if(nrow(mydata) == 0){
    print("harm null")
    next}
  
  sf <- steiger_filtering(mydata)
  
  if (inputsensi == "YES") {
    
    if(file.exists(paste0("./",outfile,"/steiger")) == TRUE){
    }else{
      dir.create(paste0("./",outfile,"/steiger"), recursive=TRUE) 
    }
    
    sf_one_file = paste0("./",outfile,"/steiger/",fileall[i], "#steiger.csv")
    write.csv(sf, file = sf_one_file,row.names = F)
  }
   
  
  sf_filtered<-data.frame(subset(sf, (steiger_dir == TRUE & mr_keep == TRUE)))
  
  if (nrow(sf_filtered)==0) {
    print("steiger null")
    next
  }
  
  sf_res <- mr(sf_filtered)
  
  sf_res<-generate_odds_ratios(sf_res)
  sf_res$b<-as.numeric(sf_res$b)
  sf_res$se<-as.numeric(sf_res$se)
  sf_res$pval<-as.numeric(sf_res$pval)
  sf_res$b_se = paste0(round(sf_res$b,3),' (',round(sf_res$b - 1.96 * sf_res$se, 3),', ', round(sf_res$b + 1.96 * sf_res$se,3),')')
  sf_res$or_ci95 = paste0(round(sf_res$or,3),' (',round(sf_res$or_lci95,3),', ',round(sf_res$or_uci95,3),')')
  sf_res$steigerfilename <- sf_one_file

  sf_ivw<-sf_res[sf_res$method=="Inverse variance weighted"|sf_res$method=="Wald ratio",]
  
  sf_all<-rbind(sf_all,sf_ivw)

}

  write.csv(sf_all,file=paste0("./",outfile,"/result_pheMR.csv"),row.names=F)
  print("MR DONE")

}




RUNMR <- function(out_Fin, inputrange,outfile,inputsensi) {

  if ("IEU" %in% inputrange) {
    file0<-list.files("./dataall/",pattern = "^IEU*")
  }else{
    file0<-c()
  }
  
  if ("UKB" %in% inputrange) {
    file1<-list.files("./dataall/",pattern = "^ukb*")
  }else{
    file1<-c()
  }
  
  
  if ("IME" %in% inputrange) {
    file2<-list.files("./dataall/",pattern = "^immune_clump*")
  }else{
    file2<-c()
  }  
  
  if ("MIC" %in% inputrange) {
    file3<-list.files("./dataall/",pattern = "^gut_clump*")
  }else{
    file3<-c()
  } 
  
  if ("MET" %in% inputrange) {
    file4<-list.files("./dataall/",pattern = "^metbolite*")
  }else{
    file4<-c()
  } 
  
  if ("UKBPPP" %in% inputrange) {
    file5<-list.files("./dataall/",pattern = "^ppp*")
  }else{
    file5<-c()
  } 
  
  if ("deCODE" %in% inputrange) {
    file6<-list.files("./dataall/",pattern = "^decode*")
  }else{
    file6<-c()
  } 
  
  if ("GEN" %in% inputrange) {
    file7<-list.files("./dataall/",pattern = "^sciImmiune*")
  }else{
    file7<-c()
  } 
  
  fileall<-c(file0,file1,file2,file3,file4,file5,file6,file7)
  
  sf_all<-data.frame()
  
  print(paste0("total number of exposure:",length(fileall)))
  
  


for (i in 1:length(fileall)) {
  
  met_exp<-fread(paste0("./dataall/",fileall[i])) 
  
  print(paste0("performing MR of ", fileall[i]))
  
  try_out<- try(out <- format_data(
    dat = out_Fin,
    type = "outcome",
    snps = met_exp$SNP,
    header = TRUE,
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos",
    eaf_col = "af_alt",
    ncase_col = "ncase",
    samplesize_col = "samplesize",
    phenotype_col = "Phenotype"
  ),silent=T)
  if('try-error' %in% class(try_out)) {
    print(paste0("error in formatting data",try_out," "))
    next
  }
  
  
  mydata <- harmonise_data(
    exposure_dat=met_exp,
    outcome_dat=out,
    action= 2
  )
  mydata <- subset(mydata, palindromic == FALSE) #去除回文
  if(nrow(mydata) == 0){
    print("harm null")
    next}
  
  sf <- steiger_filtering(mydata)
  
  if (inputsensi == "YES") {
    
    if(file.exists(paste0("./",outfile,"/steiger")) == TRUE){
    }else{
      dir.create(paste0("./",outfile,"/steiger"), recursive=TRUE) 
    }
    
    sf_one_file = paste0("./",outfile,"/steiger/",fileall[i], "#steiger.csv")
    write.csv(sf, file = sf_one_file,row.names = F)
  }
   
  
  sf_filtered<-data.frame(subset(sf, (steiger_dir == TRUE & mr_keep == TRUE)))
  
  if (nrow(sf_filtered)==0) {
    print("steiger null")
    next
  }
  
  sf_res <- mr(sf_filtered)
  
  sf_res<-generate_odds_ratios(sf_res)
  sf_res$b<-as.numeric(sf_res$b)
  sf_res$se<-as.numeric(sf_res$se)
  sf_res$pval<-as.numeric(sf_res$pval)
  sf_res$b_se = paste0(round(sf_res$b,3),' (',round(sf_res$b - 1.96 * sf_res$se, 3),', ', round(sf_res$b + 1.96 * sf_res$se,3),')')
  sf_res$or_ci95 = paste0(round(sf_res$or,3),' (',round(sf_res$or_lci95,3),', ',round(sf_res$or_uci95,3),')')
  sf_res$steigerfilename <- sf_one_file

  sf_ivw<-sf_res[sf_res$method=="Inverse variance weighted"|sf_res$method=="Wald ratio",]
  
  sf_all<-rbind(sf_all,sf_ivw)

}

  write.csv(sf_all,file=paste0("./",outfile,"/result_pheMR.csv"),row.names=F)
  print("MR DONE")

}