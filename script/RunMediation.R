
RUNMediation <- function(out_Fin,outfile,inputmediation,inputaccess) {
  Sys.setenv(OPENGWAS_JWT=inputaccess)
  ieugwasr::get_opengwas_jwt()
  ieugwasr::user()
  ref<-fread(paste0("./",outfile,"/result_pheMR.csv"))
  ref<-ref[ref$pval<0.05,]
  print(paste0("Number of significant phenotypes in primary MR analysis: ",nrow(ref)))
  
  ref$shortname<-sub("./result_test/steiger/","",ref$steigerfilename)
  
  ukb<-ref[grepl("^IEU", ref$shortname), ]
  
  if (nrow(ukb)==0) {
    print("No siginificant phenotypes in IEU, there is no need to perform mediation analysis")
    q()
  }
  
  if ("MIC" %in% inputmediation) {
    bac<-ref[grepl("^gut_", ref$shortname), ]
  }else{
    bac<-data.frame()
  }
  
  if ("IME" %in% inputmediation) {
    cell<-ref[grepl("^immune_", ref$shortname), ]
  }else{
    cell<-data.frame()
  }
  
  baccell<-rbind(bac,cell)
  
  print("Performing mediation analysis (univariable analysis)")
  sf_all<-data.frame()
  for (i in 1:nrow(baccell)) {
    
    bacfile<-baccell[i]
    bac_id<-bacfile$id.exposure
    bac_pheno<-bacfile$exposure
    
    print(paste0("i:",i,", id: ",bac_id))
    
    steigername_baccell<-bacfile$steigerfilename
    
    for (j in 1:nrow(ukb)) {
      
      ukb_file<-ukb[j]
      ukb_id<-ukb_file$id.exposure
      
      print(paste0("i:",i,",id: ",bac_id,"-j:",j, ", id: ",ukb_id))
      
      
      steigername_ukb<-ukb_file$steigerfilename
      steigername2_ukb<-steigername_ukb
      steigername_ukb<-sub("result_test/steiger","dataall",steigername_ukb)
      steigername_ukb<-sub("#steiger.csv","",steigername_ukb)
      ukb_data<-fread(steigername_ukb)
      
      try_out<- try(out_data <- extract_outcome_data(
        snps=ukb_data$SNP,
        outcomes=bac_id,
        proxies = FALSE,
        maf_threshold = 0.01),silent=T)
      if('try-error' %in% class(try_out)) {
        print(try_out)
        try_out<- try(out_data <- extract_outcome_data(
          snps=ukb_data$SNP,
          outcomes=bac_id,
          proxies = FALSE,
          maf_threshold = 0.01),silent=T)
      }
      if('try-error' %in% class(try_out)) {
        print(try_out)
        next
      }
      
      if (is.null(out_data)) {
        next
      }
      
      out_data$id.outcome <- bac_id
      out_data$outcome <- bac_pheno
      
      
      
      mydata <- harmonise_data(
        exposure_dat=ukb_data,
        outcome_dat=out_data,
        action= 2
      )
      mydata <- subset(mydata, palindromic == FALSE) #去除回文
      if(nrow(mydata) == 0){
        print("harm null")
        next}
      
      sf <- steiger_filtering(mydata)
      sf_filtered<-data.frame(subset(sf, (steiger_dir == TRUE & mr_keep == TRUE)))
      if (nrow(sf_filtered)==0) {
        print("steiger null")
        next
      }
      
      
      
      sf_res <- mr(sf_filtered)
      sf_res<-sf_res[sf_res$method %in% "Inverse variance weighted" | sf_res$method %in% "Wald ratio",]
      sf_res$ukbindex<-ukb_id
      sf_res$bacindex<-bac_id
      sf_res$source_baccell<-steigername_baccell
      sf_res$source_ukb<-steigername2_ukb
      sf_all<-rbind(sf_all,sf_res)
    }
    
  }
  
  out_data<-NULL
  
  TSMR_longres_all<-data.frame()
  sf_all<-sf_all[sf_all$pval<0.05,]
  sf_all<-generate_odds_ratios(sf_all)
  sf_all$bci<-paste0(round(sf_all$b,3)," (",round(sf_all$lo_ci,3),", ",round(sf_all$up_ci,3),")")
  sf_all$orci <- paste0(round(sf_all$or,3), " (", round(sf_all$or_lci95, 3), ", ", round(sf_all$or_uci95, 3), ")")

  print(paste0("Number of potential mediators: ",nrow(sf_all)))
  
  print("Performing mediation analysis (two step MR analysis)")
  
  for (i in 1:nrow(sf_all)) {
    print(i)
    ##beta1 exp(ukb) to med(baccell)
    row<-sf_all[i,]
    beta1<-row$b
    SE1<-row$se
    pval1<-row$pval
    bci1<-row$bci
    orci1<-row$orci
    longid<-paste0(row$ukbindex,"#",row$bacindex)
    ##bac snp 
    index_bac<-row$bacindex
    trait_bac<-row$outcome
    bac_steiger<-fread(row$source_baccell)
    ##ukb snp 
    index_ukb<-row$ukbindex
    trait_ukb<-row$exposure
    ukb_steiger<-fread(row$source_ukb)
    ##snp all
    snpall<-c(bac_steiger$SNP,ukb_steiger$SNP)
    snpall<-unique(snpall)
    
    ##beta_total
    totalrow<-ukb[ukb$steigerfilename==row$source_ukb,]
    beta_total<-totalrow$b
    SE_total<-totalrow$se
    pval_total<-totalrow$pval
    bci_total<-totalrow$b_se
    orci_total<-totalrow$or_ci95
    
    ##ukb data all
    try_ukb<- try(
      ukb_data<- extract_outcome_data(
        snps=snpall,
        outcomes=index_ukb,
        proxies = TRUE,
        maf_threshold = 0.01)
      ,silent = TRUE )
    if('try-error' %in% class(try_ukb)) {
      print(try_ukb)
    }else{
      ukb_data <- ukb_data
    }
    
    #while (is.null(ukb_data)) {
    #  try_ukb <- try(
    #    ukb_data <- extract_outcome_data(
    #      snps = snpall,
    #      outcomes = index_ukb,
    #      proxies = TRUE,
    #      maf_threshold = 0.01
    #    ),
    #    silent = TRUE
    #  )
    #  if ('try-error' %in% class(try_ukb)) {
    #    print(try_ukb)
    #  }
    #}
    
    if (is.null(ukb_data)) {
        next
      }
    ##bac data all
    try_exp<- try(
      bac_data<- extract_outcome_data(
        snps=snpall,
        outcomes=index_bac,
        proxies = TRUE,
        maf_threshold = 0.01)
      ,silent = TRUE )
    if('try-error' %in% class(try_exp)) {
      print(try_exp)
    }else{
      bac_data <- bac_data
    }
    
    #while (is.null(bac_data)) {
    #  try_exp <- try(
    #    bac_data <- extract_outcome_data(
    #      snps = snpall,
    #      outcomes = index_bac,
    #      proxies = TRUE,
    #      maf_threshold = 0.01
    #    ),
    #    silent = TRUE
    #  )
    #  if ('try-error' %in% class(try_exp)) {
    #    print(try_exp)
    #  }
    #}
    
    if (is.null(bac_data)) {
        next
      }
    ##ukb data need & bac data need -- > exposure
    data_snp_ins <- intersect(bac_data$SNP,ukb_data$SNP) 
    bac_data<-bac_data[bac_data$SNP %in% data_snp_ins,]
    ukb_data<-ukb_data[ukb_data$SNP %in% data_snp_ins,]
    
    exp_dat_bac <- format_data( 
      bac_data,
      type='exposure',
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      effect_allele_col ="effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      pval_col = "pval.outcome",
      samplesize_col = "samplesize.outcome",
      chr_col = "chr",
      pos_col = "pos",
      eaf_col = "eaf.outcome"
    )
    exp_dat_bac$id.exposure <- index_bac
    exp_dat_bac$exposure <- trait_bac
    
    
    exp_dat_ukb <- format_data( 
      ukb_data,
      type='exposure',   
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      effect_allele_col ="effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      pval_col = "pval.outcome",
      samplesize_col = "samplesize.outcome",
      chr_col = "chr",
      pos_col = "pos",
      eaf_col = "eaf.outcome"
    )
    
    exp_dat_ukb$id.exposure <- index_ukb
    exp_dat_ukb$exposure <- trait_ukb
    
    exp_dat_ukb <- exp_dat_ukb[order(exp_dat_ukb$SNP),]
    exp_dat_bac <- exp_dat_bac[order(exp_dat_bac$SNP),]
    exposure <- rbind(exp_dat_ukb, exp_dat_bac)
    lymdat<-out_Fin[out_Fin$rsids %in% data_snp_ins,]
    ##outcome 
    tf <- format_data(
      dat=lymdat,
      type = "outcome",
      snps = data_snp_ins,
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
    )
    
    if(is.null(tf) == TRUE){next}
    if(is.null(nrow(tf)) == TRUE){next}
    if(nrow(tf) == 0){next} 
    
    bac_data<-NULL
    ukb_data<-NULL
    
    tf <- tf[order(tf$SNP),]
    mvdat <- mv_harmonise_data(
      exposure_dat= exposure,
      outcome_dat=tf 
    )
    
    #TSMR
    
    TwoSampleMR_res <- mv_multiple(mvdat, pval_threshold = 1E+1)
    TwoSampleMR_res_df<-as.data.frame(TwoSampleMR_res)
    TwoSampleMR_res_df<-rename(TwoSampleMR_res_df,b=result.b,se=result.se)
    TwoSampleMR_res_df<-generate_odds_ratios(TwoSampleMR_res_df)
    TwoSampleMR_res_df$bci<-paste0(round(TwoSampleMR_res_df$b,3)," (",round(TwoSampleMR_res_df$lo_ci,3),", ",round(TwoSampleMR_res_df$up_ci,3),")")
    TwoSampleMR_res_df$orci <- paste0(round(TwoSampleMR_res_df$or,3), " (", round(TwoSampleMR_res_df$or_lci95, 3), ", ", round(TwoSampleMR_res_df$or_uci95, 3), ")")
    
    
    selected_rows_index <- which(TwoSampleMR_res_df$result.exposure %in% row$outcome)
    TwoSampleMR_res_df$index[selected_rows_index] <- index_bac
    TwoSampleMR_res_df$class[selected_rows_index] <- "medi"
    
    selected_rows_index <- which(TwoSampleMR_res_df$result.exposure %in% row$exposure)
    TwoSampleMR_res_df$index[selected_rows_index] <- index_ukb
    TwoSampleMR_res_df$class[selected_rows_index] <- "exp"
    TwoSampleMR_res_df$longid<-longid
    
    #TSMR long result
    
    row2<-TwoSampleMR_res_df[TwoSampleMR_res_df$class=="medi",]
    beta2<-row2$b
    SE2<-row2$se
    pval2<-row2$result.pval
    bci2<-row2$bci
    orci2<-row2$orci
    mediator<-row2$result.exposure
    
    rowdirect<-TwoSampleMR_res_df[TwoSampleMR_res_df$class=="exp",]
    beta_direct<-rowdirect$b
    SE_direct<-rowdirect$se
    pval_direct<-rowdirect$result.pval
    bci_direct<-rowdirect$bci
    orci_direct<-rowdirect$orci
    exposure<-rowdirect$result.exposure
    
    mediation_percentile<-(beta1*beta2)/beta_total
    mediation_percentile <- sprintf("%.2f%%", mediation_percentile * 100) 
    
    TSMR_longres<-data.frame(exposure,mediator,
                             beta_direct,SE_direct,pval_direct,bci_direct,orci_direct,
                             beta1,SE1,pval1,bci1,orci1,
                             beta2,SE2,pval2,bci2,orci2,
                             beta_total,SE_total,pval_total,bci_total,orci_total,
                             mediation_percentile,longid
    )
    
    TSMR_longres_all<-rbind(TSMR_longres_all,TSMR_longres)
    
    
    
  }
  unlink(paste0("./",outfile,"/steiger/"), recursive = TRUE)
  fwrite(TSMR_longres_all,file=paste0("./",outfile,"/result_mediationMR.csv"),row.names = F)
}

