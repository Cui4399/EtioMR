RUNSENSI <- function(outfile,inputmediation) {
  
  MRres<-fread(paste0("./",outfile,"/result_pheMR.csv"))
  MRres<-MRres[MRres$pval<0.05,]
  
  print(paste0("Total number of traits passing primary MR analysis (P < 0.05): ",length(MRres$steigerfilename)))
  
  if (length(MRres$steigerfilename)==0) {
  print("none of the MR results passed the threshold of 0.05")
}
  
  list_steiger_files = MRres$steigerfilename
  
  
  
  sf_res_all<-data.frame()
  for(i in 1:length(list_steiger_files)){
    
    print(paste0("performing sensitivity analysis for ",list_steiger_files[i]))
    
    sf = fread(list_steiger_files[i])
    
    expname<-sf$exposure[1]
    outname<-sf$outcome[1]
    
    
    # steiger_res
    sf_filtered<-data.frame(subset(sf, (steiger_dir == TRUE & mr_keep == TRUE)))
    if(nrow(sf_filtered) == 0){
      print("steiger null")
      next}
    
    
    
    sf_res <- mr(sf_filtered)
    
    #other MR methods
    
    n = length(sf_filtered$beta.exposure)
    if(1==1){
      
      is_raps_res = 1 
      raps_inf <- try(
        raps_data <- mr.raps::mr.raps(sf_filtered$beta.exposure,sf_filtered$beta.outcome,sf_filtered$se.exposure,sf_filtered$se.outcome) ,
        silent = TRUE
      )
      if('try-error' %in% class(raps_inf)){  
        is_raps_res = 0 
      }else{
        raps_res <- c(expname,outname,outname,expname,'RAPS', n, raps_data$beta.hat,raps_data$beta.se,raps_data$beta.p.value)
        sf_res <- as.data.frame(rbind(sf_res,raps_res))
      }
      
      print('raps')
      
    }
    
    
    
    ##sensitivity
    # pleio
    sf_res <- sf_res[order(match(sf_res$method, c("Inverse variance weighted","Wald ratio","Weighted median","RAPS","MR Egger",""))), ]
    
    #my_result = data.frame(subset(sf_res, (method == 'Inverse variance weighted' | method == 'Wald ratio')))
    steiger_res_rows = nrow(sf_res)
    pleio <- mr_pleiotropy_test(sf_filtered)
    if( length(pleio) == 0){
      sf_res$egger_intercept <- ''
    }else{
      sf_res$egger_intercept <-  c(paste0(round(pleio$egger_intercept,3),'(', round(pleio$pval,3), ')'), rep('', steiger_res_rows - 1))
    }
    #mr_presso  
    is_presso = TRUE 
    presso_inf <- try(
      presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = sf_filtered, NbDistribution = 1000,  SignifThreshold = 0.05),
      silent = TRUE )
    
    if('try-error' %in% class(presso_inf)) 
    {
      print(presso_inf)
      is_presso = FALSE
    }
    
    if( is_presso == TRUE){
      
      if(presso$`MR-PRESSO results`$`Global Test`$Pvalue == '<0.001'){
        sf_res$MRPRESSO_global_test <- c(paste0(round(presso$`MR-PRESSO results`$`Global Test`$RSSobs,3),
                                                '(', presso$`MR-PRESSO results`$`Global Test`$Pvalue, ')'),
                                         rep('', steiger_res_rows - 1))
      }else{
        sf_res$MRPRESSO_global_test <- c(paste0(round(presso$`MR-PRESSO results`$`Global Test`$RSSobs,3),
                                                '(', round(presso$`MR-PRESSO results`$`Global Test`$Pvalue,3), ')'),
                                         rep('', steiger_res_rows - 1))
      }
      
      
      MRPRESSO_beta_1 = presso$`Main MR results`$`Causal Estimate`[2]
      MRPRESSO_se_1 = presso$`Main MR results`$Sd[2]
      MRPRESSO_p_1 = presso$`Main MR results`$`P-value`[2]
      MRPRESSO_beta_2 = presso$`Main MR results`$`Causal Estimate`[1]
      MRPRESSO_se_2 = presso$`Main MR results`$Sd[1]
      MRPRESSO_p_2 = presso$`Main MR results`$`P-value`[1]
      #PRESSO_res <- c(my_gene,'Mature T/NK-cell lymphomas','Mature T/NK-cell lymphomas',my_gene,'MRPRESSO', n, MRPRESSO_beta_2,MRPRESSO_se_2,MRPRESSO_p_2,"","")
      #sf_res <- as.data.frame(rbind(sf_res,PRESSO_res))
      if(is.na(MRPRESSO_beta_1) == FALSE){
        MRPRESSO_or_1 = exp(MRPRESSO_beta_1)
        MRPRESSO_or_lci95_1 = exp(MRPRESSO_beta_1 - 1.96*MRPRESSO_se_1)
        MRPRESSO_or_uci95_1 = exp(MRPRESSO_beta_1 + 1.96*MRPRESSO_se_1)
        sf_res$MRPRESSO_outlier_or = c(paste0(round(MRPRESSO_or_1,3),
                                              '(', round(MRPRESSO_or_lci95_1,3), ', ', round(MRPRESSO_or_uci95_1,3), ')'),
                                       rep('', steiger_res_rows))
        #sf_res$MRPRESSO_outlier_beta = c(MRPRESSO_beta_1,rep('', steiger_res_rows))
        sf_res$MRPRESSO_outlier_p = c(MRPRESSO_p_1,rep('', steiger_res_rows)) 
        
      }else{
        sf_res$MRPRESSO_outlier_or = ''
        sf_res$MRPRESSO_outlier_p = ''
      }
      
      
    }else{
      
      sf_res$MRPRESSO_global_test <- NA
      sf_res$MRPRESSO_outlier_or <- ''
      sf_res$MRPRESSO_outlier_p <- ''
    }
    
    
    # rucker
    steiger_res_rows<-nrow(sf_res)
    is_rucker = TRUE
    rucker_inf <- try( rucker <- mr_rucker(sf_filtered), silent = TRUE )
    if('try-error' %in% class(rucker_inf)) {
      print(rucker_inf)
      is_rucker = FALSE
    }
    
    print(length(rucker))
    
    if(is_rucker == TRUE & length(rucker) == 0){
      print("rucker 0")
      is_rucker = FALSE}
    
    if(is_rucker == TRUE){
      
      rucker_Q <- rucker[[1]]$Q
      rucker_Q_ivw <- subset(rucker_Q, Method == 'Q_ivw')
      rucker_Q_egger <- subset(rucker_Q, Method == 'Q_egger')
      rucker_Q_diff <- subset(rucker_Q, Method == 'Q_diff')
      
      my_i2 = (rucker_Q_ivw$Q - rucker_Q_ivw$df) / rucker_Q_ivw$Q
      if(my_i2 < 0) { my_i2 = 0}
      
      sf_res$I2 <- c(my_i2,rep('', steiger_res_rows - 1))
      sf_res$q_ivw <- c(paste0(round(rucker_Q_ivw$Q,3),'(', round(rucker_Q_ivw$P,3), ')'), rep('', steiger_res_rows - 1))
      sf_res$q_egger <- c(paste0(round(rucker_Q_egger$Q,3),'(', round(rucker_Q_egger$P,3), ')'), rep('', steiger_res_rows - 1))
      sf_res$q_diff <- c(paste0(round(rucker_Q_diff$Q,3),'(', round(rucker_Q_diff$P,3), ')'), rep('', steiger_res_rows - 1))
      
    }else{
      
      sf_res$I2 <- NA
      sf_res$q_ivw <- ''
      sf_res$q_egger <- ''
      sf_res$q_diff <- ''
    }

    
    sf_res$b<-as.numeric(sf_res$b)
    sf_res$se<-as.numeric(sf_res$se)
    sf_res$pval<-as.numeric(sf_res$pval)
    sf_res<-generate_odds_ratios(sf_res)
    sf_res$b_se = paste0(round(sf_res$b,3),' (',round(sf_res$b - 1.96 * sf_res$se, 3),', ', round(sf_res$b + 1.96 * sf_res$se,3),')')
    sf_res$or_ci = paste0(round(sf_res$or,3),' (',round(sf_res$or_lci95, 3),', ', round(sf_res$or_uci95,3),')')
    sf_res_all <- rbind(sf_res_all,sf_res)
    
  }
  write.csv(sf_res_all,file = paste0("./",outfile,"/result_sensitivity.csv"),row.names = F)
  
  if ("NO" %in% inputmediation) {
    unlink(paste0("./",outfile,"/steiger/"), recursive = TRUE)
  }else{}
  
  
  
  
  
}