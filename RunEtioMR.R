suppressMessages(require(data.table))
suppressMessages(library("optparse"))
suppressMessages(library(TwoSampleMR))
suppressMessages(source("./script/MR.R"))
suppressMessages(source("./script/RUNSENSI.R"))
suppressMessages(source("./script/RunMediation.R"))
suppressMessages(library(MRPRESSO))
suppressMessages(library(dplyr))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--phenMR", action="store", default= NA, type='character', 
              help="Choose exposure datasets: IEU, UKB, IME, MIC, MET, UKBPPP, deCODE, GEN. Multiple traits should be separated by # [required]"),
  make_option("--sensitivity", action="store", default="NO", type='character',
              help="Whether to conduct sensitivity analysis: YES or NO, default:NO"),
  make_option("--mediation", action="store", default= "NO", type='character', 
              help="Choose mediators: IME, MIC. Multiple traits should be separated by #"),
  make_option("--IEUaccess", action="store", default= "NO", type='character', 
              help="This is required if you want to perform mediation analysis"),
  make_option("--test", action="store", default="NO", type='character',
              help="Whether to perform test sets of exposure: YES or NO, default:NO")
)

opt = parse_args(OptionParser(option_list=option_list))


if ( is.na(opt$sumstats)) {
  cat("ERROR : Can not find sumstats file")
  q()
}

if ( is.na(opt$out)) {
  cat("ERROR : Can not find out dir")
  q()
}

if ( is.na(opt$phenMR)) {
  cat("ERROR : Can not find out dir")
  q()
}




out_Fin = fread(opt$sumstats)

out_Fin<-as.data.frame(out_Fin)

inputsensi<-opt$sensitivity
inputrange<-strsplit(opt$phenMR,"#")[[1]]
inputmediation<-strsplit(opt$mediation,"#")[[1]]
inputaccess<-opt$IEUaccess

nameall<-c("IEU","UKB", "IME", "MIC", "MET", "UKBPPP", "deCODE", "GEN")
for (i in 1:length(inputrange)) {
  name<-inputrange[i]
  if (!(name %in% nameall)) {
    cat(paste0("Unknown phenotype: ",name,"\n"))
    q()
  }
}

nameall<-c("IME", "MIC","NO")
for (i in 1:length(inputmediation)) {
  name<-inputmediation[i]
  if (!(name %in% nameall)) {
    cat(paste0("Unknown phenotype: ",name,"\n"))
    q()
  }
}

if("IME" %in% inputmediation){
 if(!("IEU" %in% inputrange)){
   print("Can not perform mediation analysis without primary MR analysis of IEU")
   q()
   }
 if(!("IME" %in% inputrange)){
   print("Can not perform IME mediation analysis without primary MR analysis of IME")
   q()
   }
}

if("MIC" %in% inputmediation){
 if(!("IEU" %in% inputrange)){
   print("Can not perform mediation analysis without primary MR analysis of IEU")
   q()
   }
 if(!("MIC" %in% inputrange)){
   print("Can not perform MIC mediation analysis without primary MR analysis of MIC")
   q()
   }
}

if(inputmediation != "NO"){
  if(inputaccess == "NO"){
    print("Can not perform mediation analysis without IEU access")
    q()
  }
}

if(file.exists(opt$out) == TRUE){
  }else{
    dir.create(opt$out, recursive=TRUE) 
}

###runMR 
if (opt$test == "YES") {
  RUNMRtest(out_Fin,inputrange,opt$out,inputsensi)
}else{
  RUNMR(out_Fin,inputrange,opt$out,inputsensi)
}

###run sensitivity analysis
if (inputsensi == "YES") {
  RUNSENSI(opt$out,inputmediation)
}else if(inputsensi == "NO"){

}else{
  print("Invalid option for --sensitivity")
}

###run mediation analysis
if ("NO" %in% inputmediation){
}else {
  RUNMediation(out_Fin,opt$out,inputmediation,inputaccess)
}


