# EtioMR
The package "EtioMR" are used for exploring phenome-wide causal associations and mediation factors, which could significantly enhance the understanding of the etiology and pathogenesis of any disease, particularly rare diseases with unclear etiology, utilizing only GWAS summary statistics as input dataset. This package enables the phenome-wide MR analysis among health factors, biological traits, and molecular features; MR Mediation analysis for disease causing health factors; as well as sensitivity analysis including six MR methods, pleiotropy, and heterogeneity.

***(1) Prepare GWAS summary statistics files:***

The file can be in formats including comma-separated, tab-separated, or space-separated. 

The column names should be as follows:
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

***(2) Prepare packages needed:***
data.table
optparse
TwoSampleMR
MRPRESSO
dplyr
From Wednesday 1st May 2024, you need to prove your identity (authenticate) to use IEU service, even if you are querying a public dataset.
The following packages need to be updated as instructed:
install.packages("ieugwasr")
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

***(3) Run the test data:***
```
cd ./EtioMR

Rscript RunEtioMR.R \
--sumstats test.txt.gz \
--phenMR IEU#IME#MIC \
--sensitivity YES \
--mediation IME#MIC \
--IEUaccess Your IEU access \
--out result_test \
--test YES
```

***(4) Options:***

**--sumstats**: Path to summary statistics [required]

**--out**: Path to output files [required]

**--phenMR**: Choose exposure datasets:[required]
options:
IEU (2115 phenotypes from IEU-UK Biobank, including IEU analysis of UK Biobank phenotypes, Neale lab analysis of UK Biobank phenotypes round 1 and 2)
UKB (3400 phenotypes from UK Biobank)
IME (immune cell traits)
MIC (gut microbes and related pathways)
MET (metabolite levels and ratios)
UKBPPP (proteins from UK Biobank)
deCODE (proteins from deCODE genetics)
GEN (gene expression in 12 immune cells from DICE project by single-cell RNA-sequencing)
Multiple traits should be separated by # 

**--mediation**: Choose mediator datasets:
IME (immune cell traits)
MIC (gut microbes and related pathways)
Note that primary MR analysis must be performed first to conduct mediation analysis.
(i.e. --pheMR should contain IEU for all mediation analysis; --pheMR should contain IME for mediation analysis choosing IME; --pheMR should contain MIC for mediation analysis choosing MIC;)

**--IEUaccess**: This is required if you want to perform mediation analysis.
IEU access can be accessed at https://api.opengwas.io/

**--sensitivity**: Whether to conduct sensitivity analysis: YES or NO
default: NO

--test: Whether to perform test sets of exposure: YES or NO
default: NO
Be sure to include "--test YES" when using the test dataset (test.txt.gz)
The test sets of exposure are a small set of exposures including all categories in the --phenMR option. It is only used for a fast test of the package.

***(4) Citation***

