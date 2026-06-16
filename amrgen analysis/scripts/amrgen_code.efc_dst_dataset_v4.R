## Installing and loading latest version of R packages
# See umbrella github page: https://github.com/AMRverse

# IMPORTANT NOTE: because R functions in AMRgen and AMRrulemakeR had to be manually edit to make the code run,
# the code was downloaded from GitHub, locally stored and edited and functions loaded onto this script

# NOTE: loaded R code ending in .FC.R are the ones that have been manually edited
# Only the AMR R package will be used as installed

packageVersion("AMR")

library(AMR)
# library(AMRgen)
# library(AMRrulemakeR)

path_amrgen <- "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/AMRgen_260127_FC/R/"
path_rules  <- "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/AMRrulemakeR_260127_FC/"

# Carga de archivos del tutor
source(paste0(path_amrgen, "amr_upset.R"))
source(paste0(path_amrgen, "breakpoints.R"))
source(paste0(path_amrgen, "export_ncbi.R"))
source(paste0(path_amrgen, "import_amrfp.R"))
source(paste0(path_amrgen, "re-exports.R"))
source(paste0(path_amrgen, "AMRgen-package.R"))
source(paste0(path_amrgen, "compare_geno_pheno_id.R"))
source(paste0(path_amrgen, "get_binary_matrix.R"))
source(paste0(path_amrgen, "import_pheno.R"))
source(paste0(path_amrgen, "solo_ppv.R"))
source(paste0(path_amrgen, "as.gene.R"))
source(paste0(path_amrgen, "data.R"))
source(paste0(path_amrgen, "gtdb.R"))
source(paste0(path_amrgen, "logistic.R"))
source(paste0(path_amrgen, "utils.R"))
# source("./AMRgen_260127_FC/R/assay_distribution.R")
source(paste0(path_amrgen, "assay_distribution.FC.R"))
source(paste0(path_amrgen, "eucast_distributions.R"))
source(paste0(path_amrgen, "hAMRonization_function_v3.R"))
source(paste0(path_amrgen, "plot_estimates.R"))
load(paste0(path_amrgen, "sysdata.rda"))

# source("./AMRrulemakeR_260127_FC/R/analysis_functions.R")
source(paste0(path_rules, "R/analysis_functions.FC.R")) 
# source("./AMRrulemakeR_260127_FC/R/functions.R")
source(paste0(path_rules, "R/functions.FC.R"))
source(paste0(path_rules, "R/interpretGeno.R"))
# source("./AMRrulemakeR_260127_FC/R/makeRules.R")
source(paste0(path_rules, "R/makeRules.FC.R"))
source(paste0(path_rules, "R/organism_codes.R"))
load(paste0(path_rules, "data/organism_codes.rda"))
load(paste0(path_rules, "data/hierarchy.rda"))
load(paste0(path_rules, "data/refgene_pubmed.rda"))
load(paste0(path_rules, "data/refgene.rda"))

library(ggplot2)
library(tidyverse)
library(patchwork)
library(rvest); # needed to run read_html() within get_eucast_mic_distribution()

#############################################################################
###              LOADING AND FORMATTING INPUT PHENOTYPE DATA              ###
#############################################################################

# See, https://github.com/AMRverse/AMRrulemakeR?tab=readme-ov-file#collate-and-format-phenotype-data

# See, https://github.com/AMRverse/AMRrulemakeR?tab=readme-ov-file#collate-and-format-phenotype-data

## Loading raw input AST data
y<-read_tsv(file.path("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/efc_dst_dataset_brenda_final.tsv"))# Datos de SIR y MIC
dim(y)

## Creating a new disk variable with disk diameters
## NOTE: in the original table disk diameter values are stored within _MIC variable/columns (and _DST_method used to differentiate MIC vs. disk diameter)
## new _disk columns are created, if _DST_method is disk_diffusion then disk diameter measurements are copied to the new _disk column, for each antibiotic
## MIC values in these cases are initiated to NA
ast_method_cols = colnames(y)[which(grepl("_DST_method",colnames(y)))]
ast_disk_cols = gsub("_DST_method", "_disk", ast_method_cols)
y[ast_disk_cols] = NA; 

for(c in 1:length(ast_method_cols)){
  ast_disk_col = gsub("_DST_method", "_disk", ast_method_cols[c]); 
  ast_mic_col = gsub("_DST_method", "_MIC", ast_method_cols[c]); 
  
  # La red de seguridad que comprueba si la columna MIC existe antes de tocarla
  if(ast_mic_col %in% colnames(y)){
    ast_disk_tmp = which(grepl("disk_diffusion", unlist(y[,ast_method_cols[c]]))==TRUE); 
    if(length(ast_disk_tmp)>0){
      y[ast_disk_tmp, ast_disk_col] = y[ast_disk_tmp, ast_mic_col]; 
      y[ast_disk_tmp, ast_mic_col] = NA; 
    }
  }
}
dim(y)

## Preparing AST data into right format

# Changing antibiotic names for combinations - needed to apply pivot_longer function later
colnames(y) = gsub("quinupristin_dalfopristin", "quinupristin-dalfopristin", colnames(y))
colnames(y) = gsub("piperacillin_tazobactam", "piperacillin-tazobactam", colnames(y))

# Removing unnecessary columns
y1 <- y %>%
  select(-contains("breakpoint"))
y2 <- y1 %>%
  select(-contains("genes"))

# Selecting (id) columns
columnas_metadatos <- c("accession", "sample_id", "biosample", "clade", "study_id")

# Transformación de la tabla para el correcto funcionamiento con AMRgen
y_largo <- y2 %>%
  mutate(across(-all_of(columnas_metadatos), as.character)) %>%
  pivot_longer(
    cols = -all_of(columnas_metadatos),
    names_to = c("drug_agent", ".value"),
    names_pattern = "^(.*)_(MIC|RIS|DST_method|disk)$" 
  ) %>%
  filter(!is.na(MIC) | !is.na(RIS) | !is.na(disk)) %>%
  distinct()

y_largo$MIC<-suppressWarnings(as.mic(y_largo$MIC)) 
y_largo$RIS<-suppressWarnings(as.sir(y_largo$RIS)) 
y_largo$pheno<-y_largo$RIS 
y_largo$mic<-y_largo$MIC
y_largo$disk = suppressWarnings(as.disk(y_largo$disk))
y_largo$clade<-as.factor(y_largo$clade)
y_largo <- select(y_largo,-MIC,-RIS) 
colnames(y_largo) = gsub("DST_method", "method", colnames(y_largo))
y_largo <- y_largo %>%
  mutate(method = case_when(
    grepl("broth.*microdilution", method, ignore.case = TRUE) ~ "Broth microdilution",
    grepl("agar.*dilution", method, ignore.case = TRUE) ~ "Agar dilution",
    grepl("vitek", method, ignore.case = TRUE) ~ "Vitek 2",
    grepl("etest", method, ignore.case = TRUE) ~ "Etest",
    TRUE ~ method # Si hay algún otro método raro, que lo deje como está
  ))
y_largo$source = y_largo$study_id 
y_largo$id = y_largo$biosample

#############################################################################
###               LOADING AND FORMATTING INPUT GENOTYPE DATA              ###
#############################################################################

z <- read_tsv("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/todos_los_genes_brenda.tsv", show_col_types = FALSE)
dim(z)

# Traducción de columnas vitales
if("Element symbol" %in% colnames(z)) { z <- z %>% rename(`Gene symbol` = `Element symbol`) }
if("Element name" %in% colnames(z))   { z <- z %>% rename(`Sequence name` = `Element name`) }
if("% Coverage of reference" %in% colnames(z) & !("% Coverage of reference sequence" %in% colnames(z))) {
  z <- z %>% rename(`% Coverage of reference sequence` = `% Coverage of reference`)
}
geno<-import_amrfp(z,"Name")
# WARNING: `Element subtype` field not present in input file, cannot create marker.label.
dim(geno)

# Labeling and removing genes as not present if not with at least 90% coverage
geno$cov <- ifelse(geno$`% Coverage of reference sequence` > 90, 
                   paste("C"),
                   paste("P"))
geno2<-geno%>%
  filter(!(geno$`variation type`=='Gene presence detected' & geno$cov == "P"))
dim(geno2)

# The step below is used to assign the value drug_agent to drug_class if the latter labeled as "Other antibacterials"
geno2 <- geno2 %>%
  mutate(drug_class = ifelse(drug_class == "Other antibacterials", drug_agent,
                             drug_class))

# COLUMNAS NECESARIAS PARA EL PAQUETE ---CAMBIO---
geno2$`Element type` <- "AMR"
geno2$`Element subtype` <- "AMR"
geno2$marker.label <- ifelse(!is.na(geno2$marker), geno2$marker, "Unknown")

#############################################################################
###                        SELECT ANTIBIOTIC OF FOCUS                     ###
#############################################################################

species = "Enterococcus faecalis"

# select and run first line for individual antibiotic, and run rest of the code after

# 1. Vancomicina 
#selected_drug_agent = "vancomycin"; selected_drug_class = c("Glycopeptides", "Glycopeptide", "GLYCOPEPTIDE", "vancomycin", "VANCOMYCIN"); file_prefix = "van"
# These variables will be set manually for specific antibiotics where SIR values need to be manually called
# NOTE: for antibiotics with not available EUCAST MIC breakpoints, the breakpoints are manually chosen
# See Supplementary Table 1 and Supplementary Table 2 in coll2024 for source of chosen breakpoints.
#mic_S=NULL; mic_R=NULL; disk_S=NULL; disk_R=NULL;
# These variables will be set manually for each antibiotic to run the amrrules_analysis() with different standards depending on availability of breakpoints
#run_EUCAST = TRUE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
# These variables will be set manually for each antibiotic depending on whether 
#use_disk=TRUE; use_mic=TRUE;

# Linezolid
#selected_drug_agent = "linezolid"; selected_drug_class = "Oxazolidinones"; file_prefix = "lnz"
#mic_S=NULL; mic_R=NULL; disk_S=NULL; disk_R=NULL;
#run_EUCAST = TRUE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
#use_disk=TRUE; use_mic=TRUE;

# Eritromicina SOLO PHENO Y CLSI
#selected_drug_agent = "erythromycin"; selected_drug_class = "Macrolides/lincosamides"; file_prefix = "ery"
#mic_S=NULL; mic_R=NULL; disk_S=NULL; disk_R=NULL;
#run_EUCAST = FALSE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
#use_disk=TRUE; use_mic=TRUE;

# Cloranfenicol SOLO PHENO Y CLSI
#selected_drug_agent = "chloramphenicol"; selected_drug_class = "PHENICOL"; file_prefix = "chl"
#mic_S=8; mic_R=32; disk_S=18; disk_R=12;
#run_EUCAST = FALSE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
#use_disk=TRUE; use_mic=TRUE;

# Tetraciclina SOLO PHENO Y CLSI
#selected_drug_agent = "tetracycline"; selected_drug_class = "TETRACYCLINE"; file_prefix = "tet"
#mic_S=NULL; mic_R=NULL; disk_S=NULL; disk_R=NULL;
#run_EUCAST = FALSE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
#use_disk=FALSE; use_mic=TRUE;

# Ciprofloxacino
#selected_drug_agent = "ciprofloxacin"; selected_drug_class = "QUINOLONE"; file_prefix = "cip"
#mic_S=NULL; mic_R=NULL; disk_S=NULL; disk_R=NULL;
#run_EUCAST = TRUE; run_CLSI = TRUE; run_ECOFF = TRUE;run_chosen = FALSE
#use_disk=FALSE; use_mic=TRUE;

# Now run rest of code regardless of antibiotic selected

cat("\n--- INICIANDO ANÁLISIS PARA:", selected_drug_agent, "---\n")

geno2 <- geno2 %>% 
  mutate(drug_class = ifelse(grepl(selected_drug_agent, drug_agent, ignore.case = TRUE), 
                             selected_drug_class[1], 
                             drug_class))
  
geno_table = subset(geno2, drug_class %in% selected_drug_class)
dim(geno_table)

# Ciprofloxacino
#geno_table <- geno2 %>%
#  filter(grepl("QUINOLONE|fluoroquinolone", drug_class, ignore.case = TRUE) | 
#           grepl("QUINOLONE|fluoroquinolone", `Subclass`, ignore.case = TRUE))
#if(nrow(geno_table) > 0) {
#  geno_table$drug_class <- "QUINOLONE"
#}

# Tetraciclina SOLO PHENO Y CLSI
#geno_table <- geno2 %>%
#  filter(grepl("TETRACYCLINE|tetracycline", drug_class, ignore.case = TRUE) | 
#           grepl("TETRACYCLINE|tetracyclines", `Subclass`, ignore.case = TRUE))
#if(nrow(geno_table) > 0) {
#  geno_table$drug_class <- "TETRACYCLINE"
#}

# Cloranfenicol SOLO PHENO Y CLSI
#geno_table <- geno2 %>%
#  filter(grepl("PHENICOL|phenicol", drug_class, ignore.case = TRUE) | 
#           grepl("PHENICOL|phenicols", `Subclass`, ignore.case = TRUE))
#if(nrow(geno_table) > 0) {
#  geno_table$drug_class <- "PHENICOL"
#}

pheno_table = subset(y_largo, drug_agent == selected_drug_agent)
pheno_table$drug_agent = as.ab(pheno_table$drug_agent)

cat("Genes filtrados por clase (", paste(selected_drug_class, collapse=" / "), "):", nrow(geno_table), "\n")
cat("Pruebas clínicas encontradas:", nrow(pheno_table), "\n")

# R decide automáticamente si hay discos o MICs para este antibiótico
use_disk <- any(!is.na(pheno_table$disk))
use_mic  <- any(!is.na(pheno_table$mic))
cat("¿Datos de MIC detectados?:", use_mic, "| ¿Datos de Disco detectados?:", use_disk, "\n")

cat("Genes filtrados por clase:", nrow(geno_table), "\n")
cat("Pruebas clínicas encontradas:", nrow(pheno_table), "\n")

# interpret MIC data using EUCAST breakpoints and ECOFFs
pheno_table <- interpret_ast(pheno_table, interpret_ecoff = TRUE, interpret_eucast = TRUE, interpret_clsi = TRUE, species = species, ab = selected_drug_agent)

# NOTE: command below used to decide which value to set of: run_EUCAST run_CLSI run_ECOFF use_disk use_mic
table(pheno_table$pheno)
table(pheno_table$pheno_eucast)
table(pheno_table$pheno_eucast_disk)
table(pheno_table$pheno_clsi)
table(pheno_table$pheno_clsi_disk)
table(pheno_table$ecoff)

if(run_chosen == TRUE){
  ab_analysis_chosen <- amrrules_analysis(geno_table=geno_table, pheno_table=pheno_table, 
                                          antibiotic=as.ab(selected_drug_agent), drug_class_list=selected_drug_class,
                                          sir_col="pheno",  
                                          ecoff_col = "ecoff",
                                          species=species,
                                          use_disk=use_disk,
                                          use_mic=use_mic,
                                          call_manual=TRUE,
                                          marker_col = "marker.label",
                                          minPPV=5, mafLogReg=5, mafUpset=5,
                                          mic_S=mic_S, mic_R=mic_R, disk_S=disk_S, disk_R=disk_R, # <-- AHORA SÍ TIENE VALORES
                                          geno_sample_col="Name", pheno_sample_col="id", 
                                          sir_provided_col="pheno",
                                          info=pheno_table %>% select(id, source, method))
}

if(run_EUCAST==TRUE){
  ab_analysis_eucast <- amrrules_analysis(geno_table=geno_table, pheno_table=pheno_table, 
                                          antibiotic=as.ab(selected_drug_agent), drug_class_list=selected_drug_class,
                                          sir_col="pheno_eucast",
                                          ecoff_col = "ecoff",
                                          species=species,
                                          use_disk=use_disk,
                                          use_mic=use_mic,
                                          call_manual=FALSE,
                                          marker_col = "marker.label",
                                          minPPV=5, mafLogReg=5, mafUpset=5,
                                          geno_sample_col="Name", pheno_sample_col="id", 
                                          sir_provided_col="pheno",
                                          info=pheno_table %>% select(id, source, method))
}

if(run_CLSI==TRUE){
  ab_analysis_clsi <- amrrules_analysis(geno_table=geno_table, pheno_table=pheno_table, 
                                        antibiotic=as.ab(selected_drug_agent), drug_class_list=selected_drug_class,
                                        sir_col="pheno_clsi",
                                        ecoff_col = "ecoff",
                                        species=species,
                                        use_disk=use_disk,
                                        use_mic=use_mic,
                                        call_manual=FALSE,
                                        marker_col = "marker.label",
                                        minPPV=5, mafLogReg=5, mafUpset=5,
                                        geno_sample_col="Name", pheno_sample_col="id", 
                                        sir_provided_col="pheno", # FC NOTE: important variable, that originally used (or relabeled later) containing SIR values without available MIC or disk data
                                        info=pheno_table %>% select(id, source, method))
}

if(run_ECOFF==TRUE){
  # FC NOTE: the user should set manually the variables use_mic and use_disk based on the availability of MIC/disk data
  # but if run_ECOFF is set to TRUE but there are no ECOFFs set, then use_mic and use_disk will need to be set to FALSE
  ecoffs_mic <- getBreakpoints(species, "EUCAST 2025", selected_drug_agent, "ECOFF") %>% filter(method=="MIC");
  ecoffs_disk <- getBreakpoints(species, "EUCAST 2025", selected_drug_agent, "ECOFF") %>% filter(method=="DISK");
  
  if(nrow(ecoffs_mic)==1){ 
    mic_S=as.mic(ecoffs_mic$breakpoint_S); mic_R=as.mic(ecoffs_mic$breakpoint_R);
  } else { use_mic = FALSE; }
  
  if(nrow(ecoffs_disk)==1){ 
    disk_S=as.disk(ecoffs_disk$breakpoint_S); disk_R=as.disk(ecoffs_disk$breakpoint_R); 
  } else { use_disk = FALSE; }
  
  ab_analysis_ecoff <- amrrules_analysis(geno_table=geno_table, pheno_table=pheno_table, 
                                         antibiotic=as.ab(selected_drug_agent), drug_class_list=selected_drug_class,
                                         sir_col="ecoff",
                                         ecoff_col = "ecoff",
                                         species=species,
                                         use_disk=use_disk,
                                         use_mic=use_mic,
                                         call_manual=FALSE, ###CAMBIO###
                                         marker_col = "marker.label",
                                         minPPV=5, mafLogReg=5, mafUpset=5,
                                         #mic_S=mic_S, mic_R=mic_R, disk_S=disk_S, disk_R=disk_R,
                                         geno_sample_col="Name", pheno_sample_col="id", 
                                         sir_provided_col="pheno", # FC NOTE: important variable, that originally used (or relabeled later) containing SIR values without available MIC or disk data
                                         info=pheno_table %>% select(id, source, method))
}
# check key output plots
if(run_EUCAST==TRUE){
  print(ab_analysis_eucast$ppv_plot)
  print(ab_analysis_eucast$ppv_plot_all) 
  print(ab_analysis_eucast$upset_mic_plot)
}

if(run_CLSI==TRUE){
  print(ab_analysis_clsi$ppv_plot)
  print(ab_analysis_clsi$ppv_plot_all) 
  print(ab_analysis_clsi$upset_mic_plot)
}

if(run_chosen == TRUE){
  print(ab_analysis_chosen$ppv_plot)
  print(ab_analysis_chosen$ppv_plot_all) 
  print(ab_analysis_chosen$upset_mic_plot)
}

# save analysis tables and plots, and generate rules using the Non-meningitis breakpoint, save output to 'amrrules/' with filenames starting with 'Ciprofloxacin'
# then use the rules to predicted phenotypes from genotypes and compare to the observed phenotypes (to help try to spot issues with input data and proposed rules)
# FC NOTE: if use_disk=T but no disk data available, function makerules() within amrrules_save() will set use_disk=F
# FC NOTE: code in function makerules() edited to comment out code filtering marker combination - to output all combinations regardless of association strength

if(run_EUCAST==TRUE){
  ab_rules <- amrrules_save(ab_analysis_eucast, dir_path="amrrules", file_prefix=paste(file_prefix,"_eucast",sep = ""), guide="EUCAST 2025", minObs=5)
}
if(run_CLSI==TRUE){
  ab_rules <- amrrules_save(ab_analysis_clsi, dir_path="amrrules", file_prefix=paste(file_prefix,"_clsi",sep = ""), guide="CLSI 2025", minObs=5)
}
if(run_ECOFF==TRUE){
  ab_rules <- amrrules_save(ab_analysis_ecoff, dir_path="amrrules", file_prefix=paste(file_prefix,"_ecoff",sep = ""), guide="ECOFF 2025", minObs=5,
                            mic_S=mic_S, mic_R=mic_R, disk_S=disk_S, disk_R=disk_R, use_mic=use_mic, use_disk=use_disk)
}

if(run_chosen == TRUE){
  # Guardado general para cualquier antibiótico con run_chosen = TRUE
  ab_rules <- amrrules_save(ab_analysis_chosen, dir_path="amrrules", file_prefix=paste(file_prefix,"_chosen",sep = ""), guide="Chosen Custom", minObs=5,
                            mic_S=mic_S, mic_R=mic_R, disk_S=disk_S, disk_R=disk_R, use_mic=use_mic, use_disk=use_disk)
}


# view the proposed rules, in AMRrules specification format, with quantitative fields added
if(exists("ab_rules")) { View(ab_rules$rules) }

cat("\n--- 1. FENOTIPO ---\n")
table(pheno_table$pheno, useNA = "always")

cat("\n--- 2. EUCAST 2025 ---\n")
table(pheno_table$pheno_eucast, useNA = "always")

cat("\n--- 3. ECOFF ---\n")
table(pheno_table$ecoff, useNA = "always")

cat("\n--- 4. CLSI 2025 ---\n")
table(pheno_table$pheno_clsi, useNA = "always")
