
library(dplyr)
library(tidyr)
library(stringr)
library(epiR)

# 1. DATA ACQUISITION
# Load the definitive phenotypic dataset and the raw AMRFinderPlus genomic outputs
fenotipo <- read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/efc_dst_dataset_brenda_final.tsv", sep="\t", header=TRUE)
genes_raw <- read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/todos_los_genes_brenda.tsv", sep="\t", header=TRUE)

# 2. GENOTYPIC DATA NORMALIZATION
# Transform the vertical multi-class gene layout into a binary presence/absence matrix
tabla_genotipica <- genes_raw %>%
  separate_rows(Class, sep = "/") %>%
  mutate(Class = str_trim(Class)) %>%
  select(Name, Class) %>%
  distinct() %>%
  mutate(Genoma_Predice = "R") %>%
  pivot_wider(names_from = Class, values_from = Genoma_Predice, values_fill = "S") %>%
  rename(accession = Name) 

# 3. DATA INTEGRATION
# Perform a Left Join to retain all clinical isolates, including pan-susceptible strains
datos_completos <- fenotipo %>%
  left_join(tabla_genotipica, by = "accession")

# 4. CLINICAL AST MIC RESCUE 
# Helper function to strip non-numeric symbols from MIC character expressions
limpiar_mic <- function(x) { as.numeric(gsub("[<>= ]", "", x)) }

# Core function to infer categorical RIS values from raw MIC numbers when missing
rescatar_datos <- function(df, col_ris, col_mic, umbral_R, umbral_S) {
  if (col_ris %in% colnames(df) & col_mic %in% colnames(df)) {
    df[[col_ris]] <- case_when(
      !is.na(df[[col_ris]]) ~ df[[col_ris]], # Prioritize pre-existing qualitative categories
      grepl(">", df[[col_mic]]) & limpiar_mic(df[[col_mic]]) >= umbral_S ~ "R", # Right-censored adjustments
      grepl("<", df[[col_mic]]) & limpiar_mic(df[[col_mic]]) <= umbral_R ~ "S", # Left-censored adjustments
      limpiar_mic(df[[col_mic]]) >= umbral_R ~ "R",
      limpiar_mic(df[[col_mic]]) <= umbral_S ~ "S",
      TRUE ~ df[[col_ris]] 
    )
  }
  return(df)
}

# Apply localized clinical breakpoints adjusted to current diagnostic guidelines
datos_completos <- rescatar_datos(datos_completos, "vancomycin_RIS", "vancomycin_MIC", 8, 4) 
datos_completos <- rescatar_datos(datos_completos, "chloramphenicol_RIS", "chloramphenicol_MIC", 32, 8)
datos_completos <- rescatar_datos(datos_completos, "doxycycline_RIS", "doxycycline_MIC", 16, 4) 
datos_completos <- rescatar_datos(datos_completos, "quinupristin.dalfopristin_RIS", "quinupristin.dalfopristin_MIC", 4, 1)
datos_completos <- rescatar_datos(datos_completos, "tylosin_RIS", "tylosin_MIC", 32, 16)
datos_completos <- rescatar_datos(datos_completos, "gentamicin_RIS", "gentamicin_MIC", 128, 127.9)
datos_completos <- rescatar_datos(datos_completos, "streptomycin_RIS", "streptomycin_MIC", 128, 127.9)
datos_completos <- rescatar_datos(datos_completos, "kanamycin_RIS", "kanamycin_MIC", 128, 127.9)
datos_completos <- rescatar_datos(datos_completos, "erythromycin_RIS", "erythromycin_MIC", 8, 0.5) 
datos_completos <- rescatar_datos(datos_completos, "tetracycline_RIS", "tetracycline_MIC", 16, 4)
datos_completos <- rescatar_datos(datos_completos, "ciprofloxacin_RIS", "ciprofloxacin_MIC", 4, 1)

# 5. DIAGNOSTIC ACCURACY STATISTICAL EVALUATION 
calcular_precision_con_CI <- function(df, antibiotico_clinico, familia_genomica) {
  
  if (!antibiotico_clinico %in% colnames(df)) df[[antibiotico_clinico]] <- NA
  if (!familia_genomica %in% colnames(df)) df[[familia_genomica]] <- "S"
  
  # Exclude missing phenotypes to isolate the evaluable binary population
  df_eval <- df %>% filter(!is.na(!!sym(antibiotico_clinico)))
  if (nrow(df_eval) == 0) return(NULL)
  
  clinica <- df_eval[[antibiotico_clinico]]
  genoma <- df_eval[[familia_genomica]]
  genoma[is.na(genoma)] <- "S" # Genomes with no hits are designated as susceptible
  
  # Calculate contingency matrix components (Ignoring Intermediate "I" classifications)
  TP <- sum(clinica == "R" & genoma == "R", na.rm = TRUE)
  TN <- sum(clinica == "S" & genoma == "S", na.rm = TRUE)
  FP <- sum(clinica == "S" & genoma == "R", na.rm = TRUE)
  FN <- sum(clinica == "R" & genoma == "S", na.rm = TRUE)
  Total <- TP + TN + FP + FN
  
  # Compute statistical indices using the 95% Wilson Score Interval method
  stats <- epi.tests(c(TP, FP, FN, TN), conf.level = 0.95, method = "wilson")
  
  # Structure extractor to handle structure variations across package updates
  formatear_metrica <- function(metrica_nombre) {
    if ("statistic" %in% colnames(stats$detail)) {
      fila <- stats$detail[stats$detail$statistic == metrica_nombre, ]
      est <- round(fila$est * 100, 1)
      low <- round(fila$lower * 100, 1)
      upp <- round(fila$upper * 100, 1)
    } else {
      df_metrica <- as.data.frame(stats$detail[[metrica_nombre]])
      est <- round(df_metrica[1, 1] * 100, 1)
      low <- round(df_metrica[1, 2] * 100, 1)
      upp <- round(df_metrica[1, 3] * 100, 1)
    }
    return(list(est = est, ci = paste0("(", low, "-", upp, ")")))
  }
  
  sen <- formatear_metrica("se")
  spe <- formatear_metrica("sp")
  vpp <- formatear_metrica("pv.pos")
  
  return(data.frame(
    Antibiotic = antibiotico_clinico,
    Genomic_Class = familia_genomica,
    Total_Isolates = Total,
    Resistant_Isolates = TP + FN,
    Susceptible_Isolates = TN + FP,
    Sensitivity = sen$est,
    Sensitivity_95_CI = sen$ci,
    Specificity = spe$est,
    Specificity_95_CI = spe$ci,
    PPV = vpp$est,
    PPV_95_CI = vpp$ci
  ))
}

# 6. PIPELINE EXECUTION FOR THE SELECTED ANTIBIOTICS
lista_antibioticos <- list(
  c("vancomycin_RIS", "GLYCOPEPTIDE"), c("linezolid_RIS", "OXAZOLIDINONE"),
  c("erythromycin_RIS", "MACROLIDE"), c("chloramphenicol_RIS", "PHENICOL"),
  c("tetracycline_RIS", "TETRACYCLINE"), c("ciprofloxacin_RIS", "QUINOLONE"),
  c("quinupristin.dalfopristin_RIS", "STREPTOGRAMIN"), c("tylosin_RIS", "MACROLIDE"),
  c("doxycycline_RIS", "TETRACYCLINE"), c("gentamicin_RIS", "AMINOGLYCOSIDE"),
  c("streptomycin_RIS", "AMINOGLYCOSIDE"), c("kanamycin_RIS", "AMINOGLYCOSIDE")
)

# Compile results and export tidy structured dataset
tabla_final_TFG <- bind_rows(lapply(lista_antibioticos, function(x) calcular_precision_con_CI(datos_completos, x[1], x[2])))

print(tabla_final_TFG)
write.table(tabla_final_TFG, "accuracy_stats.csv", sep=";", row.names=FALSE, dec=",")