# =========================================================================== #
# SCRIPT DE SIMPLIFICACIÓN - E. FAECALIS (VERSIÓN EXTRACCIÓN AUTOMÁTICA TOTAL)
# =========================================================================== #

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# 1. Cargar tu tabla original
y_raw <- read_excel("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/amrgen_analysis/bacterias_1088_final.xlsx", col_types = "text")
colnames(y_raw) <- trimws(colnames(y_raw))

# 2. LECTURA DE LOS ARCHIVOS DE GENES
ruta_amr <- "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/amrfinder/amrfinder_results_efc_brenda/"
archivos_genes <- list.files(ruta_amr)
ids_genes <- trimws(sub("\\..*", "", archivos_genes))

# 3. EL CRUCE INTELIGENTE DUAL-CHECK (ERR vs SAM)
y_raw <- y_raw %>%
  mutate(
    run_clean = trimws(as.character(run_accession)),
    sam_clean = trimws(as.character(sample_accession)),
    ID_DEFINITIVO = case_when(
      run_clean %in% ids_genes ~ run_clean,
      sam_clean %in% ids_genes ~ sam_clean,
      TRUE ~ run_clean
    )
  )

# 4. CONSTRUCCIÓN DE LA TABLA "TIPO FAECIUM" (Metadatos)
y_limpia <- data.frame(
  accession = y_raw$ID_DEFINITIVO,
  sample_id = y_raw$internal_id,
  biosample = y_raw$ID_DEFINITIVO,
  clade     = y_raw$host,
  study_id  = y_raw$internal_study_name,
  stringsAsFactors = FALSE
)
y_limpia$study_id[is.na(y_limpia$study_id)] <- "Brenda_TFG"

# 5. DETECCIÓN AUTOMÁTICA DE TODOS LOS ANTIBIÓTICOS EN EL EXCEL
# Buscamos todas las columnas que indican medidas o fenotipos y extraemos el nombre del fármaco
cols_claves <- grep("_(measurement|resistance_phenotype)(|\\.\\.\\.[0-9]+)$", colnames(y_raw), value = TRUE)
atbs <- unique(sub("_(measurement|resistance_phenotype).*", "", cols_claves))

cat("\n¡Detectados", length(atbs), "antibióticos automáticamente en tu Excel!\n")

# FUNCIÓN MÁGICA: Fusiona las columnas repetidas tapando huecos vacíos
fusionar_columnas <- function(df, patron) {
  cols <- grep(patron, colnames(df), value = TRUE, ignore.case = TRUE)
  if(length(cols) > 0) return(do.call(coalesce, as.list(df[cols])))
  return(rep(NA_character_, nrow(df)))
}

# 6. FUSIÓN Y EXTRACCIÓN DINÁMICA
for(atb in atbs) {
  patron_mic  <- paste0("^", atb, "_measurement(|\\.\\.\\.[0-9]+)$")
  patron_sign <- paste0("^", atb, "_measurement_sign(|\\.\\.\\.[0-9]+)$")
  patron_ris  <- paste0("^", atb, "_resistance_phenotype(|\\.\\.\\.[0-9]+)$")
  patron_meth <- paste0("^", atb, "_laboratory_typing_method(|\\.\\.\\.[0-9]+)$")
  
  # Fusionar MICs y Signos
  valores_mic <- fusionar_columnas(y_raw, patron_mic)
  signos      <- fusionar_columnas(y_raw, patron_sign)
  
  if(!all(is.na(valores_mic))) {
    signos_limpios <- ifelse(is.na(signos) | signos %in% c("=", "=="), "", trimws(signos))
    # Filtro para evitar "NANA"
    valores_finales <- ifelse(is.na(valores_mic) | trimws(valores_mic) == "NA", NA_character_, paste0(signos_limpios, trimws(valores_mic)))
    y_limpia[[paste0(atb, "_MIC")]] <- valores_finales
  }
  
  # Fusionar RIS (Fenotipos)
  valores_ris <- fusionar_columnas(y_raw, patron_ris)
  if(!all(is.na(valores_ris))) {
    y_limpia[[paste0(atb, "_RIS")]] <- case_match(valores_ris, 
                                                  "Susceptible" ~ "S", "S" ~ "S",
                                                  "Resistant" ~ "R", "Resistent" ~ "R", "R" ~ "R",
                                                  "INT" ~ "I", "I" ~ "I",
                                                  .default = NA_character_)
  }
  
  # Fusionar Métodos
  valores_meth <- fusionar_columnas(y_raw, patron_meth)
  if(!all(is.na(valores_meth))) y_limpia[[paste0(atb, "_DST_method")]] <- valores_meth
}

# 7. GUARDADO DEL ARCHIVO DEFINITIVO
write.table(y_limpia, "efc_dst_dataset_brenda_final.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=================================================================\n")
cat("¡FUSIÓN TOTAL COMPLETADA!\n")
cat("Archivo generado: efc_dst_dataset_brenda_final.tsv\n")
cat("=================================================================\n")