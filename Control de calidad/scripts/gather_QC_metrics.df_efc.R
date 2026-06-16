# This R script gathers all QC metrics for dataset df_efc
# The following steps are implemented
## 1. Load run accessions and link them to ENA accessions (bioproject, biosample, etc.) and study ids
## 3. Load Quast information
## 4. Load checkM2 information

## 1. Load run accessions and link them to ENA accessions 
# Read run_accessions "assemblies" (SRR/ERR)
run_accessions1 = as.vector(as.matrix(read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/textos/df_efc.run_accessions_assemblies.paired.txt", sep = "\n", header = F)))

# Read run_accessions "ATB_assemblies" (SAMN/SAMEA)
run_accessions2_raw = as.vector(as.matrix(read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/textos/df_efc.run_accessions.atb_assemblies.txt", sep = "\n", header = F)))

# Clean extension ".fa.gz" ATB
run_accessions2 = gsub(".fa.gz", "", run_accessions2_raw)
run_accessions2 = gsub(".fa", "", run_accessions2)

run_accessions_combinados = c(run_accessions1, run_accessions2)
print(length(run_accessions_combinados))

# Read ENA
ena_accessions = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/tablas/df_efc.bioproject_accessions.ena_info.csv", sep = "\t", header = T)

#column "ID_UNICO" 
ena_accessions$ID_UNICO = ena_accessions$run_accession 
for(id in run_accessions2){
  idx = which(ena_accessions$sample_accession == id | ena_accessions$secondary_sample_accession == id)
  if(length(idx) > 0){
    ena_accessions$ID_UNICO[idx] = id
  }
}

tmp = match(run_accessions_combinados, ena_accessions$ID_UNICO)
print(paste("Faltan en ENA tras la traducción automática:", length(which(is.na(tmp)))))

qc_table2 = ena_accessions

## 3. Load Quast information

quast_results1 = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/tablas/df_efc_assemblies.quast.csv", sep = "\t", header = T)
quast_results2 = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/tablas/df_efc.atb_assemblies.quast.csv", sep = "\t", header = T)

quast_results2$sample = gsub(".fa.gz", "", quast_results2$sample)
quast_results2$sample = gsub(".fa", "", quast_results2$sample)

quast_results_combinados = rbind(quast_results1, quast_results2)

# Make sure all run accessions present in quast_results
tmp = match(ena_accessions$ID_UNICO, quast_results_combinados$sample)
missing_quast_stats = ena_accessions$ID_UNICO[which(is.na(tmp) & ena_accessions$instrument_platform=="ILLUMINA")]
write.table(missing_quast_stats, file = "df_efc.missing_quast_stats.run_accessions.txt", sep = "\n", col.names = F, row.names = F, quote = F)

# Delete duplicates
quast_results_combinados2 = quast_results_combinados[-which(duplicated(quast_results_combinados$sample)),]
if(length(which(duplicated(quast_results_combinados$sample))) == 0) { quast_results_combinados2 = quast_results_combinados }

# Pegar las tablas usando nuestro adaptador ID_UNICO
qc_table3 = merge(qc_table2, quast_results_combinados2, by.x = "ID_UNICO", by.y = "sample", all.x = T)


## 4. Load checkM2 information

checkm2_results1 = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/tablas/df_efc_assemblies.checkm2.csv", sep = "\t", header = T)
checkm2_results2 = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/procesamiento datos/assemblies/tablas/df_efc.atb_assemblies.checkm2.csv", sep = "\t", header = T)

# Limpiar los nombres en la tabla CheckM2 de ATB
checkm2_results2$Name = gsub(".fa.gz", "", checkm2_results2$Name)
checkm2_results2$Name = gsub(".fa", "", checkm2_results2$Name)

checkm2_results_combinados = rbind(checkm2_results1, checkm2_results2)

# Make sure all run accessions present in checkM2
tmp = match(ena_accessions$ID_UNICO, checkm2_results_combinados$Name)
missing_checkm2_stats = ena_accessions$ID_UNICO[which(is.na(tmp) & ena_accessions$instrument_platform=="ILLUMINA")]
write.table(missing_checkm2_stats, file = "df_efc.missing_checkm2_stats.run_accessions.txt", sep = "\n", col.names = F, row.names = F, quote = F)

# Eliminar genomas duplicados si los hay
checkm2_results_combinados2 = checkm2_results_combinados[-which(duplicated(checkm2_results_combinados$Name)),]
if(length(which(duplicated(checkm2_results_combinados$Name))) == 0) { checkm2_results_combinados2 = checkm2_results_combinados }

# Pegar las tablas usando nuestro adaptador ID_UNICO
qc_table4 = merge(qc_table3, checkm2_results_combinados2, by.x = "ID_UNICO", by.y = "Name", all.x = T)

## KEEP FINAL FILES

#Borramos la columna adaptadora ID_UNICO para que la tabla final quede super limpia
qc_table4$ID_UNICO = NULL

output = "df_efc.all_qc_metrics.txt" 
write.table(qc_table4, file = output, sep = "\t", col.names = T, row.names = F, quote = F)

# Guardar IDs de los estudios
study_ids = sort(unique(unlist(strsplit(unique(qc_table4$study_accession),";"))))
output2 = "df_efc.studies.txt"
write.table(study_ids, file = output2, sep = "\n", col.names = F, row.names = F, quote = F)
