
# This R script is used to apply QC thresholds and keep run accessions that passed QC

# Qualibact cutoffs downloaded from here: 
# https://happykhan.github.io/qualibact/Enterococcus/Enterococcus_faecalis/
# https://happykhan.github.io/qualibact/methods/
# File: Enterococcus_faecalis_metrics.csv

ruta_biblio = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/búsqueda bibliográfica/results/df_cruzado_definitiva.csv"
biblio_raw = read.delim(ruta_biblio, sep = "\t", header = T)
biblio_faecalis = biblio_raw[biblio_raw$scientific_name_ena == "Enterococcus faecalis", ]
biblio_faecalis = biblio_faecalis[!is.na(biblio_faecalis$run_accession_ena) & biblio_faecalis$run_accession_ena != "", ]

input = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 1/results/df_efc.all_qc_metrics.txt"
qc_stats = read.delim(input, sep = "\t", header = T)
dim(qc_stats)

## 1. Filter Ecf 
qc_stats_faecalis = qc_stats[qc_stats$scientific_name == "Enterococcus faecalis", ]
qc_stats_faecalis = qc_stats_faecalis[!is.na(qc_stats_faecalis$run_accession), ]

# Filter NA
qc_stats = qc_stats_faecalis[!is.na(qc_stats_faecalis$Completeness_Specific), ]
qc_stats_faecalis = qc_stats_faecalis[!is.na(qc_stats_faecalis$run_accession), ]

# Filter by length and GC content respect to E.faecalis 
qc_stats_faecalis = qc_stats_faecalis[
  as.numeric(qc_stats_faecalis$GC_Content) >= 0.36 & 
    as.numeric(qc_stats_faecalis$GC_Content) <= 0.38 &
    as.numeric(qc_stats_faecalis$Genome_Size) >= 2600000 & 
    as.numeric(qc_stats_faecalis$Genome_Size) <= 3600000, 
]

print(paste("Muestras tras el filtro taxonómico estricto:", nrow(qc_stats_faecalis)))

## 2. APPLYING QUALIBACT CUTOFFS

qc_cutoffs = read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/Enterococcus_faecalis_metrics.csv", sep = ",", header = T)
dim(qc_cutoffs)
qc_cutoffs

## QUAST N50
# NOTE: there is a steep decline before 46,000 N50, we could be more strict
plot(sort(as.numeric(qc_stats$N50)), ylim = c(0,100000))
qc_stats$passed_qualibact_N50 = "no"
qc_stats$passed_qualibact_N50[which(as.numeric(qc_stats$N50)>46000)] = "yes"
table(qc_stats$passed_qualibact_N50)

## QUAST number of contigs
# NOTE: there is a steep increase after 230, we could be more strict
plot(sort(as.numeric(qc_stats$number)), ylim = c(0,1000))
qc_stats$passed_qualibact_no_of_contigs = "no"
qc_stats$passed_qualibact_no_of_contigs[which(as.numeric(qc_stats$number)<230)] = "yes"
table(qc_stats$passed_qualibact_no_of_contigs)

## checkM2 GC_Content 
# NOTE: GC_Content is a rather stable metric with 
plot(sort(as.numeric(qc_stats$GC_Content)), ylim = c(0.30,0.40))
qc_stats$passed_qualibact_gc_content = "no"
qc_stats$passed_qualibact_gc_content[which(as.numeric(qc_stats$GC_Content)>=0.36 & as.numeric(qc_stats$GC_Content)<=0.38)] = "yes"
table(qc_stats$passed_qualibact_gc_content)

## checkM2 Completeness 
plot(sort(as.numeric(qc_stats$Completeness_Specific)), ylim = c(80,100))
qc_stats$passed_qualibact_completeness = "no"
qc_stats$passed_qualibact_completeness[which(as.numeric(qc_stats$Completeness_Specific)>96)] = "yes"
table(qc_stats$passed_qualibact_completeness)

## checkM2 Contamination 
plot(sort(as.numeric(qc_stats$Contamination)), ylim = c(0,100))
qc_stats$passed_qualibact_contamination = "no"
qc_stats$passed_qualibact_contamination[which(as.numeric(qc_stats$Contamination)<3)] = "yes"
table(qc_stats$passed_qualibact_contamination)

## checkM2 Total_Coding_Sequences
plot(sort(as.numeric(qc_stats$Total_Coding_Sequences)), ylim = c(1000,5000))
qc_stats$passed_qualibact_total_coding_sequences = "no"
qc_stats$passed_qualibact_total_coding_sequences[which(as.numeric(qc_stats$Total_Coding_Sequences)>=2400 & as.numeric(qc_stats$Total_Coding_Sequences)<=3600)] = "yes"
table(qc_stats$passed_qualibact_total_coding_sequences)

## checkM2 Genome_Size
plot(sort(as.numeric(qc_stats$Genome_Size)), ylim = c(2000000,4000000))
qc_stats$passed_qualibact_total_genome_size = "no"
qc_stats$passed_qualibact_total_genome_size[which(as.numeric(qc_stats$Genome_Size)>=2600000 & as.numeric(qc_stats$Genome_Size)<=3600000)] = "yes"
table(qc_stats$passed_qualibact_total_genome_size)

## Applying all Qualibact QC cut-offs

qualibact_columns = colnames(qc_stats)[grepl("passed_qualibact_", colnames(qc_stats))]
passed_qualibact_yes = apply(qc_stats[,qualibact_columns], 1, function(x) sum(x == "yes"))
table(passed_qualibact_yes)

qc_stats$passed_qualibact_qc = "no"
qc_stats$passed_qualibact_qc[which(passed_qualibact_yes==7)] = "yes"
table(qc_stats$passed_qualibact_qc)

## 3. SPLIT FAILED AND PASSED 

dataset_aprobadas = qc_stats[qc_stats$passed_qualibact_qc == "yes", ]
dataset_suspensas = qc_stats[qc_stats$passed_qualibact_qc == "no", ]

#output_aprobadas = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 2/df_efc.all_qc_metrics.qualibact_qc.txt"
#write.table(dataset_aprobadas, file = output_aprobadas, sep = "\t", col.names = T, row.names = F, quote = F)

output_suspensas = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 2/df_efc_SUSPENSAS.txt"
write.table(dataset_suspensas, file = output_suspensas, sep = "\t", col.names = T, row.names = F, quote = F)

## 4. OUTPUT OF SAMPLES NOT ANALYSED IN QC_STATS

bacterias_ensambladas_ids = unique(qc_stats$sample_accession)
tabla_perdidas = biblio_faecalis[!(biblio_faecalis$sample_accession %in% bacterias_ensambladas_ids), ]

output_perdidas = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 2/df_efc_PERDIDAS_ENSAMBLAJE.txt"
#write.table(tabla_perdidas, file = output_perdidas, sep = "\t", col.names = T, row.names = F, quote = F)
print(paste("Tabla guardada: Muestras perdidas en ensamblaje =", nrow(tabla_perdidas)))

## 5. OUTPUT OF DUPLICATED SAMPLES

frecuencia_aprobadas = table(dataset_aprobadas$sample_accession)
ids_duplicados_ok = names(frecuencia_aprobadas[frecuencia_aprobadas > 1])

dataset_duplicados_ok = dataset_aprobadas[dataset_aprobadas$sample_accession %in% ids_duplicados_ok, ]
dataset_duplicados_ok = dataset_duplicados_ok[order(dataset_duplicados_ok$sample_accession), ]

output_duplicados_ok = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 2/df_efc_DUPLICADOS.txt"
write.table(dataset_duplicados_ok, file = output_duplicados_ok, sep = "\t", col.names = T, row.names = F, quote = F)
print(paste("Tabla guardada: Duplicados que pasaron el QC (filas) =", nrow(dataset_duplicados_ok)))

## 6. FINAL OUTPUT

# Ordenamos las aprobadas primero por muestra y luego por N50 (de mayor a menor calidad)
dataset_aprobadas_ordenado = dataset_aprobadas[order(dataset_aprobadas$sample_accession, -as.numeric(dataset_aprobadas$N50)), ]

# Nos quedamos con la primera aparición de cada bacteria (la que tiene mejor N50)
dataset_definitivo_unico = dataset_aprobadas_ordenado[!duplicated(dataset_aprobadas_ordenado$sample_accession), ]

output_definitivo = "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/control de calidad/step 2/df_efc.all_qc_metrics.qualibact_qc.txt"
write.table(dataset_definitivo_unico, file = output_definitivo, sep = "\t", col.names = T, row.names = F, quote = F)