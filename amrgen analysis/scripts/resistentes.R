# -------------------------------------------------------------------------
# TABLA DE AISLAMIENTOS (R e I) CON METADATOS CLAVE
# -------------------------------------------------------------------------
cat("\n[+] Generando tabla resumen de cepas R e I...\n")

# 1. Filtramos las cepas de interés y seleccionamos solo tus 5 columnas
cepas_R_I <- pheno_table %>%
  filter(pheno %in% c("R", "I")) %>%
  select(
    biosample = id,     # El ID interno del programa viene de aquí
    sample_id,          
    accession,
    study_id = source,  # El programa renombra study_id como source
    pheno
  )

# 2. Agrupamos los genes por cada cepa
genes_agrupados <- geno_table %>%
  group_by(Name) %>%
  summarise(`genes asociados` = paste(unique(marker.label), collapse = ", "))

# 3. Cruzamos los datos usando el biosample como llave
tabla_cepas_genes <- cepas_R_I %>%
  left_join(genes_agrupados, by = c("biosample" = "Name")) %>%
  # Seleccionamos y ordenamos exactamente como has pedido
  select(
    accession,
    sample_id,
    biosample,
    study_id,
    pheno, # Añadimos el fenotipo para que sepas por qué están en esta tabla (R o I)
    `genes asociados`
  ) %>%
  arrange(pheno, accession)

# 4. Rellenamos el vacío si alguna cepa resistente no tiene genes conocidos
tabla_cepas_genes <- tabla_cepas_genes %>%
  mutate(`genes asociados` = ifelse(is.na(`genes asociados`), "Sin determinantes", `genes asociados`))

# 5. Guardamos en la carpeta y abrimos la pestaña
write_tsv(tabla_cepas_genes, "amrrules/STR1/tabla_cepas_str_metadatos.tsv")
View(tabla_cepas_genes)

cat("\n[+] ¡Tabla generada! Revisa la nueva pestaña.\n")
# -------------------------------------------------------------------------