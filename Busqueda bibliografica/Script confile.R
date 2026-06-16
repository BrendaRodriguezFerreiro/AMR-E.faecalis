# 1.  SET UP ####
install.packages(c("tidyverse", "readxl", "data.table", "httr"))
## 1.1 Cargar la librería tidyverse y readxl ####
library(tidyverse)
library(readxl)
library(data.table)
library(httr)
library(dplyr)

## 1.2 Definir la ruta base ####
#TODO! base_path Argumento de entrada al script. args[1]
#base_path <- "//nasibv05/Datos/Grupos/UAMG/rferragud/5.eco_amr_pred/collections/RFF_revised"

base_path <- "/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/búsqueda bibliográfica/articles/revised"

# Validar si la ruta base existe
if (!dir.exists(base_path)) {
  stop("La ruta base no existe o no es accesible: ", base_path)
}

# Crear archivo de log
logfile <- file.path(base_path, paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S_"), "log.txt"))

log_connection <- file(logfile, open = "a")

# Iniciar logging (a archivo y consola)
sink(log_connection, append = TRUE, split = TRUE)
sink(log_connection, append = TRUE, type = "message")

log_msg <- function(study_name, msg) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", study_name, " - ", msg, "\n")
}

# Redirigir las advertencias a la salida de log en tiempo real
options(warn = 1)

#log_msg("====== INICIO DEL SCRIPT ======")


# 2. SCRIPT ####

all_summaries <- data.frame("isolate_id_publication" = as.character())

# create Logfile

## 2.1 Crear un data.frame vacío con nombres de columnas únicos ####
summary_df <- 
  data.frame(
    isolate_id_publication = character(),
    isolate_name_alias = character(),
    internal_study_name = character(),
    bioproject_accession = character(),
    sample_accession = character(),
    secondary_sample_accession = character(),
    sample_alias = character(),
    run_accession = character(),
    geographic_location_country = character(),
    geographic_location_region = character(),
    geographic_location_location = character(),
    latitude_and_longitude = character(),
    collection_date = character(),
    collected_by = character(),
    host = character(),
    isolation_source = character(),
    host_disease = character(),
    host_tissue_sampled = character(), 
    host_subject_id = character(),
    additional_info = character()
  )

## 2.2 Obtener la lista de estudios (subdirectorios en "collections") ####

studies <- 
  list.dirs(base_path, recursive = FALSE)

for (study in studies) {
  # Validar si el directorio del estudio existe
  if (!dir.exists(study)) {
    log_msg(study_name, "⚠️ El directorio del estudio no existe:")
    next 
  } 
  
  tryCatch({
    
    #study = "C:/Users/rferragud/Documents/Projecte/RFF_revised/lang2019"
    
    #study = "//nasibv05/Datos/Grupos/UAMG/rferragud/5.eco_amr_pred/collections/RFF_revised/mancini2019"
    
    # set working directory
    setwd(study)
    
    #getwd()
    
    study_name <- basename(study)  # Extrae solo el nombre del directorio
    
    log_msg(study_name, "Inicio del estudio")
    
    ## 2.3 Buscar las rutas de los archivos de configuración (partA) ####
    
    # Buscar archivos PartA y partB usando Sys.glob para manejar regular expressions
    partA_files <- 
      Sys.glob(file.path(study, paste0("configuration_file.", study_name, "*.partA.v0.*.xlsx")))     ##CAMBIO DE CSV A XLSX
    
    partB_files <- 
      Sys.glob(file.path(study, paste0("configuration_file.", study_name, "*.partB.v0.*.xlsx")))    ##CAMBIO DE CSV A XLSX
    
    if (length(partA_files) == 0 || length(partB_files) == 0) { 
      log_msg(study_name, "⚠️ ADVERTENCIA: No se encontraron archivos PartA o PartB para el estudio")
      next } else {
        
    # 3. partA ####
    
    # Lista para almacenar todos los summary_df_partA que se generen en el bucle for 
    all_summary_df_partA <- list() 
    
      #la partB seguiria fent-se no?
      
      for (partA_file in partA_files) { # iterar sobre las partsA del estudio
        
        log_msg(study_name, paste("Leyendo partA:", basename(partA_file)))
        
        tryCatch({
          partA <- as.data.frame(read_excel(path = partA_file)) %>% # leer partA file                                   ##CAMBIO CSV A excel
            mutate(across(everything(), str_trim))  #quitar los espacios de delante y detras 
        }, error = function(e) {
          log_msg(study_name, paste("💥 ERROR leyendo PartA:", partA_file, "-", conditionMessage(e)))
          next
        })
        
        # Verificar si la celda en 'value' de la fila 'bioproject_accession' está vacía
        if (partA$value[partA$standard_field == "bioproject_accession"] == "" | 
            is.na(partA$value[partA$standard_field == "bioproject_accession"])) {
          log_msg(study_name, paste("⚠️ ADVERTENCIA La columna 'value' en la fila 'bioproject_accession' está vacía."))
        }
        
        if (partA$value[partA$standard_field == "internal_study_name"] == "" | 
            is.na(partA$value[partA$standard_field == "internal_study_name"])) {
          log_msg(study_name, paste("⚠️ ADVERTENCIA La columna 'value' en la fila 'internal_study_name' está vacía."))
        }
        
        if (partA$value[partA$standard_field == "isolate_id_publication"] == "" | 
            is.na(partA$value[partA$standard_field == "isolate_id_publication"])) {
          log_msg(study_name, paste("💥 ERROR La columna 'value' en la fila 'isolate_id_publication' está vacía."))
        }
        

        ## 3.2. Obtener las tablas suplementarias e iterar sobre estas ####
        
        # obtener número de tablas suplementarias
        supp_data_tables <- 
          partA %>% 
          filter(!source_field %in% c("missing", "not applicable", "")) %>% 
          select(source_field) %>% 
          distinct() %>% 
          pull(source_field)
        
        all_supp_data_tables <-  data.frame("isolate_id_publication" = as.character())
        
        # iterar sobre las tablas suplementarias 
        for (supp_table in supp_data_tables) {
          #log_msg(study_name, "for (supp_table in supp_data_tables)")
      
          # comprobar tipo de archivo (csv o excel)
          file_extension <- tools::file_ext(supp_table)
          
          ### 3.2.1 Extraer y leer los campos de interés en csv ####
          
          if (file_extension == "csv") {  
            
            log_msg(study_name, paste("Leyendo csv: ", supp_table))
            
            # Correspondencia standard_field/value de los campos de la supp_table (sin espacios)
            supp_table.fields2extract.temp <-  
              partA %>% 
              filter(source_field == supp_table) %>% 
              select(standard_field, value) %>%
              mutate(value = str_squish(value)) %>% #quitar espacios antes y despues de las palabras
              mutate(value = str_replace_all(value, " ", "\\.")) #quitar los puntos que queden entre las palabras 
            # Si usamos make.names quita las ; 
            
            # Obtener los nombres de las columnas 
            col2extract <- 
              supp_table.fields2extract.temp %>% 
              pull(value) %>% 
              strsplit(";") %>% 
              unlist()   %>%
              map_chr(str_trim) %>%
              str_replace_all( "^\\.", "")
            
            log_msg(study_name, "col2extract")
            
            # Leer las columnas de interés de las supp_table
            tryCatch({
              
              suptab2 <- 
                fread(supp_table, check.names = T) %>%
                select(all_of(col2extract))
              
            }, error = function(e) {
              log_msg(study_name, paste("💥 ERROR leyendo tabla suplementaria:", supp_table, "-", conditionMessage(e)))
              next
            })
          
            ## 3.2.2 Extraer y leer los campos de interés en excel ####
            
          }  else if (file_extension %in% c("xls", "xlsx", "XLSX")) { 
            
            log_msg(study_name, paste("Leyendo excel", supp_table))
            
            # comprobar cuantas tabs hay 
            supp_table.tabs <-  
              partA %>% 
              filter(source_field == supp_table) %>% 
              select(tab_number) %>% 
              mutate(tab_number = as.integer(tab_number)) %>% 
              distinct() %>% 
              pull(tab_number)
            
            isolate_id_publication.0 <- 
              partA %>%
              filter(standard_field == "isolate_id_publication")%>%
              select(standard_field, value)
            
            isolate_id_publication.1 <- 
              partA %>%
              filter(standard_field == "isolate_id_publication") %>%  
              pull(value)  # Extrae el valor
            
            #crear un df vacío con una sola columna isolate_id_publication.1 para después poder unir la info de las tabs
            suptab2 <- setNames(data.frame(as.character()), isolate_id_publication.1)

            for (tab in supp_table.tabs) {
              tab.t=tab
              # Damos por hecho que en todas las tabs hay las mismas isolate_id_publication
              # para cada tab hay un isolate_id_publication column 
              
              #log_msg(study_name, "for (tab in supp_table.tabs)")
              
              # 3.3 Obtener la tabla suplementaria donde se encuentran los sample ids ####
             
              # Extraemos primer los isolate_id_publication de todas las tablas y tabs
              isolate_id_publication.suptable.tab <- 
                partA %>% 
                filter(source_field == supp_table & tab_number == tab.t &
                        str_starts(standard_field, "isolate_id_publication")) %>% 
                pull(value)
              
              isolate_id_publication.suptable.tab
              
              # Correspondencia standard_field/value de los campos de la supp_table (con espacios)
              supp_table.fields2extract.temp <-  
                partA %>% 
                filter(source_field == supp_table & tab_number == tab.t) %>% 
                select(standard_field, value) 
              
              # Obtener los nombres de las columnas 
              col2extract <-
                supp_table.fields2extract.temp %>% 
                pull(value) %>% 
                strsplit(";") %>%
                unlist() %>% 
                map_chr(str_trim)
              
              log_msg(study_name, "col2extract")
              
              # Leer las columnas de interés de las supp_table
              tryCatch({
                suptab2.temp <- 
                  as.data.frame(read_excel(supp_table, sheet = as.integer(tab.t), col_types = "text", progress = FALSE)) %>% 
                  select(all_of(col2extract)) %>% 
                  rename_with(~ isolate_id_publication.1, !!isolate_id_publication.suptable.tab)
                }, error = function(e) {
                log_msg(study_name, paste("💥 ERROR leyendo tabla suplementaria:", supp_table, "-", conditionMessage(e)))
                next
              })
              
              # Unir las info de las tabs al df vacío
              suptab2 <- 
                full_join(
                  suptab2, suptab2.temp, 
                  by= isolate_id_publication.1)
             
              
            } #  for (tab in supp_table.tabs)
            
            #print(head(suptab2))
            
           supp_table.fields2extract.temp <-  
             partA  %>%
             filter(source_field == supp_table) %>% 
             select(standard_field, value) %>% 
             mutate(value =str_replace_all(value, "^\\.", ""))
           
           if (length(supp_data_tables) > 1){
             supp_table.fields2extract.temp <- 
               bind_rows(supp_table.fields2extract.temp, isolate_id_publication.0) %>%
               mutate(value = make.names(value))
           }
            
           #log_msg("supp_table.fields2extract.temp")
            
          } # else if (file_extension %in% c("xls", "xlsx"))
          
          ## 3.4 Obtener la correspondencia entre standard_field y value #### 
          
          # correspondencia entre standard_field y value con todos los elementos separados (sin ;)
          supp_table.fields2extract <- 
            supp_table.fields2extract.temp %>%
            separate_rows(value, sep = ";") %>%
            mutate(value = str_trim(value)) %>%
            mutate(value = str_replace_all(value, "^\\.", ""))
          
          if (any(str_detect(supp_table.fields2extract.temp$value, ";"))) {
          # Crear identificadores únicos para la columna standard_field
          supp_table.fields2extract$standard_field <- #poner en el if??
            ave(
              supp_table.fields2extract$standard_field,
              supp_table.fields2extract$standard_field,
              FUN = function(x) {
                if (length(x) > 1) {
                  paste0(x, "_", seq_along(x))
                } else {
                  x
                }
              }
            )
          
          #log_msg("supp_table.fields2extract")
          
         #print(head(supp_table.fields2extract))
          
         ## 3.5 Cambiar los nombres de las variables por los normalizados ####
          
          # Crear el vector de correspondencia 
          new_colnames <- setNames(supp_table.fields2extract$standard_field, 
                                   supp_table.fields2extract$value)
          
          # Renombrar las columnas de suptab2 usando la correspondencia
          colnames(suptab2) <- new_colnames[colnames(suptab2)]
          
          suptab2 <- 
            suptab2 %>% 
            filter(
              if_any(starts_with("isolate_id_publication"), ~ !is.na(.))
            )
          
          log_msg(study_name, "tabla normalizada")
          #print(head(suptab2))
          
          ## 3.6 Unir las columnas que aportan info de la misma variable ####
          
          #encontrar el standard_field para el cual el value tiene ; 
         # if (any(str_detect(supp_table.fields2extract.temp$value, ";"))) {
            standard_fields_with_semicolon <- supp_table.fields2extract.temp %>%
              filter(str_detect(value, ";")) %>%
              pull(standard_field)
            
            #hacer una lista que contiene los standard_fields con ; seguidos de _número
            matching_columns.lists <- map(standard_fields_with_semicolon, function(prefix) {
              names(suptab2)[str_detect(names(suptab2), paste0("^", prefix, "_\\d+$"))]
            })
            
            #poner como nombre a las listas el standard_field correspondiente
            names(matching_columns.lists) <- standard_fields_with_semicolon
            
            #bucle en el cual para cada nombre de las listas se unen las columnas 
            for (new_col_name in names(matching_columns.lists)) {
              # Extraer las columnas a unir
              cols_to_unite <- matching_columns.lists[[new_col_name]]
              # Aplicar unite
              suptab2 <- suptab2 %>%
                unite(!!new_col_name, all_of(cols_to_unite), sep = ";", remove = TRUE)
            }
            
            suptab_united <- suptab2 %>%
              mutate(across(everything(), ~ str_remove(., ";$")))
            
            tail(suptab_united)
            
          } else {
            
            # Crear el vector de correspondencia 
            new_colnames <- setNames(supp_table.fields2extract$standard_field, 
                                     supp_table.fields2extract$value)
            
            # Renombrar las columnas de suptab2 usando la correspondencia
            colnames(suptab2) <- new_colnames[colnames(suptab2)]
            
            suptab2 <- 
              suptab2 %>% 
              filter(
                if_any(starts_with("isolate_id_publication"), ~ !is.na(.))
              )
            
            
            suptab_united <- suptab2
          }
          
          #log_msg("suptab_united")
          #print(head(suptab_united))
          
          ## 3.7 Unir la info obtenida al summary_df ####
          
          # Obtener los nombres de las columnas de ambos dataframes
          suptab_united.columns <- colnames(suptab_united)
          summary_df.columns <- colnames(summary_df)
          
          # Setdiff indica las columnas de summary_df_columns que no están en suptab2_columns
          missing_columns <- setdiff(summary_df.columns, suptab_united.columns)
          
          # Crear un dataframe con las columnas faltantes en suptab2, llenando con NA (numero de filas igual que suptab2)
          missing_data <- data.frame(matrix(NA, ncol = length(missing_columns), nrow = nrow(suptab_united)))
          colnames(missing_data) <- missing_columns
          
          # Unir los datos de suptab2 con los missing_data
          combined_df <- bind_cols(suptab_united, missing_data)
          
          # Poner las columnas ordenadas como en summary_df
          summary_df_suptab <- combined_df %>%
            select(all_of(summary_df.columns))
          
          # Imprimir el dataframe final
          
          #print("summary_df_suptab")
          #print(head(summary_df_suptab))
          
          ## 3.8 Unir valores del manuscript al summary_df ####
          
          manuscript_fields <-  
            partA %>% 
            filter(source_type == "manuscript") %>% 
            select(standard_field, value)
        
          for (i in 1:nrow(manuscript_fields)) {
            # 1. Extraemos el nombre de la columna del dataframe 'summary_df' que queremos actualizar.
            column_to_update <- manuscript_fields$standard_field[i]
            
            # 2. Extraemos el valor que queremos asignar a dicha columna.
            value_to_assign <- manuscript_fields$value[i]
            
            # 3. Actualizamos la columna correspondiente en 'summary_df'.
            summary_df_suptab[[column_to_update]] <- ifelse(
              is.na(summary_df_suptab[[column_to_update]]),  # Si el valor en la columna es NA...
              value_to_assign,                              # Asignamos el nuevo valor.
              summary_df_suptab[[column_to_update]]          # Si no, mantenemos el valor existente.
            )
          }
          
          log_msg(study_name, "manuscript_fields unidos a la tabla normalizada")
          #print(head(summary_df_suptab))
          
          summary_df_suptab <- summary_df_suptab %>%
            mutate(isolate_id_publication = as.character(isolate_id_publication))
          
          all_supp_data_tables2 <-
            full_join(all_supp_data_tables, summary_df_suptab,
                  by = "isolate_id_publication")
          
          } #supp_table
          
          # Add study name
          summary_df_partA <- 
            all_supp_data_tables2 %>% 
            mutate(internal_study_name = !!study_name)
          
          rm(summary_df_suptab)
          
          # Añadir a la lista los df con partA 
          all_summary_df_partA <- append(all_summary_df_partA, list(summary_df_partA))
          
          log_msg(study_name, "all_summary_df_partA")
          
      } # for (partA_file in partA_files)
    
    # Unir todos los summary_df_partA que hay en la lista all_summary_df_partA en un solo data frame
    final_summary_df_partA <- bind_rows(all_summary_df_partA)
    
    # Imprimir la tabla combinada
    #print(head(final_summary_df_partA))
    
    #view(final_summary_df_partA)
    
    
    # 4. Part B ====================================================================
    
    ## 4.1 Iterar sobre todas las partsB del estudio ####
    
    #Lista para almacenar todos los summary_df_partB que se generen en el bucle for 
    
    all_summary_df_partB <- list()
    
    for (partB_file in partB_files) { # iterar sobre las partsB del estudio
      
      log_msg(study_name, paste("Leyendo partB:", basename(partB_file)))
      
      tryCatch({
        partB <- as.data.frame(read_excel(path = partB_file)) %>% # leer partB file                                 ##CAMBIO CSV A XLSX
          mutate(across(everything(), str_trim))%>% #quitar los espacios de delante y detras 
          filter(standard_subfield != "internal_study_name" )
      }, error = function(e) {
        log_msg(study_name, paste("💥 ERROR leyendo PartB:", partA_file, "-", conditionMessage(e)))
        next
      })
      
        ## 4.2 Limpiar el archivo de configuración ####
        
        #Eliminar las filas que no nos dan info ("missing") 
        partB_tidy <- 
          partB %>% 
          filter((standard_field == "to all" & source_type != "missing") | (standard_field != "to all")) %>%
          mutate(across(everything(), str_trim)) %>% 
          mutate(standard_field = tolower(standard_field))
        
        #print(head(partB_tidy))
        
        atb_name_list <- 
          partB_tidy %>% 
          select(standard_field) %>% 
          filter(standard_field != "to all") %>% 
          distinct()
        
        valid_names <- read.delim("/Users/brendaaracelirodriguezferreiro/Library/CloudStorage/OneDrive-UPV/TFG/2026_Efc/búsqueda bibliográfica/articles/revised/antibiotics.txt", sep = "\t") %>%
          pull(1) 
        
        # Verificar cuáles no coinciden
        invalid_names <- atb_name_list %>%
          filter(!standard_field %in% valid_names) %>% 
          distinct() %>% 
          pull(standard_field)
        
        # Lanzar warning si hay nombres inválidos
        if (length(invalid_names) > 0) {
          log_msg(study_name, paste("⚠️ ADVERTENCIA Se encontraron valores inválidos: ", 
                  paste(invalid_names, collapse = ", ") 
                  ))
        }
        
        
        ## 4.3 Duplicar las filas "to all" para todos los antibióticos ####
        
        # Filas base "to all"
        to_all_rows <- partB_tidy %>% filter(standard_field == "to all")
        
        # Filas con antibióticos
        antibiotic_rows <- partB_tidy %>% filter(standard_field != "to all")
        
        # Reemplazar las rows con la info especifica de cada atb
        to_all_rows.standard_subfield <-
        to_all_rows %>% 
          select(standard_subfield) %>% 
          distinct() %>% 
          pull(standard_subfield)

        atb_to_remove <-
        antibiotic_rows %>% 
          filter(standard_subfield %in% to_all_rows.standard_subfield)
        
        # Expandir las filas "to all" para cada antibiótico
        new_rows <- antibiotic_rows %>% 
          distinct(standard_field) %>% 
          cross_join(to_all_rows) %>% 
          mutate(standard_field = standard_field.x) %>% 
          select(-standard_field.x, -standard_field.y) %>% 
        anti_join(atb_to_remove, by = c("standard_field", "standard_subfield"))
      
        
        # Combinar los "to all" de los antibióticos con el df original (tiene info propia del atb) 
        partB_tidy_atb <- bind_rows(partB_tidy, new_rows) %>% 
          filter(standard_field != "to all") %>% 
          arrange(standard_field)
        
        #print(head(partB_tidy_atb))
        
        #unir standard_field con standard_subfield por _ 
        partB_tidy_atb.suptab <- partB_tidy_atb %>%
          unite(standard_subfield_atb, standard_field, standard_subfield, sep="_")
        
        log_msg(study_name, "información to all combinada con los antibióticos")
        
        ## 4.4 Extraer la información de la tabla suplementaria  ####
        
        # Obtener la tabla suplementaria de la partB 
        suptabB <- 
          partB_tidy_atb %>% 
          filter(!source_field %in% c("missing", "not applicable", "")) %>% 
          select(source_field) %>% 
          distinct() %>% 
          pull(source_field)
        
        # comprobar tipo de archivo (csv o excel)
        file_extension <- tools::file_ext(suptabB)
        
        ### 4.4.1 Extraer y leer los campos de interés en csv ####
        
        if (file_extension == "csv") {
          log_msg(study_name, paste("Leyendo csv: ", suptabB))
          
          #Obtener los diferentes id de las muestras de la partA. Sirve para unir con la partA   
          isolate_id_publication.partA <-  
            partA %>% 
            filter(standard_field %in% c("isolate_id_publication", "sample_accession", "secondary_sample_accession")) %>% 
            select(standard_field, value) %>% 
            mutate(
              value = str_trim(value)) %>% 
            filter(!value %in% c("missing", "")) %>%
            rename(standard_subfield_atb = standard_field)
          
          #Obtener las columnas donde está la información en la tabla suplementaria de la partB
          suptabB.fields2extract <-  
            partB_tidy_atb.suptab %>% 
            filter(source_field == suptabB) %>% 
            select(standard_subfield_atb, value) %>%
            mutate(value = str_replace_all(value, " ", "\\.")) 
          
          log_msg(study_name, "suptabB.fields2extract")
          
          #unir la info de la partB con en isolate_id_publication de la partA
          suptabB.fields2extract_id <- 
            bind_rows(suptabB.fields2extract, isolate_id_publication.partA) %>%
            mutate(value = make.names(value)) #quita espacios y caracteres especailes (guiones)
          
          #info + isolate_id_publications value
          col2extract_partB <- suptabB.fields2extract_id %>% 
            pull(value)
          
          #Leer los campos de interés en la tabla suplementaria
          tryCatch({
            suptabB_info <- fread(suptabB, check.names = T) %>%
              select(all_of(col2extract_partB))
          }, error = function(e) {
            log_msg(study_name, paste("💥 ERROR leyendo tabla suplementaria:", suptabB, "-", conditionMessage(e)))
            next
          })
          
          #suptabB_info <-
          #  read_delim(suptabB, delim = NULL) %>% 
          #  select(all_of(col2extract_partB))

         ### 4.4.2 Extraer y leer los campos de interés en excel #### 
          
        } else if (file_extension %in% c("xls", "xlsx", "XLSX")) {
          log_msg(study_name, paste("Leyendo excel: ", suptabB))
          
          #Buscar las tabs de la suptabB
          tab.B <- 
            partB_tidy %>%
            filter(source_field == suptabB) %>% 
            select(tab_number) %>%
            distinct() %>%
            pull(tab_number)
          
          # Buscar el isolate_id_publication que se encuentra en la tab de los AST
          isolate_id_publication.partB <-
            partA %>% 
            filter(source_field == suptabB & tab_number == tab.B &
                     str_starts(standard_field, "isolate_id_publication")) %>%
            select(standard_field, value) %>% 
            rename(standard_subfield_atb = standard_field)
          
          isolate_id_publication.partB.value <-
            isolate_id_publication.partB %>%
            pull(value)
          
          # Correspondencia de los campos a extraer
          suptabB.fields2extract <-  
            partB_tidy_atb.suptab %>% 
            filter(source_field == suptabB) %>% 
            select(standard_subfield_atb, value)
          
          log_msg(study_name, "suptabB.fields2extract")
          
          # Correspondencia de los campos a extraer con el isolate_id_publication
          suptabB.fields2extract_id.temp <- bind_rows(suptabB.fields2extract, isolate_id_publication.partB)
          
          # Nombres de las columnas a extraer
          col2extract_partB <- suptabB.fields2extract_id.temp %>% 
            pull(value) 
          
          # Obtener isolate_id_publication.partA para después poder unirlo con la partB
          isolate_id_publication.partA <- 
            partA %>% 
            filter(standard_field == "isolate_id_publication") %>%
            select(standard_field, value)%>%
            rename(standard_subfield_atb = standard_field)
          
          suptabB.fields2extract_id <- bind_rows(suptabB.fields2extract, isolate_id_publication.partA)
          
          # Leer los campos de interés de la tabla suplementaria
          tryCatch({
            suptabB_info <- 
              as.data.frame(read_excel(supp_table, sheet = as.integer(tab.B), col_types = "text", progress = FALSE)) %>%
              select(all_of(col2extract_partB)) %>% 
              rename_with(~ isolate_id_publication.1, !!isolate_id_publication.partB.value) 
            # renombramos el isolate_id_publication.partB.value para poder unir todas las partB con el df vacío
          }, error = function(e) {
            log_msg(study_name, paste("💥 ERROR leyendo tabla suplementaria:", supp_table, "-", conditionMessage(e)))
            next
          })
          
        } else {
          
          log_msg(study_name, "Formato de archivo no soportado. Solo se admiten .csv, .xls y .xlsx")
        }  
        
        ## 4.5 Cambiar los nombres de las variables de la tabla suplementaria por los normalizados ####
        
        # Crear el vector de correspondencia nuevamente
        new_colnames <- setNames(suptabB.fields2extract_id$standard_subfield_atb, 
                                 suptabB.fields2extract_id$value)
        
        # Renombrar las columnas de suptab2 usando la correspondencia
        colnames(suptabB_info) <- new_colnames[colnames(suptabB_info)]
        
        log_msg(study_name, "tabla normalizada")
        #print(head(suptabB_info))
        
        #suptabB_info %>%  colnames
        suptabB_info <- 
          suptabB_info %>%
          filter(
            if_any(starts_with("isolate_id_publication"), ~ !is.na(.))
          )
        
        ## 4.6 Unir valores del manuscript al summary_df_partB ####
        
        # Obtener los values y standard_subfield_atb de los campos que vienen del manuscript
        manuscript_fields.partB <-  
          partB_tidy_atb.suptab %>% 
          filter(source_type == "manuscript") %>% 
          select(standard_subfield_atb, value)
        
        log_msg(study_name, "manuscript_fields.partB")
        
        # Pivotar para obtener el df con el standard_subfield_atb y los valores 
        # correspondinetes y poner los valores para todas las filas del summary_df_partA_suptabB 
        expanded_manuscript_fields.partB <- 
          manuscript_fields.partB %>%
          pivot_wider(names_from = standard_subfield_atb, values_from = value) %>% 
          slice(rep(1, nrow(suptabB_info)))

        # Unir info extraída de la tabla suplementaria (suptabB_info) y info del manuscript 
        summary_df_partB <- bind_cols(suptabB_info , expanded_manuscript_fields.partB)
        
        # Ver el resultado
        #print(head(summary_df_partB))
        
        # Unir todas las partsB
        all_summary_df_partB <- append(all_summary_df_partB, list(summary_df_partB))
        
        log_msg(study_name, "all_summary_df_partB")
        
    } #for (partB_file in partB_files)
        
        final_summary_df_partB <- bind_rows(all_summary_df_partB)
        
        final_summary_df_partB %>%  colnames()
        
        #view(final_summary_df_partB)
      
    
    # 5. Unir info de partA y partB ####
    
    # Realizar la unión con merge
    
      summary_df_partA_partB <-
        merge(final_summary_df_partA, final_summary_df_partB,
               by = intersect(names(final_summary_df_partA), names(final_summary_df_partB)), 
               all = TRUE)
    
      summary_df_partA_partB_2 <- 
        summary_df_partA_partB %>% 
        filter(!is.na(isolate_id_publication))%>%
        mutate(across(where(is.list), ~ sapply(., toString)))
      
      muestras <- length(summary_df_partA_partB_2$isolate_id_publication)
      
      sapply(summary_df_partA_partB_2, class)
      
      #print(head(summary_df_partA_partB_2))
    
    # 6. Unir todos los estudios #### 
      all_summaries <-
      merge(all_summaries, summary_df_partA_partB_2,
             by = intersect(names(all_summaries), names(summary_df_partA_partB_2)),
             all = TRUE)
      
      all_summaries_id <- all_summaries %>%
        mutate(internal_id = paste0("ID_", sample(1e6:9e6, n(), replace = FALSE))) %>%
        select(internal_id, everything())
      
      
      #log_msg(study_name, "all_summarise\n")
      log_msg(study_name, paste("Se añadieron", muestras, "muestras\n"))
    
      #view(all_summaries)
      
      dim(all_summaries_id)
    
      #sapply(all_summaries, class)
    
      } #if (length(partA_files) == 0 || length(partB_files) == 0)    
    
  }, error = function(e) {
  log_msg(study_name, paste("💥 ERROR procesando el estudio: ",  conditionMessage(e)))
  }, warning = function(w) {
    log_msg(study_name, paste("⚠️ ADVERTENCIA durante el procesamiento: ", conditionMessage(w)))
  })

  
}#for study in studies

numero_estudios <- length(studies)
total_muestras <- length(all_summaries_id$isolate_id_publication)

log_msg("Resumen", paste("Se analizaron", numero_estudios, "estudios"))
log_msg("Resumen", paste("Se han obtenido un total de", total_muestras, "muestras"))

dim(all_summaries_id)
#warnings()

#view(all_summaries)

# 7. Extraer la información del ENA a través de los bioprojects ####
#Extract ENA data of the bioprojects included in the summaries

#obtener las bioproject accessions de all_summaries
bioproject_accessions <-
  all_summaries_id %>% 
  select(bioproject_accession) %>% 
  separate_rows(bioproject_accession, sep = ";") %>%
  distinct() %>% 
  mutate(across(everything(), str_trim)) %>%   #quitar los espacios de delante y detras 
  filter(!is.na(bioproject_accession)) %>%
  pull()

numero_bioprojects <- length(bioproject_accessions)
log_msg("Resumen", paste("Se obtuvieron muestras procedentes de", numero_bioprojects, "BioProjects"))

# crear un df vacío para guardar todos los datos 
all_data <- tibble()

# bucle para obtener la información del ENA de cada bioproject
for (accession in bioproject_accessions) {
  print(accession)
  
  url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", accession,
    "&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,",
    "experiment_accession,run_accession,sample_alias,sample_title,scientific_name,tax_id,",
    "experiment_alias,run_alias,instrument_platform,instrument_model,library_layout,",
    "collection_date,host,tissue_type,isolation_source,strain,country,collected_by",
    "&format=tsv&download=true"
  )
  
  
  # descargar la información de todos los bioprojects 
  response <- GET(url)
  
  # sacar el contenido de response para obtener el df all_data con toda la info de los bioprojects 
  # se elimina la info de la descarga GET(url) que no nos interesa 
  if (status_code(response) == 200) {
    data <- read_tsv(content(response, "text"), col_types = cols(.default = col_character())) %>% 
      rename_with(~ paste0("ENA_", .)) #todas las columnas obtenidas tendran el prefijo ENA_
    
    data <- mutate(data, ENA_bioproject_accession = accession)  # Add bioproject accession column
    
    all_data <- bind_rows(all_data, data) 
    # Append to dataframe
  } else {
    warning(paste("Failed to download data for", accession))
  }
}

# Print final dataframe
print(all_data)

# crear una función para que las columnas columns_to_check de all_summaries 
# busquen un match en all_data para poder cruzar los datos  
find_corresponding_columns_with_counts <- function(df1, df2, columns_to_check) {
  # Initialize a list to store results
  corresponding_columns <- list()
  
  # Loop through specified columns in df1
  for (col1 in columns_to_check) {
    if (col1 %in% names(df1)) {
      for (col2 in names(df2)) {
        # Check for overlapping values between df1[col1] and df2[col2]
        overlaps <- intersect(na.omit(df1[[col1]]), na.omit(df2[[col2]]))
        
        if (length(overlaps) > 0) {
          # Count how many rows in df1[col1] have matches in df2[col2]
          matching_rows <- sum(na.omit(df1[[col1]]) %in% na.omit(df2[[col2]]))
          # Store the correspondence and the count
          corresponding_columns[[col1]] <- append(
            corresponding_columns[[col1]],
            list(list(
              column_in_df2 = col2,
              matching_rows = matching_rows
            ))
          )
        }
      }
    } else {
      cat(sprintf("Column '%s' not found in df1.\n", col1))
    }
  }
  
  # Return the results as a list or a message if no matches are found
  if (length(corresponding_columns) == 0) {
    return("No corresponding columns found between df1 and df2.")
  } else {
    return(corresponding_columns)
  }
}

# llamar la función para las columnas que son identificadores de muestras  
columns_to_check <- c("isolate_id_publication", "isolate_name_alias", "sample_accession",
                      "secondary_sample_accession", "sample_alias",                                                              
                      "run_accession")

results <- find_corresponding_columns_with_counts(all_summaries_id, all_data, columns_to_check)

# Print els resultados para que queden dispuestos de forma más intuitiva 
if (is.character(results)) {
  cat(results, "\n")
} else {
  for (col in names(results)) {
    cat(sprintf("Column '%s' in df1:\n", col))
    for (info in results[[col]]) {
      cat(sprintf(
        "  Corresponds to column '%s' in df2 with %d matching rows.\n",
        info$column_in_df2, info$matching_rows
      ))
    }
    cat("\n")
  }
}

# Column 'isolate_id_publication' in df1:
#   Corresponds to column 'ENA_run_accession' in df2 with 5248 matching rows.
# Corresponds to column 'ENA_sample_accession' in df2 with 69 matching rows.
# Corresponds to column 'ENA_secondary_sample_accession' in df2 with 254 matching rows.
# Corresponds to column 'ENA_sample_alias' in df2 with 5917 matching rows.
# Corresponds to column 'ENA_sample_title' in df2 with 1894 matching rows.
# Corresponds to column 'ENA_tax_id' in df2 with 16 matching rows.
# Corresponds to column 'ENA_experiment_alias' in df2 with 4654 matching rows.
# Corresponds to column 'ENA_run_alias' in df2 with 14 matching rows.
# Corresponds to column 'ENA_isolation_source' in df2 with 2 matching rows.
# Corresponds to column 'ENA_strain' in df2 with 3325 matching rows.
# 
# Column 'sample_accession' in df1:
#   Corresponds to column 'ENA_sample_accession' in df2 with 1539 matching rows.
# Corresponds to column 'ENA_sample_alias' in df2 with 49 matching rows.
# 
# Column 'secondary_sample_accession' in df1:
#   Corresponds to column 'ENA_secondary_sample_accession' in df2 with 3723 matching rows.
# 
# Column 'run_accession' in df1:
#   Corresponds to column 'ENA_run_accession' in df2 with 8734 matching rows.
# Corresponds to column 'ENA_secondary_sample_accession' in df2 with 487 matching rows.

## cruzar los df 

cross_names <- c("sample_accession", "secondary_sample_accession", "run_accession", "sample_alias")

 # Aplicar cada regex por separado
all_summaries_accessions <- all_summaries_id %>%
  mutate(
    match_sam = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "SAM(E|D|N)[A-Z]?[0-9]+")
    }),
    match_ers = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RS[0-9]{6,}")
    }),
    match_erx = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RX[0-9]{6,}")
    }),
    match_err = apply(., 1, function(row) {
      str_extract(paste(row, collapse = " "), "(E|D|S)RR[0-9]{6,}")
    })
  )

dim(all_summaries_id) 
dim(all_summaries_accessions)

  # buscar las entradas con el formato del ENA y cruzarlas

  # ver cuántas entradas hay con cada formato 

dim_ERS <-  all_summaries_accessions %>%
  filter(!is.na(match_ers))%>%
  dim()
  
dim_ERX <-  all_summaries_accessions %>%
  filter(!is.na(match_erx))%>%
  dim()

dim_ERR <-  all_summaries_accessions %>%
  filter(!is.na(match_err)) %>%
  distinct()%>%
  dim()

dim_SAM <- all_summaries_accessions %>%
  filter(!is.na(match_sam)) %>%
  distinct()%>%
  dim()

# empezamos cruzando las match_err con run_accession pq tiene mayor numero de muestras 

data_ERR <- all_summaries_accessions %>%
  filter(!is.na(match_err)) %>% 
  distinct()%>%
  inner_join(all_data, by = c("match_err" = "ENA_run_accession")) %>%
  rename(ENA_run_accession = match_err) 

dim(data_ERR)

# obtener los match_err cruzados 
err <- data_ERR %>% pull(ENA_run_accession)

# cruzar match_ers con ENA_secondary_sample_accession (sin tener en cuenta las match_err)
data_ERS <-     
  all_summaries_accessions %>%
  filter(!is.na(match_ers)) %>% 
  distinct()%>%
  filter(!match_err %in% err)%>%
  inner_join(all_data, by = c("match_ers" = "ENA_secondary_sample_accession")) %>%
  rename(ENA_secondary_sample_accession = match_ers)

dim(data_ERS)

# obtener los match_ers cruzados 
ers <- data_ERS %>% pull(ENA_secondary_sample_accession)

# cruzar match_sam con ENA_experiment_accession (sin tener en cuenta las muestras ya cruzadas)
data_SAM <-     
  all_summaries_accessions %>%
  filter(!is.na(match_sam)) %>% 
  distinct()%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  inner_join(all_data, by = c("match_sam" = "ENA_sample_accession")) %>%
  rename(ENA_sample_accession = match_sam)

dim(data_SAM)

# obtener los match_sam cruzados 
sam <- data_SAM %>% pull(ENA_sample_accession)

# cruzar match_erx con ENA_experiment_accession (sin tener en cuenta las muestras ya cruzadas)
data_ERX <-     
  all_summaries_accessions %>%
  filter(!is.na(match_erx)) %>% 
  distinct()%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  inner_join(all_data, by = c("match_erx" = "ENA_experiment_accession"))%>%
  rename(ENA_experiment_accession = match_erx)

dim(data_ERX)

# obtener los match_erx cruzados 
erx <- data_ERX %>% pull(ENA_experiment_accession)

#crurzar isolate_id con sample_alias (sin tener en cuenta las muestras ya cruzadas)

# Paso 1: hacer un join solo por sample alias
data_sample_alias.temp <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  filter(!match_sam %in% sam)%>%
  filter(!match_erx %in% erx)%>%
  inner_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_sample_alias")
  )


# Paso 2: filtrar por si el ENA_bioproject_accession es el mismo que el bioproject_accession del estudio 
# así aseguramos que no se han cruzado muestras diferentes con el mismo nombre 
data_sample_alias <- data_sample_alias.temp %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

dim(data_sample_alias)

# obtener los sample_alias cruzados 
names_sample_alias <- data_sample_alias %>% pull(isolate_id_publication)

#crurzar isolate_id con sample_title (sin tener en cuenta las muestras ya cruzadas)

# Paso 1: hacer un join solo por sample title 
data_sample_title.temp <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  filter(!match_sam %in% sam)%>%
  filter(!match_erx %in% erx)%>%
  filter(!isolate_id_publication %in% names_sample_alias)%>%
  inner_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_sample_title")
  )


# Paso 2: filtrar por si el ENA_bioproject_accession es el mismo que el bioproject_accession del estudi 
# así aseguramos que no se han cruzado muestras diferentes con el mismo nombre 
data_sample_title <- data_sample_title.temp %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

# obtener los sample_title cruzados 
names_sample_title <- data_sample_title %>% pull(isolate_id_publication)
dim(data_sample_title)

#comprobar pq no sale ningun sample_title 
data_sample_title.temp2 <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  inner_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_sample_title")
  )

data_sample_title.2 <- data_sample_title.temp2 %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

colnames(data_sample_title.2)


#crurzar isolate_id con experiment_alias (sin tener en cuenta las muestras ya cruzadas)

# Paso 1: hacer un join solo por experiment alias
data_experiment_alias.temp <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  filter(!match_sam %in% sam)%>%
  filter(!match_erx %in% erx)%>%
  filter(!isolate_id_publication %in% names_sample_alias)%>%
  filter(!isolate_id_publication %in% names_sample_title)%>%
  inner_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_experiment_alias")
  )

# Paso 2: filtrar por si el ENA_bioproject_accession es el mismo que el bioproject_accession del estudi 
# así aseguramos que no se han cruzado muestras diferentes con el mismo nombre 
data_experiment_alias <- data_experiment_alias.temp %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

dim(data_experiment_alias)

# obtener los experiment_alias cruzados 
names_experiment_alias <- data_experiment_alias %>% pull(isolate_id_publication)


#___________________________________________________

# comprobar que no hay experiment_alias cruzadas pq esas muestras ya se han cruzado
data_experiment_alias.temp2 <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  left_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_experiment_alias")
  )

data_experiment_alias.2 <- data_experiment_alias.temp2 %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

colnames(data_sample_title.2)

#___________________________________________________

#crurzar isolate_id con ENA_strain (sin tener en cuenta las muestras ya cruzadas)

# Paso 1: hacer un join solo por strain
data_strain.temp <- all_summaries_accessions %>%
  filter(!is.na(isolate_id_publication))%>%
  filter(!match_err %in% err)%>%
  filter(!match_ers %in% ers)%>%
  filter(!match_sam %in% sam)%>%
  filter(!match_erx %in% erx)%>%
  filter(!isolate_id_publication %in% names_sample_alias)%>%
  filter(!isolate_id_publication %in% names_sample_title)%>%
  filter(!isolate_id_publication %in% names_experiment_alias)%>%
  inner_join(
    all_data,
    by = c("isolate_id_publication" = "ENA_strain")
  )

# Paso 2: filtrar por si el ENA_bioproject_accession es el mismo que el bioproject_accession del estudi 
# así aseguramos que no se han cruzado muestras diferentes con el mismo nombre 
data_strain <- data_strain.temp %>%
  filter(
    pmap_lgl(
      list(bioproject_accession, ENA_bioproject_accession), 
      ~ str_detect(.x, fixed(.y))
    )
  )

# obtener los strain cruzados 
names_strain <- data_experiment_alias %>% pull(isolate_id_publication)
dim(data_strain)

#unir todas las muestras que se han podido cruzar en df_cruzado
df_cruzado <- 
  bind_rows(data_ERR, data_ERS, data_SAM, data_ERX, data_sample_alias, data_sample_title, data_experiment_alias, data_strain) 

dim(df_cruzado)

colnames(df_cruzado)

#view(df_cruzado)


#------------------------------------------------------------------#

#buscar las 4653 muestras que no se han cruzado con el ENA 

no_ENA <- anti_join(all_summaries_accessions, df_cruzado, by = "internal_id")

muestras_no_cruzadas <- length(no_ENA$isolate_id_publication)

no_ENA %>% pull(internal_study_name)%>% unique()

log_msg("Resumen", paste("No se pudieron cruzar", muestras_no_cruzadas, "muestras"))

#------------------------------------------------------------------#


# 8. Devolver tabla final ####
write.csv(df_cruzado, paste0(base_path, "/df_cruzado_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

#log_msg("====== FIN DEL SCRIPT ======")
#sink(type = "message")
#sink()

#rm(list=ls())

#sink(NULL)
#close(log_connection)

#colnames(all_summaries_id)

# ====== FIN DEL SCRIPT ======
sink(type = "message")
sink()
sink(NULL)

# 1. Primero cerramos el archivo log correctamente
close(log_connection)

# 2. Después puedes mirar tus columnas (porque la tabla aún existe)
colnames(all_summaries_id)

# 3. Y SOLO AL FINAL DEL TODO limpiamos la memoria (si quieres)
rm(list=ls())


