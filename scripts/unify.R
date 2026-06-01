unificar_featurecounts <- function(directorio_entrada, archivo_salida) {
  cat("Buscando archivos en:", directorio_entrada, "...\n")

  # 1. Buscamos el patrón correcto de los archivos que generaste
  archivos <- list.files(
    directorio_entrada,
    pattern = "_counts\\.txt$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(archivos) == 0) {
    stop("No encontré archivos terminados en _counts.txt en esa ruta.")
  }

  # 2. Leer el primer archivo (Tomamos Geneid y TODAS las columnas de conteos)
  # featureCounts guarda los conteos a partir de la columna 7
  tabla_maestra <- read.delim(
    archivos[1],
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  tabla_maestra <- tabla_maestra[, c(1, 7:ncol(tabla_maestra))]

  # 3. Bucle para pegar los demás grupos
  if (length(archivos) > 1) {
    for (i in 2:length(archivos)) {
      temp <- read.delim(
        archivos[i],
        comment.char = "#",
        stringsAsFactors = FALSE
      )
      temp <- temp[, c(1, 7:ncol(temp))]

      # Merge seguro usando el Geneid como ancla
      tabla_maestra <- merge(tabla_maestra, temp, by = "Geneid", all = TRUE)
    }
  }

  # 4. Limpieza dinámica de nombres de columnas
  # Extrae estrictamente la cadena que empieza con SRR seguida de números,
  # ignorando toda la ruta de carpetas que R haya leído del archivo.
  colnames(tabla_maestra) <- gsub(
    ".*(SRR[0-9]+).*",
    "\\1",
    colnames(tabla_maestra)
  )

  # 5. Acomodar rownames
  rownames(tabla_maestra) <- tabla_maestra$Geneid
  tabla_maestra$Geneid <- NULL

  write.table(
    tabla_maestra,
    archivo_salida,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
  )
  cat(
    "✔ Matriz unificada guardada en:",
    archivo_salida,
    "(",
    ncol(tabla_maestra),
    "muestras)\n"
  )
}

# ====================================================================
# EJECUCIÓN
# ====================================================================

unificar_featurecounts(
  directorio_entrada = "/export/space4/users/pedropm/proyecto_transcriptomica/featurecounts",
  archivo_salida = "/export/space4/users/pedropm/proyecto_transcriptomica/featurecounts/unified_regeneration.txt"
)
