#!/bin/bash
#$ -N fastqc_raw
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -q default


# Definir las rutas (Ajusta según tu entorno)
# RAW_DIR debe apuntar a la carpeta "madre" donde tienes separados los modelos
RAW_DIR="/export/space3/users/pedropm/proyecto_transcriptomica/rawdata" 
OUT_DIR="/export/space3/users/pedropm/proyecto_transcriptomica/fastqc_raw_results"
ADAPTERS_FILE="/export/space3/users/pedropm/proyecto_transcriptomica/adapters.txt"

# Crear el directorio de salida si no existe
mkdir -p "$OUT_DIR"

echo "Iniciando análisis FastQC de datos crudos..."

# Buscar recursivamente todos los .fastq o .fastq.gz y ejecutar FastQC
find "$RAW_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fastq" \) | while read FILE; do
    echo "Procesando: $(basename "$FILE")"
    
    # Tu comando optimizado con la bandera de adaptadores
    fastqc -t 8 -a "$ADAPTERS_FILE" -o "$OUT_DIR" "$FILE"
done

echo "FastQC completado para todas las muestras crudas."


