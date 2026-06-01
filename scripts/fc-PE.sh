#!/bin/bash
#$ -N featurecounts
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 30
#$ -q default

# Aborta el script inmediatamente si cualquier comando falla.
set -e

echo "------------------------------------"
echo "Job-ID: $JOB_ID"
echo "Hostname: $(hostname)"
echo "Fecha: $(date)"
echo "Directorio de trabajo: $(pwd)"
echo "------------------------------------"

# 1) Inicializa Mamba/Conda(para deepTools)
echo "Activando el entorno de Conda..."
source /etc/bashrc
conda activate subread
echo "Entorno de Conda activado con éxito."

# 2) Verifica que los comandos clave existan
echo "Verificando la ubicación de los comandos..."
which conda
which parallel
which featureCounts
echo "Verificación de comandos completa."


if [ -z "$1" ];
then
    echo "ERROR: No se proporcionó la ruta del directorio de entrada (Argumento $1)."
    exit 1
fi
INPUT_DIR=$1
echo "Directorio de entrada: $INPUT_DIR"
echo "------------------------------------"

if [ -z "$2" ];
then
    echo "ERROR: No se proporcionó la ruta del directorio de salida (Argumento $2)."
    exit 1
fi
OUTPUT_FILE=$2
echo "Archivo de salida: $OUTPUT_FILE"
echo "------------------------------------"

if [ -z "$3" ];
then
    echo "ERROR: No se proporcionó la lista de muestras (Argumento $3)."
    exit 1
fi
GTF_FILE=$3
echo "Archivo GTF utilizado: $GTF_FILE"
echo "------------------------------------"


featureCounts -T 30 \
              -Q 10 \
              -p -B \
              -t exon \
              -g gene_id \
              -a "$GTF_FILE" \
              -o "$OUTPUT_FILE" \
              "$INPUT_DIR"/*_Aligned.out.bam
echo "------------------------------------"
echo "Cuantificación PE finalizada. Matriz generada en: $OUTPUT_FILE"