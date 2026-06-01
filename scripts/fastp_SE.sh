#!/bin/bash
#$ -N fastp_SE
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 25
#$ -q default

# Aborta el script inmediatamente si cualquier comando falla.
set -e

echo "------------------------------------"
echo "Job-ID: $JOB_ID"
echo "Hostname: $(hostname)"
echo "Fecha: $(date)"
echo "Directorio de trabajo: $(pwd)"
echo "------------------------------------"

# 1) Inicializa Mamba/Conda(para fastp)
echo "Activando el entorno de Conda..."
source /etc/bashrc
conda activate fastp
echo "Entorno de Conda activado con éxito."

# 2) Verifica que los comandos clave existan
echo "Verificando la ubicación de los comandos..."
which conda
which parallel
which fastp
echo "Verificación de comandos completa."

if [ -z "$1" ]; then
    echo "ERROR: No se proporcino la lista de muestras (Argumento $1)."
    exit 1
fi
SAMPLE_LIST="$1"
echo "Ruta del archivo de muestras: $SAMPLE_LIST"
echo "------------------------------------"

if [ -z "$2" ];
then
    echo "ERROR: No se proporcionó la ruta del directorio de entrada (Argumento $2)."
    exit 1
fi
INPUT_DIR=$2
echo "Directorio de entrada: $INPUT_DIR"
echo "------------------------------------"

if [ -z "$3" ];
then
    echo "ERROR: No se proporcionó la ruta del directorio de salida (Argumento $3)."
    exit 1
fi
OUTPUT_DIR=$3
echo "Directorio de salida: $OUTPUT_DIR"
echo "------------------------------------"

ADAPTER_FASTA="/export/space3/users/pedropm/proyecto_transcriptomica/adpaters.fa"
echo "Archivo de adaptadores: $ADAPTER_FASTA"
echo "------------------------------------"
export INPUT_DIR OUTPUT_DIR ADAPTER_FASTA

# 3) Ejecuta el proceso principal

echo "Iniciando limpieza con fastp en paralelo..."

parallel -j 5 "
    echo 'Procesando muestra: {}'
    fastp -i $INPUT_DIR/{}.fastq.gz \
          -o $OUTPUT_DIR/{}_clean.fastq.gz \
          --adapter_fasta $ADAPTER_FASTA \
          -3 -5 -M 20 -q 20 -l 25 -g -x -p -P 15 \
          --thread 5 \
          -h $OUTPUT_DIR/{}_fastp.html \
          -j $OUTPUT_DIR/{}_fastp.json
" < "$SAMPLE_LIST"

echo "Limpieza SE completada con éxito."