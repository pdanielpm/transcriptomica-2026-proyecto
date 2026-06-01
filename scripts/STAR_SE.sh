#!/bin/bash
#$ -N STAR_pe
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


# 2) Verifica que los comandos clave existan
echo "Verificando la ubicación de los comandos..."
which conda
which parallel
echo "Verificación de comandos completa."

# 1) Inicializa Mamba/Conda(para STAR)
echo "Activando el entorno de Conda..."
source /etc/bashrc
conda activate star
echo "Entorno de Conda activado con éxito."


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

if [ -z "$4" ];
then
    echo "ERROR: No se proporcionó la ruta del índice de STAR (Argumento $4)."
    exit 1
fi
STAR_INDEX=$4
echo "Índice de STAR: $STAR_INDEX"
echo "------------------------------------" 

export SAMPLE_LIST INPUT_DIR OUTPUT_DIR STAR_INDEX

echo "Iniciando alineamiento SE con STAR..."

parallel -j 3 "
    echo 'Alineando muestra: {}'
    
    STAR --runThreadN 10 \\
         --genomeLoad NoSharedMemory \\
         --runMode alignReads \\
         --quantMode TranscriptomeSAM GeneCounts \\
         --outSAMtype BAM Unsorted \\
         --outFileNamePrefix \"$OUTPUT_DIR/{}_\" \\
         --genomeDir $STAR_INDEX \\
         --readFilesIn $INPUT_DIR/{}_clean.fastq.gz \\
         --readFilesCommand zcat \\
         --outFilterMultimapNmax 15 \\
         --outFilterMismatchNoverLmax 0.06 \\
         --outFilterMatchNmin 16
" < "$SAMPLE_LIST"

echo "Alineamiento SE completado."