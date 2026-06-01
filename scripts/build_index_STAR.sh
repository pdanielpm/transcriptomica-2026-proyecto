#!/bin/bash
#$ -N STAR_Genome_Index
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 15
#$ -q default
#$ -j y

set -e

echo "Iniciando generación de índices para STAR"

# Carga de módulos (Ajusta a los equivalentes de tu servidor o actívalos vía conda)
source /etc/bashrc
conda activate star


# 1. Definir rutas absolutas de tus archivos de referencia
FASTA_FILE="/export/space3/users/pedropm/proyecto_transcriptomica/Danio_rerio.GRCz11.dna.primary_assembly.fa"
GTF_FILE="/export/space3/users/pedropm/proyecto_transcriptomica/Danio_rerio.GRCz11.115.chr.gtf"
BASE_DIR="$PWD/STAR_Indices"

mkdir -p "$BASE_DIR"

# 2. Bucle para generar los índices de 50, 75 y 150 pb
# La sintaxis es "TAMAÑO_READ:OVERHANG"
for CONFIG in "50:49" "75:74" "150:149"; do
    
    # Extraer el nombre de la carpeta y el valor de overhang
    SIZE="${CONFIG%%:*}"
    OVERHANG="${CONFIG##*:}"
    
    OUT_DIR="${BASE_DIR}/${SIZE}genomeDirectory"
    mkdir -p "$OUT_DIR"
    
    echo "--------------------------------------------------------"
    echo "Generando índice para lecturas de ${SIZE}pb (Overhang: ${OVERHANG})"
    
    STAR --genomeDir "$OUT_DIR" \
         --runThreadN 15 \
         --runMode genomeGenerate \
         --genomeFastaFiles "$FASTA_FILE" \
         --sjdbGTFfile "$GTF_FILE" \
         --genomeSAindexNbases 10 \
         --sjdbOverhang "$OVERHANG"
         
    echo "Índice de ${SIZE}pb completado."
done
