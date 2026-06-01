library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)
library(limma)
library(sva)

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    BiocManager::install("EnhancedVolcano", ask = FALSE)
}
library(EnhancedVolcano)

filtrar_degs <- function(res_obj, nombre_contraste, lfc_cutoff, fdr_cutoff ) {
  # Convertir a data frame y quitar NAs
  res_df <- as.data.frame(res_obj)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  # UP-regulated: FDR < 0.05 y LFC > 1
  up <- res_df[res_df$padj < fdr_cutoff & res_df$log2FoldChange > lfc_cutoff, ]
  # DOWN-regulated: FDR < 0.05 y LFC < -1
  down <- res_df[res_df$padj < fdr_cutoff & res_df$log2FoldChange < -lfc_cutoff, ]

  cat(nombre_contraste, "DESeq2 — UP:", nrow(up), "| DOWN:", nrow(down), "\n")
  
  # Guardar en CSV
  write.csv(up, paste0("UP_", nombre_contraste, ".csv"))
  write.csv(down, paste0("DOWN_", nombre_contraste, ".csv"))
  
  cat(nombre_contraste, "- UP:", nrow(up), "| DOWN:", nrow(down), "\n")
  return(list(UP = up, DOWN = down))
}

plot_volcano <- function(res_obj, titulo) {

  etiquetas_personalizadas <- c('NS',
                                'Solo Log2 FC',
                                'Solo Significancia (FDR)',
                                'Significativo (FDR y Log2 FC)')
  
  contrastvolcano <- EnhancedVolcano(res_obj,
    lab = NA,                     
    x = 'log2FoldChange',
    y = 'padj',                  
    title = titulo,
    subtitle = "Umbrales: FDR < 0.05 | |Log2FC| > 1",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 1.2,              # Puntos más finos para reducir el caos visual
    labSize = 0,                  # Tamaño de etiqueta en 0 por seguridad
    col = c("grey70", "forestgreen", "royalblue", "firebrick3"), 
    colAlpha = 0.5,              
    legendPosition = 'bottom',    
    legendLabels = etiquetas_personalizadas,
    legendLabSize = 11,
    legendIconSize = 4.0,
    drawConnectors = FALSE,       
    gridlines.major = TRUE,       
    gridlines.minor = FALSE,      
    border = 'partial',
    borderWidth = 1.0,
    borderColour = 'black'
  )
  
  return(contrastvolcano)
}

# ====================================================================
# 1. METADATOS: CONSTRUCCIÓN Y CORRECCIÓN
# ====================================================================
metadata <- data.frame(
  row.names = c(
    # -- Ablation (SE) --
    "SRR11294120", "SRR11294121", "SRR11294122", "SRR3143894", "SRR3143897",
    # -- Amputated (PE) --
    "SRR12554274", "SRR12554275", "SRR12554276", "SRR8867570", "SRR8867571", "SRR8867572",
    # -- Cryoinjury (SE) --
    "SRR5809472", "SRR5809473", "SRR5809474", "SRR5809475", "SRR6910771", "SRR6910772", "SRR6910773",
    # -- Sham (PE) --
    "SRR12554271", "SRR12554272", "SRR12554273",
    # -- Sham (SE) --
    "SRR5809464", "SRR5809465", "SRR5809466", "SRR5809467", "SRR6910774", "SRR6910775", "SRR6910776",
    # -- Uninjured PE (LAS DEL ERROR DE GEO) --
    "SRR11033311", "SRR11033312", "SRR11033313",
    # -- Uninjured SE --
    "SRR11294126", "SRR11294127", "SRR11294128", "SRR3143899", "SRR3143916"
  ),
  Condition_Raw = c(
    rep("Ablation", 5),
    rep("Amputation", 6),
    rep("Cryoinjury", 7),
    rep("Sham", 3),
    rep("Sham", 7),
    rep("Uninjured", 3), 
    rep("Uninjured", 5)
  ),
  Batch = c(
    rep("GSE146859", 3), rep("GSE75894", 2),
    rep("GSE157170", 3), rep("GSE129499", 3),
    rep("GSE100892", 4), rep("GSE112452", 3),
    rep("GSE157170", 3),
    rep("GSE100892", 4), rep("GSE112452", 3),
    rep("GSE144831", 3),
    rep("GSE146859", 3), rep("GSE75894", 2)
  ),
  Platform = c(
    rep("BGI-SEQ-500", 3), rep("Genome Analyzer II", 2),
    rep("NovaSeq 6000", 3), rep("HiSeq X Ten", 3),
    rep("NextSeq 500", 4), rep("BGI-SEQ-500", 3),
    rep("NovaSeq 6000", 3),
    rep("NextSeq 500", 4), rep("BGI-SEQ-500", 3),
    rep("HiSeq 2000", 3),
    rep("BGI-SEQ-500", 3), rep("Genome Analyzer II", 2)
  ),
  Layout = c(
    rep("Single", 3), rep("Single", 2),
    rep("Paired", 3), rep("Paired", 3),
    rep("Single", 4), rep("Single", 3),
    rep("Paired", 3),
    rep("Single", 4), rep("Single", 3),
    rep("Paired", 3),
    rep("Single", 3), rep("Single", 2)
  ),
  stringsAsFactors = FALSE
)


# Mislabeling GSE144831
metadata$Condition <- metadata$Condition_Raw
metadata$Condition[metadata$Batch == "GSE144831" & metadata$Condition_Raw == "Uninjured"] <- "Amputation"

# Factores y Baseline (Sham)
metadata$Condition <- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "Sham")
metadata$Layout <- factor(metadata$Layout)



# ====================================================================
# 2. FILTRADO Y CORRECCIÓN DE LOTE CON COMBAT-SEQ 
# ====================================================================

counts_data <- read.delim("/export/space4/users/pedropm/proyecto_transcriptomica/featurecounts/unified_regeneration.txt", 
                          row.names = 1, stringsAsFactors = FALSE)
counts_data <- counts_data[, rownames(metadata)]
stopifnot(all(rownames(metadata) == colnames(counts_data)))

cat("\n--- PASO A: Filtrado con edgeR (Igual al paper) ---\n")
# Creamos un objeto temporal de edgeR solo para aprovechar su filtro inteligente
y_temp <- DGEList(counts = counts_data, group = metadata$Condition)

# filterByExpr evalúa el tamaño de librería y pide expresión consistente por grupo
keep_edgeR <- filterByExpr(y_temp, group = metadata$Condition)

# Cortamos la matriz ANTES de ComBat-Seq
counts_filtrados <- counts_data[keep_edgeR, ]
cat("Genes retenidos tras el filtro:", nrow(counts_filtrados), "de", nrow(counts_data), "\n")

cat("\n--- PASO B: Ejecutando ComBat-Seq en la matriz limpia ---\n")
# Ahora sí, a ComBat-Seq solo entra la señal real
matriz_limpia_num <- as.matrix(counts_filtrados)

adjusted_counts <- ComBat_seq(counts = matriz_limpia_num, 
                              batch = metadata$Platform, 
                              group = metadata$Condition)
cat("✔ Matriz ajustada con éxito.\n")

# ====================================================================
# 3. EL MODELO DESEQ2 DEFINITIVO (Usando conteos ajustados)
# ====================================================================
cat("\n--- Generando el modelo DESeq2 ---\n")

dds_combat <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                     colData = metadata,
                                     design = ~ Condition)



# Ajuste del modelo matemático
dds_combat <- DESeq(dds_combat)


# 1. Extraer resultados crudos
res_resection_raw  <- results(dds_combat, contrast = c("Condition", "Amputation", "Sham"))
res_ablation_raw   <- results(dds_combat, contrast = c("Condition", "Ablation", "Sham"))
res_uninjured_raw  <- results(dds_combat, contrast = c("Condition", "Uninjured", "Sham"))
res_cryoinjury_raw <- results(dds_combat, contrast = c("Condition", "Cryoinjury", "Sham"))

# 2. Aplicar el lfcShrink con el método "ashr" (como dicta el paper)
res_resection  <- lfcShrink(dds_combat, contrast = c("Condition", "Amputation", "Sham"), res = res_resection_raw, type = "ashr")
res_ablation   <- lfcShrink(dds_combat, contrast = c("Condition", "Ablation", "Sham"), res = res_ablation_raw, type = "ashr")
res_uninjured  <- lfcShrink(dds_combat, contrast = c("Condition", "Uninjured", "Sham"), res = res_uninjured_raw, type = "ashr")
res_cryoinjury <- lfcShrink(dds_combat, contrast = c("Condition", "Cryoinjury", "Sham"), res = res_cryoinjury_raw, type = "ashr")


dds_raw <- DESeqDataSetFromMatrix(countData = counts_data,
                                  colData = metadata,
                                  design = ~ 1)

vsd_raw <- vst(dds_raw, blind = TRUE)

# 1. Extraer los datos del PCA de la matriz cruda (vsd_raw)
pca_data_raw <- plotPCA(vsd_raw, intgroup = c("Condition", "Batch"), returnData = TRUE)
percentVar_raw <- round(100 * attr(pca_data_raw, "percentVar"))

pca_raw_plot <- ggplot(pca_data_raw, aes(x = PC1, y = PC2, color = Condition, shape = Batch)) +
  geom_point(size = 5, alpha = 0.8) +    
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8, 9)) + 
  xlab(paste0("PC1: ", percentVar_raw[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_raw[2], "% variance")) +
  ggtitle("PCA de Raw Data", 
          subtitle = "Matriz cruda (Previo a corrección de lote con ComBat-Seq)") +
  theme_bw(base_size = 14) +           
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )


print(pca_raw_plot)
ggsave("PCA_CrudoFinal.pdf", plot = pca_raw_plot, width = 8, height = 6)

# Heatmap de los 500 genes más variables (Crudo)
# 1. Calcular la varianza para cada gen a lo largo de todas las muestras
topVarGenes <- head(order(rowVars(assay(vsd_raw)), decreasing = TRUE), 500)

# 2. Extraer la matriz de expresión solo para esos top 500 genes
mat_raw <- assay(vsd_raw)[ topVarGenes, ]

# 3. Centrar los datos gen por gen (Z-score por fila) para que los colores contrasten bien
mat_raw <- mat_raw - rowMeans(mat_raw)

# 4. Preparar las anotaciones de las columnas (barras de colores arriba del heatmap)
df_col <- as.data.frame(colData(vsd_raw)[, c("Condition", "Platform")])

# 5. Generar y guardar el Heatmap directamente en PDF
pdf("Heatmap_Crudo_Top500.pdf", width = 10, height = 8)
pheatmap(mat_raw,
         annotation_col = df_col,
         show_rownames = FALSE,     
         cluster_cols = TRUE,        
         cluster_rows = TRUE,        
         scale = "row",             
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Perfil de Expresión de los 500 Genes Más Variables (Sin corrección)")
dev.off()


vsd_combat <- vst(dds_combat, blind = FALSE)

pcaData_combat <- plotPCA(vsd_combat, intgroup=c("Condition", "Batch"), returnData=TRUE)
percentVar_combat <- round(100 * attr(pcaData_combat, "percentVar"))

pca_final <- ggplot(pcaData_combat, aes(PC1, PC2, color=Condition, shape=Batch)) +
  geom_point(size=5, alpha=0.85) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8)) +
  xlab(paste0("PC1: ", percentVar_combat[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_combat[2], "% variance")) +
  theme_bw() + theme(text = element_text(size=14, face="bold")) +
  ggtitle("PCA Con Corrección de Lote", 
          subtitle = "Matriz ajustada con ComBat-Seq")

print(pca_final)
ggsave("PCA_CombatSeq.pdf", plot = pca_final, width = 8, height = 6, device = "pdf")

lfc_cutoff = 1
fdr_cutoff = 0.05

degs_resection  <- filtrar_degs(res_resection, "Amputation_vs_Sham", lfc_cutoff, fdr_cutoff)
degs_ablation   <- filtrar_degs(res_ablation, "Ablation_vs_Sham", lfc_cutoff, fdr_cutoff)
degs_cryoinjury <- filtrar_degs(res_cryoinjury, "Cryoinjury_vs_Sham", lfc_cutoff, fdr_cutoff)
degs_uninjured  <- filtrar_degs(res_uninjured, "Uninjured_vs_Sham", lfc_cutoff, fdr_cutoff)

volcano_resection  <- plot_volcano(res_resection, "Amputation vs Sham")
volcano_ablation   <- plot_volcano(res_ablation, "Ablation vs Sham")
volcano_cryoinjury <- plot_volcano(res_cryoinjury, "Cryoinjury vs Sham")
volcano_uninjured  <- plot_volcano(res_uninjured, "Uninjured vs Sham")

ggsave("Volcano_Amputation_vs_Sham.pdf", plot = volcano_resection, width = 8, height = 6)
ggsave("Volcano_Ablation_vs_Sham.pdf", plot = volcano_ablation, width = 8, height = 6)
ggsave("Volcano_Cryoinjury_vs_Sham.pdf", plot = volcano_cryoinjury, width = 8, height = 6)
ggsave("Volcano_Uninjured_vs_Sham.pdf", plot = volcano_uninjured, width = 8, height = 6)

# ====================================================================
# EXTRACCIÓN DEL BACKGROUND PARA ORA
# ====================================================================

background_genes <- rownames(dds_combat)

# Guardamos la lista en un archivo de texto simple 
write.table(background_genes, file = "Background_Expressed_Genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(" Background guardado. Total de genes :", length(background_genes), "\n")


# ====================================================================
# EXTRACCIÓN DEL CORE REGENERATIVO Y FIRMAS ESPECÍFICAS
# ====================================================================

cat("\n--- Extrayendo IDs de todas las condiciones ---\n")

# 1. Extraer los IDs (rownames) de las listas UP
up_amp  <- rownames(degs_resection$UP)
up_cryo <- rownames(degs_cryoinjury$UP)
up_abla <- rownames(degs_ablation$UP)
up_unin <- rownames(degs_uninjured$UP)

# 2. Extraer los IDs (rownames) de las listas DOWN
down_amp  <- rownames(degs_resection$DOWN)
down_cryo <- rownames(degs_cryoinjury$DOWN)
down_abla <- rownames(degs_ablation$DOWN)
down_unin <- rownames(degs_uninjured$DOWN)


cat("\n--- Calculando Listas Crudas (Restando Uninjured) ---\n")

# ---------------------------------------------------------
# A. CORE REGENERATIVO (Intersección pura de las 3 lesiones)
# ---------------------------------------------------------
# UP
core_up_raw <- Reduce(intersect, list(up_amp, up_cryo, up_abla))
core_up_puro <- setdiff(core_up_raw, up_unin)

# DOWN
core_down_raw <- Reduce(intersect, list(down_amp, down_cryo, down_abla))
core_down_puro <- setdiff(core_down_raw, down_unin)

cat("Core UP (Crudo -> Puro):", length(core_up_raw), "->", length(core_up_puro), "genes.\n")
cat("Core DOWN (Crudo -> Puro):", length(core_down_raw), "->", length(core_down_puro), "genes.\n")


# ---------------------------------------------------------
# B. FIRMAS ESPECÍFICAS (Exclusivos de cada lesión)
# ---------------------------------------------------------
# Amputación
amp_spec_up_raw <- setdiff(up_amp, union(up_cryo, up_abla))
amp_spec_up_puro <- setdiff(amp_spec_up_raw, up_unin)

amp_spec_down_raw <- setdiff(down_amp, union(down_cryo, down_abla))
amp_spec_down_puro <- setdiff(amp_spec_down_raw, down_unin)

# Criolesión
cryo_spec_up_raw <- setdiff(up_cryo, union(up_amp, up_abla))
cryo_spec_up_puro <- setdiff(cryo_spec_up_raw, up_unin)

cryo_spec_down_raw <- setdiff(down_cryo, union(down_amp, down_abla))
cryo_spec_down_puro <- setdiff(cryo_spec_down_raw, down_unin)

# Ablación
abla_spec_up_raw <- setdiff(up_abla, union(up_amp, up_cryo))
abla_spec_up_puro <- setdiff(abla_spec_up_raw, up_unin)

abla_spec_down_raw <- setdiff(down_abla, union(down_amp, down_cryo))
abla_spec_down_puro <- setdiff(abla_spec_down_raw, down_unin)

cat("Amputación UP (Crudo -> Puro):", length(amp_spec_up_raw), "->", length(amp_spec_up_puro), "genes.\n")
cat("Criolesión UP (Crudo -> Puro):", length(cryo_spec_up_raw), "->", length(cryo_spec_up_puro), "genes.\n")
cat(" Ablación UP (Crudo -> Puro):", length(abla_spec_up_raw), "->", length(abla_spec_up_puro), "genes.\n")


# ====================================================================
# EXPORTAR LAS LISTAS DEFINITIVAS PARA CYTOSCAPE / SHINYGO
# ====================================================================
cat("\n--- Exportando archivos de texto plano (.txt) ---\n")

# Core
write.table(core_up_puro, "Core_Regeneration_UP_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(core_down_puro, "Core_Regeneration_DOWN_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Específicos Amputación
write.table(amp_spec_up_puro, "Amputation_Specific_UP_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(amp_spec_down_puro, "Amputation_Specific_DOWN_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Específicos Criolesión
write.table(cryo_spec_up_puro, "Cryoinjury_Specific_UP_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cryo_spec_down_puro, "Cryoinjury_Specific_DOWN_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Específicos Ablación
write.table(abla_spec_up_puro, "Ablation_Specific_UP_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(abla_spec_down_puro, "Ablation_Specific_DOWN_Definitivo.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("✔ Archivos purificados y exportados con éxito. Listos para ShinyGO/Cytoscape.\n")


# ====================================================================
# VISUALIZACIÓN DE INTERSECCIONES (UpSet Plots)
# ====================================================================
library(ComplexHeatmap)

cat("\n--- Generando UpSet Plots para UP y DOWN ---\n")

# 1. Agrupar las listas de genes UP
list_up <- list(
  Amputation = up_amp,
  Cryoinjury = up_cryo,
  Ablation   = up_abla,
  Uninjured  = up_unin
)

# 2. Agrupar las listas de genes DOWN
list_down <- list(
  Amputation = down_amp,
  Cryoinjury = down_cryo,
  Ablation   = down_abla,
  Uninjured  = down_unin
)

# 3. Crear las matrices de combinación
mat_up   <- make_comb_mat(list_up)
mat_down <- make_comb_mat(list_down)

# 4. Generar y exportar el UpSet Plot de genes UP
pdf("UpSet_Plot_UP.pdf", width = 10, height = 6)
upset_up <- UpSet(mat_up, 
                  set_order = c("Amputation", "Cryoinjury", "Ablation", "Uninjured"),
                  comb_order = order(comb_size(mat_up), decreasing = TRUE),
                  pt_size = unit(5, "mm"), lwd = 3,
                  top_annotation = upset_top_annotation(mat_up, add_numbers = TRUE),
                  right_annotation = upset_right_annotation(mat_up, add_numbers = TRUE),
                  column_title = "Intersección de Genes UP-regulados")
draw(upset_up)
dev.off()

# 5. Generar y exportar el UpSet Plot de genes DOWN
pdf("UpSet_Plot_DOWN.pdf", width = 10, height = 6)
upset_down <- UpSet(mat_down, 
                    set_order = c("Amputation", "Cryoinjury", "Ablation", "Uninjured"),
                    comb_order = order(comb_size(mat_down), decreasing = TRUE),
                    pt_size = unit(5, "mm"), lwd = 3,
                    top_annotation = upset_top_annotation(mat_down, add_numbers = TRUE),
                    right_annotation = upset_right_annotation(mat_down, add_numbers = TRUE),
                    column_title = "Intersección de Genes DOWN-regulados")
draw(upset_down)
dev.off()

cat(" UpSet plots guardados exitosamente como PDFs.\n")