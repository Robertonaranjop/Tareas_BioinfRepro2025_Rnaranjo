###########################################
# Testing Differential Expression Testing
# ======= ============ ========== =======
#
#           Ricardo A. Verdugo
#
#                  2011
#
# Data description
# ---- -----------
# 
# This is an expression profile experiment done in the
# Illumina Mouse-Ref8 platform.
#
# Objective: To assess the effect of genetic variation
# in mouse chromosome Y on the size of cardiomyocytes.
# 
# Experimental design: Eight adult male mice from 
# two strains were profiled, C57BL/6J and 
# C57BL/6J-chrY<A/J/NaJ>, referred to as B and BY herein.
# From each strain (genotype), four animals were 
# castrated and four were sham operated. RNA was hybridized
# to Illumina MouseRef-8 v2.0 Bead-Chips that carry eight
# microarrays each, containing 25,697 probes. Only 5000 probes
# were arbitrarily selected for this tutorial. 
#
# Aims: 
#       1) To determine differential expression between genotypes
#       2) To determine differential expression between treatments
#       3) To assess differences in the response to treatment between
#          the two genotypes
#
# References:
#
# For more information about the experimental samples:
#   Llamas, Bastien, Ricardo Verdugo, Gary Churchill, and Christian Deschepper. 2009.
#   Chromosome Y variants from different inbred mouse strains are linked to differences in the
#   morphologic and molecular responses of cardiac cells to postpubertal testosterone. BMC
#   Genomics 10, no. 1 (April 7): 150. doi:10.1186/1471-2164-10-150.
# 
# For informations about the analysis of this data: Verdugo, Ricardo A., Christian F.
#   Deschepper, Gloria Munoz, Daniel Pomp, and Gary A. Churchill. 2009. Importance of
#   randomization in microarray experimental designs with Illumina platforms. Nucl. Acids Res.
#   37, no. 17 (September): 5610-8. doi:10.1093/nar/gkp573.
# 
# To speed up the tutorial, only the first 5000 of the probes in the microarray
# are included. The full dataset is available from the GEO database by id GSE15354.
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15354



# Preliminaries
# =============

source("data/install_packages.R")

## Define some constants (it is a good practice to declare them at the beginning of the scripts)
setwd("C:/Users/rnara/Desktop/Repositorio/Tareas_BioinfRepro2025_Rnaranjo/Unid4/sesion1")
outdir     <- "C:/Users/rnara/Desktop/Repositorio/Tareas_BioinfRepro2025_Rnaranjo/Unid4/sesion1/results"
fdr_th     <- 0.19 # Proportion of false discoveries that are acceptable

# Define some useful functions (it is a good idea to always document what the intend to do)
source("data/Rfxs.R")

# Create an output directory
if(!file.exists(outdir)) {
  dir.create(outdir, mode = "0755", recursive=T)
}

# Read microarray data in Illumina format
# ==== ========== ==== == ======== ======
Data.Raw  <- read.delim("data/Illum_data_sub5000.txt")
signal    <- grep("AVG_Signal", colnames(Data.Raw)) # vector of columns with data
detection <- grep("Detection.Pval", colnames(Data.Raw)) # vector of columns with p-values

# Read probe annotations
# ==== ===== ===========
annot     <- read.delim("data/MouseRef-8_annot_full.txt")
# Make sure that probes are in same order as in the dataset
annot <- annot[match(Data.Raw$PROBE_ID, annot$ProbeID), ]
# Not all probes have the same quality from sequence alignment to the genome.
table(annot$ProbeQuality)
# We will group 'Bad' with 'No match' as bad and everything else as good.
# For details see Nucleic Acids Res. 38:e17
probe_qc <- ifelse(annot$ProbeQuality %in% c("Bad", "No match"), "Bad probes", "Good probes")

# Read hybridization design
# ==== ============= ======
design <- read.csv("data/YChrom_design.csv", as.is=F)
# Lets take a look
print(design)

# Quality Control
# ======= =======

# Boxplots
palette(rainbow(4))
alabel <- sprintf("Array%02i", 1:length(signal))
afact  <- factor(rep(alabel, each=nrow(Data.Raw)))
qcfact <- factor(rep(probe_qc, length(signal)))

# Color by genotype
png(file.path(outdir,"boxplot_raw_probe_qc.png"), width=6, height=3, unit="in", res=150)
  par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7, las=3)
  boxplot(unlist(log2(Data.Raw[,signal]))~qcfact+afact, horiz=T, main="Raw log2 values Boxplot",
          col=rep(1:2, length(signal)), axes=F, varwidth=TRUE, ylab="log2(intensity)")
  axis(1, at=seq(1, length(signal)*2, by=2)+.5, labels=alabel)
  axis(2)
  legend("top", legend=levels(qcfact), fill=1:2, ncol=2, xjust=.5, bty="n", inset=-.15)
dev.off()

# Color by treatment
png(file.path(outdir,"boxplot_raw_treatment.png"), width=6, height=3, unit="in", res=150)
  par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7)
  boxplot(as.data.frame(log2(Data.Raw[,signal])), horiz=T, main="Raw log2 values Boxplot", las=1,
          col=design$Treatment, names=design$Sentrix_Position, cex.axis=.9,
          xlab="Sentrix position", ylab="log2(intensity)")
  legend("top", legend=levels(design$Treatment), fill=1:2, ncol=2, xjust=.5, bty="n", inset=-.15)
dev.off()

# Scatter plots of raw data in log2 scale.
png(file.path(outdir,"Pairs_scatter_log2.png"), width=8, height=8, unit="in", res=150)
  par(cex=.2, mar=c(2.1,2.1,2.1,1.1))
  pairs(log2(Data.Raw[,signal]), main="Log2 Raw Intensity Values", pch=".",  gap=.5, cex.labels=.5)
dev.off()

# Probe filtering by QC
# ===== ========= == ==
# Since probe with bad score tend to have a lower signal than good probes, it is 
# recommended to remove them. 4706 probes are left in the dataset.
Data.Raw <- Data.Raw[probe_qc %in% "Good probes",]
annot    <- annot[probe_qc %in% "Good probes",]

# Create matrix of raw data
# ====== ====== == === ====
rawdata           <- as.matrix(Data.Raw[,signal])
rownames(rawdata) <- Data.Raw$PROBE_ID
colnames(rawdata) <- design$Sample_Name


# Data Normalization
# ==== =============

# Load functions for normalization

library(preprocessCore)

normdata           <- normalize.quantiles(rawdata) 
colnames(normdata) <- colnames(rawdata)
rownames(normdata) <- rownames(rawdata)


# Probe Filtering
# ===== =========
# This step aims to removing probes that did not detect a transcript
# in any of the experimental groups. Note that this step can be optional.
#
# Create a vector or P/A calls for each probe using detection probabilities calculated by BeadStudio
probe_present      <- Data.Raw[,detection] < 0.04
detected_per_group <- t(apply(probe_present, 1, tapply, design$Group, sum))
present  <- apply(detected_per_group >= 1, 1, all)
normdata <- normdata[present,]
annot    <- annot[present, ]

# Testing for differential expression
# ======= === ============ ==========

library(limma)

# 1) Verificación básica: columnas en el mismo orden que el diseño
stopifnot(all(colnames(normdata) == design$Sample_Name))

# 2) Medias y errores estándar por grupo (equivalente a lo que hacías con madata$data)
Means <- t(apply(normdata, 1, tapply, design$Group, mean))
SEs   <- t(apply(normdata, 1, tapply, design$Group, function(x) sqrt(var(x)/length(x))))

colnames(Means) <- paste("Mean", colnames(Means), sep=":")
colnames(SEs)   <- paste("SE", colnames(SEs), sep=":")

# 3) Modelo por grupo (4 grupos)
design_limma <- model.matrix(~ 0 + design$Group)
colnames(design_limma) <- levels(design$Group)

##
# =========================
# Contrasts definition
# =========================
# Order of groups MUST be:
# B.C  B.I  BY.C  BY.I

cmat <- rbind(
  Geno     = c( 1,  1, -1, -1 ) * 0.5,
  Trt      = c( 1, -1,  1, -1 ) * 0.5,
  Int      = c( 1, -1, -1,  1 ),
  Geno_I   = c( 0,  1,  0, -1 ),
  Geno_C   = c( 1,  0, -1,  0 ),
  Trt_B    = c( 1, -1,  0,  0 ),
  Trt_BY   = c( 0,  0,  1, -1 ),
  B.C_BY.I = c( 1,  0,  0, -1 ),
  B.I_BY.C = c( 0,  1, -1,  0 )
)

# Nombres de columnas = niveles del factor Group
colnames(cmat) <- levels(design$Group)


# 4) Matriz de contrastes: conservar tu cmat, pero alineada al orden del diseño
cmat2 <- cmat[, colnames(design_limma)]

# 5) Ajuste limma + contrasts
fit  <- lmFit(normdata, design_limma)

C <- t(cmat2)
dim(C)
fit2 <- contrasts.fit(fit, C)

fit2 <- eBayes(fit2)

##
colnames(design_limma)
colnames(cmat2)      # deberían ser iguales
rownames(C)          # deberían ser los 4 grupos
colnames(C)          # deberían ser los contrastes (Geno, Trt, Int, ...)


# 6) Fold-change (equivalente a logDiffs y logdiff2FC del tutorial)
#    - logDiffs: matriz (probes x contrasts)
logDiffs <- Means %*% t(cmat2)
FC <- apply(logDiffs, 2, logdiff2FC)

# 7) Armar tabla de resultados "tipo tutorial"
results <- data.frame(annot, Means, SEs)

# Guardar por contraste: P y FDR (BH), y también logFC de limma
for (nm in colnames(C)) {
  tt <- topTable(fit2, coef = nm, number = Inf, adjust.method = "BH", sort.by = "none")
  results[[paste0("P.", nm)]]     <- tt$P.Value
  results[[paste0("FDR.", nm)]]   <- tt$adj.P.Val
  results[[paste0("logFC.", nm)]] <- tt$logFC
}

# Agregar FC (como en el tutorial: una columna por contraste)
# (FC ya tiene nombres de columnas = nombres de contrastes)
for (nm in colnames(FC)) {
  results[[paste0("FC.", nm)]] <- FC[, nm]
}

# 8) Exportar resultados
write.table(results, file=file.path(outdir,"DE_results.csv"), sep=",", row.names=FALSE)

# 9) (Opcional pero útil) Histograma de p-values para los 3 contrastes principales (Geno/Trt/Int)
if (all(c("P.Geno","P.Trt","P.Int") %in% names(results))) {
  png(file.path(outdir,"P-values_Hist.png"), width=6, height=4, unit="in", res=150)
  par(mfrow=c(1,3), cex=.8)
  hist(results$P.Geno, main="P-values: Geno", xlab="P", breaks=30)
  hist(results$P.Trt,  main="P-values: Trt",  xlab="P", breaks=30)
  hist(results$P.Int,  main="P-values: Int",  xlab="P", breaks=30)
  dev.off()
}

# =========================
# Count genes + Venn diagrams
# =========================

# Crear GeneID como en el tutorial (Entrez si existe, si no ProbeID)
results$GeneID <- results$EntrezID
results$GeneID[is.na(results$GeneID)] <- results$ProbeID[is.na(results$GeneID)]

# Probes seleccionadas por DE (FDR <= 0.19)
Probes.DE <- results[, c("FDR.Geno", "FDR.Trt", "FDR.Int")] <= fdr_th

# Regla pedida: un gen se selecciona SOLO si TODAS sus sondas están seleccionadas
Genes.DE <- apply(Probes.DE, 2, tapply, results$GeneID, all)

library(limma)

Counts.DE <- vennCounts(Genes.DE)
print(Counts.DE)

png(file.path(outdir, "vennDiagram_DiffExprs.png"), width=3.5, height=3, unit="in", res=150)
par(cex=.7)
vennDiagram(Counts.DE, names=c("Geno", "Trt", "Int"),
            main="\n\n\nDifferentially Expressed Genes")
dev.off()

# --- Interacción: particionar por Geno_I/Geno_C y Trt_B/Trt_BY (solo en genes con Int significativo)
idx_int <- results$FDR.Int <= fdr_th

Probes.Int_Geno <- results[idx_int, c("FDR.Geno_I", "FDR.Geno_C")] <= fdr_th
Genes.Int_Geno  <- apply(Probes.Int_Geno, 2, tapply, results$GeneID[idx_int], all)

Probes.Int_Trt <- results[idx_int, c("FDR.Trt_B", "FDR.Trt_BY")] <= fdr_th
Genes.Int_Trt  <- apply(Probes.Int_Trt,  2, tapply, results$GeneID[idx_int], all)

Counts.Int_Geno <- vennCounts(Genes.Int_Geno)
Counts.Int_Trt  <- vennCounts(Genes.Int_Trt)

print(Counts.Int_Geno)
print(Counts.Int_Trt)

png(file.path(outdir, "vennDiagram_Int.png"), width=6.5, height=3, unit="in", res=150)
par(mfrow=c(1,2), cex=.7)

vennDiagram(Counts.Int_Geno, names=c("I", "C"),
            main="\n\n\nGenes Responding to Genotype\nin a Treatment Dependent Manner")

vennDiagram(Counts.Int_Trt, names=c("B", "BY"),
            main="\n\n\nGenes Responding to Treatment\nin a Genotype Dependent Manner")
dev.off()

# =========================
# Functional testing (topGO)
# =========================

# topGO necesita EntrezID (saca sondas sin EntrezID)
results_go <- results[!is.na(results$EntrezID), ]

# Genes seleccionados por interacción (a nivel PROBE)
probes.int <- results_go$FDR.Int <= fdr_th

# Regla pedida: gen seleccionado solo si TODAS sus sondas están seleccionadas
genes.int <- tapply(probes.int, results_go$EntrezID, all)

# topGO espera un factor de 0/1
genes.int <- ifelse(genes.int, 1, 0)
genes.int <- as.factor(genes.int)

#######
# 1) Asegura BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2) Instala dependencias que suelen faltar en Windows
install.packages(c("DBI", "RSQLite"), dependencies = TRUE)

# 3) Instala AnnotationDbi (Bioconductor) que usa RSQLite
BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)

# 4) Ahora instala los annotation packages + topGO
BiocManager::install(c("GO.db", "org.Mm.eg.db", "topGO"),
                     ask = FALSE, update = FALSE)

####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")


#########

library(topGO)
library(org.Mm.eg.db)

GOdata <- new("topGOdata",
              ontology="BP",
              description="Genes DE by Trt by Geno Interaction",
              allGenes=genes.int,
              nodeSize=5,
              annotationFun=annFUN.org,
              mapping="org.Mm.eg.db",
              ID="entrez")

resultFisher.classic <- runTest(GOdata, algorithm="classic", statistic="fisher")
resultFisher.elim    <- runTest(GOdata, algorithm="elim", statistic="fisher")

GO_BP_Table <- GenTable(GOdata,
                        Fisher.classic=resultFisher.classic,
                        Fisher.elim=resultFisher.elim,
                        orderBy="Fisher.elim",
                        ranksOf="Fisher.classic",
                        topNodes=20)

print(GO_BP_Table)

write.table(GO_BP_Table, file.path(outdir, "GO_BP_Table.csv"),
            sep=",", row.names=FALSE)


