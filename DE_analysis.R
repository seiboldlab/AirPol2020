library(DESeq2)

### load count matrix
cm <- read.table("AirPol_raw_counts.txt", sep='\t', header=T)

### load design table
design <- read.table("AirPol_design.txt", sep='\t', header=T)
rownames(design) <- design$Full_ID
design$ID <- factor(design$ID)

### filter out genes with low expression
cm_filt <- cm[which(rowSums(cm > 6) >= 10),]

### filter out MT and ribosomal genes
cm_filt <- cm_filt[-grep("MTAT|MT-|MTCO|MTCY|MTERF|MTND|MTRF|MTRN|MRPL|RMRP|MRPS|RPL|RPS",rownames(cm_filt)),]

donors_low_and_high <- c("5337", "5609", "5613", "5633", "5697")
donors_NIST <- c("5337", "5609", "5613", "5633")

#==============================================================#
# Differential expression analysis between low OE dose vs Ctrl #
#==============================================================#

design_low <- design[design$Treatment %in% c("0.75OE", "OECtrl") & design$ID %in% donors_low_and_high,]
design_low$Treatment <- factor(design_low$Treatment, levels=c("OECtrl", "0.75OE"))
cm_low <- cm_filt[,rownames(design_low)]

dds_low <- DESeqDataSetFromMatrix(
       countData = cm_low,
       colData = design_low,
       design = ~ ID + Treatment)

dds_low <- DESeq(dds_low)
res_dds_low <- results(dds_low)
res_dds_low <- res_dds_low[order(res_dds_low$padj),]

low_dose_DEGs <- res_dds_low[which(!is.na(res_dds_low$padj) & res_dds_low$padj < 0.05),]

#===================================================================#
# Differential expression analysis between moderate OE dose vs Ctrl #
#===================================================================#

design_mod <- design[design$Treatment %in% c("7.5OE", "OECtrl"),]
design_mod$Treatment <- factor(design_mod$Treatment, levels=c("OECtrl", "7.5OE"))
cm_mod <- cm_filt[,rownames(design_mod)]

dds_mod <- DESeqDataSetFromMatrix(
       countData = cm_mod,
       colData = design_mod,
       design = ~ ID + Treatment)

dds_mod <- DESeq(dds_mod)
res_dds_mod <- results(dds_mod)
res_dds_mod <- res_dds_mod[order(res_dds_mod$padj),]

mod_dose_DEGs <- res_dds_mod[which(!is.na(res_dds_mod$padj) & res_dds_mod$padj < 0.05),]

#===============================================================#
# Differential expression analysis between high OE dose vs Ctrl #
#===============================================================#

design_high <- design[design$Treatment %in% c("75OE", "OECtrl") & design$ID %in% donors_low_and_high,]
design_high$Treatment <- factor(design_high$Treatment, levels=c("OECtrl", "75OE"))
cm_high <- cm_filt[,rownames(design_high)]

dds_high <- DESeqDataSetFromMatrix(
       countData = cm_high,
       colData = design_high,
       design = ~ ID + Treatment)

dds_high <- DESeq(dds_high)
res_dds_high <- results(dds_high)
res_dds_high <- res_dds_high[order(res_dds_high$padj),]

high_dose_DEGs <- res_dds_high[which(!is.na(res_dds_high$padj) & res_dds_high$padj < 0.05),]

#=====================================================================#
# Differential expression analysis between NIST vs Water Extract Ctrl #
#=====================================================================#

design_NIST <- design[design$Treatment %in% c("NIST", "WECtrl") & design$ID %in% donors_NIST,]
design_NIST$Treatment <- factor(design_NIST$Treatment, levels=c("WECtrl", "NIST"))
cm_NIST <- cm_filt[,rownames(design_NIST)]

dds_NIST <- DESeqDataSetFromMatrix(
       countData = cm_NIST,
       colData = design_NIST,
       design = ~ ID + Treatment)

dds_NIST <- DESeq(dds_NIST)
res_dds_NIST <- results(dds_NIST)
res_dds_NIST <- res_dds_NIST[order(res_dds_NIST$padj),]

NIST_DEGs <- res_dds_NIST[which(!is.na(res_dds_NIST$padj) & res_dds_NIST$padj < 0.05),]

