suppressMessages(library('edgeR'))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))

counts <- read.delim(snakemake@input[['counts']], comment.char='#')
rownames(counts) <- counts$Geneid

dcounts <- counts %>% select(c(
  'deeptools_input.NPC_9sca_ctrl_RNA_rep1.bam',
  'deeptools_input.NPC_9sca_ctrl_RNA_rep2.bam',
  'deeptools_input.NPC_9sca_ko_RNA_rep1.bam',
  'deeptools_input.NPC_9sca_ko_RNA_rep2.bam'
))

genotype <- factor(c("ctrl", "ctrl", "ko", "ko"))

# Build DGEList and normalize (TMM)
dge <- DGEList(counts = dcounts, group = genotype)
dge <- calcNormFactors(dge)

# Design matrix — "ctrl" is reference level (alphabetical)
design <- model.matrix(~ genotype)

# Estimate dispersion and fit quasi-likelihood model
dge    <- estimateDisp(dge, design)
fit    <- glmQLFit(dge, design)
qlf    <- glmQLFTest(fit)

# Retrieve all genes; FDR via Benjamini-Hochberg (matches DESeq2's default)
res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")$table

# Split into down / up / non-DE  (FDR ≡ padj,  logFC ≡ log2FoldChange)
down <- res %>%
  filter(FDR  <  snakemake@params[['padj']]) %>%
  filter(logFC < -snakemake@params[['l2fc']])

up <- res %>%
  filter(FDR  <  snakemake@params[['padj']]) %>%
  filter(logFC >  snakemake@params[['l2fc']])

none <- res %>%
  filter(
    FDR >= snakemake@params[['padj']] | abs(logFC) <= snakemake@params[['l2fc']]
  )

write.table(down, snakemake@output[['down']],   sep='\t', quote=F)
write.table(up,   snakemake@output[['up']],     sep='\t', quote=F)
write.table(none, snakemake@output[['nonde']],  sep='\t', quote=F)
