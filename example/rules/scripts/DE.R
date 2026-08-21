suppressMessages(library('edgeR'))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
set.seed(123)

padj_cutoff <- as.double(snakemake@params[['padj']])
l2fc_cutoff <- as.double(snakemake@params[['l2fc']])

counts <- read.delim(snakemake@input[['counts']], comment.char='#')
rownames(counts) <- counts$Geneid

dcounts <- counts %>% select(c(
  'deeptools_input.NPC_9sca_ctrl_RNA_rep1.bam',
  'deeptools_input.NPC_9sca_ctrl_RNA_rep2.bam',
  'deeptools_input.NPC_9sca_ko_RNA_rep1.bam',
  'deeptools_input.NPC_9sca_ko_RNA_rep2.bam'
))

genotype <- factor(c("ctrl", "ctrl", "ko", "ko"))
keep <- filterByExpr(dcounts, group = genotype)
dcounts <- dcounts[keep, ]

dge <- DGEList(counts = dcounts, group = genotype)
dge <- calcNormFactors(dge)

design <- model.matrix(~ genotype)
dge    <- estimateDisp(dge, design)
fit    <- glmQLFit(dge, design)
qlf    <- glmQLFTest(fit)

res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")$table

coef_name <- colnames(design)[2]
shrunk <- predFC(dge, design, prior.count = 5, dispersion = dge$trended.dispersion)
res$logFC <- shrunk[rownames(res), coef_name]

down <- res %>%
  filter(FDR  <  padj_cutoff) %>%
  filter(logFC < -l2fc_cutoff)

up <- res %>%
  filter(FDR  <  padj_cutoff) %>%
  filter(logFC >  l2fc_cutoff)

none <- res %>%
  filter(FDR > padj_cutoff) %>%
  filter(abs(logFC) < l2fc_cutoff) %>%
  slice_sample(n = 1000)

write.table(down, snakemake@output[['down']],   sep='\t', quote=F)
write.table(up,   snakemake@output[['up']],     sep='\t', quote=F)
write.table(none, snakemake@output[['nonde']],  sep='\t', quote=F)
write.table(res,  snakemake@output[['res']],    sep='\t', quote=F)
