library('DESeq2')
library("tidyr")
library("dplyr")

counts <- read.delim()(snakemake@input[['counts']], coment.char='#')
counts %>% colnames
rownames(counts) <- counts$Geneid
dcounts <- counts %>% select(c('deeptools_input.NPC_9sca_ctrl_RNA_rep1.bam', 'deeptools_input.NPC_9sca_ctrl_RNA_rep2.bam', 'deeptools_input.NPC_9sca_ko_RNA_rep1.bam', 'deeptools_input.NPC_9sca_ko_RNA_rep2.bam'))
metadf <- data.frame(
  'genotype' = factor(c("ctrl", "ctrl", "ko", "ko"))
)
rownames(metadf) <- dcounts %>% colnames()
metadf
dds <- DESeqDataSetFromMatrix(countData = dcounts,
                              colData = metadf,
                              design = ~ genotype)

dds <- DESeq(dds)
res <- results(dds)
up <- res %>% data.frame() %>% filter(padj < snakemake@params[['padj']]) %>% filter(log2FoldChange > snakemake@params[['l2fc']]) 
down <- res %>% data.frame() %>% filter(padj < snakemake@params[['padj']]) %>% filter(log2FoldChange < -snakemake@params[['l2fc']]) 
write.table(up, snakemake@output[['up']], sep='\t', quote=F)
write.table(down, snakemake@output[['down']], sep='\t', quote=F)