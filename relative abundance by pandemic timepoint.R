library(DESeq2)
 
physeqR4 <- qza_to_phyloseq(features="Rtable10tp4.qza", tree="recrooted-tree.qza", taxonomy="Rtaxonomy.qza", metadata="R4map.tsv")
 
library(qiime2R)
 
alpha = 0.01
 
physeqR4 <- tax_glom(physeqR4, taxrank = "Genus")
 
R4covid <- phyloseq_to_deseq2(physeqR4, ~ batch)
 
R4covid <- estimateSizeFactors(R4covid, type="poscount")
 
R4covid <- DESeq(R4covid)
 
resR4covid <- results(R4covid)
 
sigtabR4covid <- resR4covid[which(resR4covid$padj < alpha), ]
 
sigtabR4covid = cbind(as(sigtabR4covid, "data.frame"), as(tax_table(physeqR4)[rownames(sigtabR4covid), ], "matrix"))
 
View(sigtabR4covid)
