
factors <- c("ENSG00000117713", "ENSG00000117318", "ENSG00000137337")
names(factors)<- mapIds(org.Hs.eg.db, keys=factors, keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )
rownames(dds_genes) <- sapply(strsplit(rownames(dds_genes) ,".", fixed=TRUE),`[`, 1)
dir.create(file.path(base_results.dir,"Genecounts"))

for(gene in 1:length(factors)){ 
  Count <- plotCounts(dds_genes, gene = gene,intgroup = c("group", "genotype", "replicate", "group_uncomplete2"), returnData = TRUE, normalized=TRUE)
  
  pdf(file.path(base_results.dir,"Genecounts",paste0(names(factors[gene]),".pdf")), height = 10, width = 15)
  print(ggboxplot(Count, x="group_uncomplete2", y="count", color="genotype", facet.by = "genotype",merge=FALSE,combine=TRUE,yscale = "log2", 
            xlab="Genotype", ylab=("Normalized counts"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
}



fpkm <- readRDS(file.path(base.dir,"data", "FPKM.rds"))
rownames(fpkm) <- sapply(strsplit(rownames(fpkm) ,".", fixed=TRUE),`[`, 1)

for(gene in 1:length(factors)){ 
  Count <- as.data.frame(t(fpkm[gene, ]))
  colnames(Count)<- c("count")
  Count$sample_name <- rownames(Count)
  Count <- left_join(Count, as.data.frame(colData(dds_genes)))

  pdf(file.path(base_results.dir,"Genecounts",paste0(names(factors[gene]),"_FPKM.pdf")), height = 10, width = 15)
  print(ggboxplot(Count, x="group_uncomplete2", y="count", color="genotype", facet.by = "genotype",merge=FALSE,combine=TRUE,#yscale = "log2", 
                  xlab="Genotype", ylab=("Normalized counts"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
}

for(gene in 1:length(factors)){ 
  Count <- as.data.frame(t(fpkm[gene, ]))
  colnames(Count)<- c("count")
  Count$sample_name <- rownames(Count)
  Count <- left_join(Count, as.data.frame(colData(dds_genes)))
  
  pdf(file.path(base_results.dir,"Genecounts",paste0(names(factors[gene]),"_per_sample_FPKM.pdf")), height = 10, width = 15)
  print(ggboxplot(Count, x="sample_name", y="count", color="genotype", merge=FALSE,combine=TRUE,#yscale = "log2", 
                  xlab="sample_name", ylab=("Normalized counts"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
}

