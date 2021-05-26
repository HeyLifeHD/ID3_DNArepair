#Anand GO 
library(KEGG.db)
library(goseq)

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))

DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#functions
goseq_wrapper = function(assayedGenes, deGenes, bn = NULL){
  gene.vector=as.integer(assayedGenes %in% deGenes)
  names(gene.vector)= assayedGenes
  
  pwf = goseq::nullp(gene.vector,"hg19","ensGene", plot.fit = FALSE)
  GO = goseq::goseq(pwf,"hg19","ensGene", test.cats = c("GO:BP", "GO:MF")#, method = "Hypergeometric"
  )
  
  GO.kegg = goseq::goseq(pwf,"hg19","ensGene", test.cats = 'KEGG')
  rownames(GO.kegg) = GO.kegg$category
  GO.kegg$pathName = unlist(mget(x = as.character(GO.kegg$category), envir = KEGGPATHID2NAME, ifnotfound = NA))
  GO.kegg = GO.kegg[order(GO.kegg$over_represented_pvalue),]
  GO.kegg$fdr = p.adjust(GO.kegg$over_represented_pvalue)
  
  GO$fdr = p.adjust(GO$over_represented_pvalue)
  
  if(!is.null(bn)){
    write.table(GO, paste0(bn, '_GO.xls'), sep='\t', quote = FALSE, row.names = FALSE)
    write.table(GO.kegg, paste0(bn, '_KEGG.xls'), sep='\t', quote = FALSE, row.names = FALSE)
  }
  #return(GO)
  return(list(go = GO, kegg = GO.kegg))
}


#This function plots top ten GO process from GOSEq results
plotgo = function(gores, top = c(20, 10), bn = NULL){
  gores.go = gores$go
  gores.go$term = factor(x = gores.go$term, levels = rev(gores.go$term))
  g.gg = ggplot(data = gores.go[1:top[1],], aes(x = term, y = -log10(fdr)))+geom_bar(stat = 'identity')+coord_flip()+cowplot::theme_cowplot(font_size = 12, line_size = 1.3)+cowplot::background_grid(major = 'xy')+xlab('')
  gores.kegg = gores$kegg
  gores.kegg$pathName = factor(x = gores.kegg$pathName, levels = rev(gores.kegg$pathName))
  k.gg = ggplot(data = gores.kegg[1:top[2],], aes(x = pathName, y = -log10(fdr)))+geom_bar(stat = 'identity')+coord_flip()+cowplot::theme_cowplot(font_size = 12, line_size = 1.3)+cowplot::background_grid(major = 'xy')+xlab('')
  
  gg = cowplot::plot_grid(g.gg, k.gg, nrow = 2, rel_heights = c(2, 1), labels = c('GO', 'KEGG'))
  
  if(!is.null(bn)){
    write.table(gores.go, file = paste0(bn, '_GO.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    colnames(gores.kegg)[c(6,7)] = c('KEGG_ID', 'Pathway')
    write.table(gores.kegg, file = paste0(bn, '_KEGG.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    cowplot::save_plot(filename = paste0(bn, '_GO.pdf'), plot = gg, base_height = 7, base_width = 8)
  }
  
  print(gg)
}

#analysis
cutoff <- 0.05
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff),]
  x <- x$ensembl
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]

res <- list()
for (i in names(genes2plot)){
  res[[i]] <- goseq_wrapper(DEG_results_list[[1]]$ensembl, genes2plot[[i]])
  plotgo(res[[i]], bn= paste0(PostDE.dir, "/", i,"/OR"))
}

##GSEA of DNA damage
dna_damage <- as.character(read.table(file.path(data.dir, "DNArepairgenelist.csv"))$V1)