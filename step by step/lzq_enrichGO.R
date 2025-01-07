lzq_enrichGO <- function(gene,qvalueCutoff,GO2GENE = GO2GENE,GOdetial=GOdetial){
  x <- enricher(gene,minGSSize = 1,maxGSSize = 50000,
                pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = qvalueCutoff,TERM2GENE = GO2GENE)
  x@result <- x@result%>%
    dplyr::filter(qvalue<qvalueCutoff,pvalue<0.05)
  x@result <- merge(GOdetial,x@result[,-2],by=1)%>%
    dplyr::rename('Description'='Term','ID'='go_id','ONTOLOGY'='Ontology')%>%
    arrange(p.adjust)
  return(x)
}