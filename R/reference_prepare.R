generate_ref_celltype_gep<-function(sc_count,label,normalize=T){
  sc<-CreateSeuratObject(sc_count)
  sc$celltype<-sc_label
  type<-as.character(unique(label))
  ref_celltype_gep<-matrix(nrow = nrow(sc),ncol = length(type))
  row.names(ref_celltype_gep)<-row.names(sc)
  colnames(ref_celltype_gep)<-type
  if (normalize==T){
    sc<-NormalizeData(sc)
    for (i in type) {
      ref_celltype_gep[,i]<-apply(sc@assays$RNA$data[,which(sc$celltype==i)],1,mean)
    }
  } else {
    for (i in type) {
      ref_celltype_gep[,i]<-apply(sc@assays$RNA$counts[,which(sc$celltype==i)],1,mean)
    }
  }

  return(ref_celltype_gep)
}

generate_marker_celltype_gep<-function(ref_celltype_gep,sc_count,label,max_p=1e-30,top_count=100){
  sc<-CreateSeuratObject(sc_count)
  sc$celltype<-sc_label
  type<-as.character(unique(label))
  Idents(sc)<-sc$celltype
  sc<-NormalizeData(sc)
  markers<-FindAllMarkers(sc,only.pos = T)
  select_markers<-markers %>% group_by(cluster) %>% top_n(n = top_count, wt = avg_log2FC)
  select_markers<-select_markers[which(select_markers$p_val_adj<max_p),]
  marker_ref_celltype_gep<-ref_celltype_gep[select_markers$gene,]

  return(marker_ref_celltype_gep)
}

