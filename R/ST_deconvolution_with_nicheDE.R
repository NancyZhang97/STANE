generate_dist_mtx<-function(coord_df,use_rank=T){
  dist_mat<-pointDistance(coord_df, lonlat=FALSE)
  row.names(dist_mat)<-row.names(coord_df)
  colnames(dist_mat)<-row.names(coord_df)
  dist_mtx<-dist_mat
  if (use_rank==T){
    rank_mat<-round(dist_mat)
    rank_mat<-apply(rank_mat, 1, dense_rank)
    row.names(rank_mat)<-row.names(coord_df)
    dist_mtx<-rank_mat
  }

  return(dist_mtx)
}

neighbor_abundance<-function(feature_mat,dist_mat,spotID,c=1,consider_self=T){
  x<-dist_mat[spotID,]
  x<-x[intersect(names(x),row.names(feature_mat))]
  neighbor_ID<-names(x)
  if (consider_self==F){
    neighbor_ID<-setdiff(neighbor_ID,spotID)
  }
  prop_mat<-feature_mat[neighbor_ID,]
  output_abundance<-apply(prop_mat,2,function (t) weighted.mean(t,exp(-c*x)))
  return(output_abundance)
}

generate_spot_exp<-function(alpha=0.2,sample,sig_beta_list,marker_ref_celltype_gep,neighbor_prop){
  type<-colnames(marker_ref_celltype_gep)
  cellcomb<-c()
  for (i in type) {
    for (j in setdiff(type,i)) {
      cellcomb<-c(cellcomb,paste0(i,"-",j))
    }
  }
  new_marker_ref_celltype_count<-marker_ref_celltype_gep
  sig_beta_mtx<-matrix(data = 0, nrow = nrow(marker_ref_celltype_gep),ncol = length(cellcomb))
  row.names(sig_beta_mtx)<-row.names(marker_ref_celltype_gep)
  colnames(sig_beta_mtx)<-cellcomb
  for (i in (1:nrow(sig_beta_list))) {
    if (sig_beta_list$geneName[i]%in%row.names(sig_beta_mtx)){
      sig_beta_mtx[sig_beta_list$geneName[i],sig_beta_list$interaction[i]]<-sig_beta_list$coefficient[i]
    }
  }
  for (i in type) {
    data<-sig_beta_mtx[,paste0(i,"-",setdiff(type,i))]
    #data<-apply(data, 2, function(x) ifelse(x==0,0,1/x))
    data<-apply(data, 2, function(x) ifelse(x>0,1,ifelse(x<0,-1,0)))
    for (j in setdiff(type,i)) {
      data[,paste0(i,"-",j)]<-data[,paste0(i,"-",j)]*neighbor_prop[sample,j]
    }
    niche_w<-apply(data, 1, sum)
    new_marker_ref_celltype_count[,i]<-new_marker_ref_celltype_count[,i]*(1+alpha*niche_w)
  }
  return(new_marker_ref_celltype_count)
}

determine_cell_exist<-function(prop_vec,cutoff=0.1,cell_cut_list=NULL){
  cell_exist<-prop_vec
  if (is.null(cell_cut_list)==F){
    for (i in names(prop_vec)) {
      cell_exist[i]<-ifelse(prop_vec[i]>cell_cut_list[i],1,0)
    }
  } else {
    for (i in names(prop_vec)) {
      cell_exist[i]<-ifelse(prop_vec[i]>cutoff,1,0)
    }
  }

  return(cell_exist)
}

calculate_nnls_cut<-function(spot_count,marker_ref_celltype_gep,cutoff=0.1,cell_cut_list=NULL){
  type<-colnames(marker_ref_celltype_gep)
  y<-spot_count[row.names(marker_ref_celltype_gep)]
  mod1 <- nnls(marker_ref_celltype_gep,y)
  x<-mod1$x
  x<-x/sum(x)
  names(x)<-type
  exist_vector<-x
  if (is.null(cell_cut_list)==F){
    for (i in names(x)) {
      exist_vector[i]<-ifelse(x[i]>cell_cut_list[i],1,0)
    }
  } else {
    for (i in names(x)) {
      exist_vector[i]<-ifelse(x[i]>cutoff,1,0)
    }
  }
  data0<-marker_ref_celltype_gep[,which(exist_vector==1),drop=F]
  mod1 <- nnls(data0,y)
  y<-rep(0,length(type))
  y[which(exist_vector==1)]<-mod1$x
  names(y)<-type
  return(y)
}
