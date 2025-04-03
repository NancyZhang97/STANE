calculate_cell_exp<-function(spot_count,spot_nnls_output,spot_ref){
  x<-apply(spot_ref, 1, function(t) as.numeric(spot_nnls_output)*t/sum(as.numeric(spot_nnls_output)*t))
  x<-t(x)
  cell_count<-x*spot_count
  return(cell_count)
}

calculate_expected_celltype_exp<-function(filepath=getwd(),st_count,st_coords,ref_celltype_gep,loop=2,c=0.05,consider_self=T,alpha=0.2,p_cut=1e-5,max_cores=8){
  mat_neighbor_abundance<-readRDS(paste0(filepath,"/deconv_result/neighbor_abundance_loop_",loop,".rds"))
  pre_nicheDE<-read.csv(paste0(filepath,"/nicheDE_result/nicheDE_res_loop_",loop,".csv"))
  sig_beta_data<-pre_nicheDE[which(pre_nicheDE$padj<p_cut),]
  nnls_output<-readRDS(paste0(filepath,"/deconv_result/cell_deconvolute_res_loop_",loop,".rds"))
  common_gene<-intersect(row.names(st_count),row.names(ref_celltype_gep))
  st_count<-st_count[common_gene,]
  ref_celltype_gep<-ref_celltype_gep[common_gene,]

  if (max_cores > 1){
    numCores = parallel::detectCores();
    if(parallel::detectCores() > max_cores)
      numCores <- max_cores
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="")
    doParallel::registerDoParallel(cl)
    environ = c('calculate_cell_exp','generate_spot_exp','st_count','ref_celltype_gep','nnls_output','alpha','mat_neighbor_abundance','sig_beta_data')
    cell_count_list0 <- foreach::foreach(i = (1:ncol(st_count)), .export = environ) %dopar% { #
      new_ref_celltype_gep=generate_spot_exp(alpha = alpha,sample = colnames(st_count)[i],marker_ref_celltype_gep = ref_celltype_gep,
                                                    neighbor_prop = mat_neighbor_abundance,sig_beta_list = sig_beta_data)
      res<-calculate_cell_exp(spot_count=st_count[,i],spot_ref=new_ref_celltype_gep,spot_nnls_output = nnls_output[i,])
      res
    }
    cell_count_list<-list()
    for (k in colnames(ref_celltype_gep)) {
      count<-lapply(cell_count_list0,function(t) t[,k])
      count<-t(as.matrix(do.call(rbind, count)))
      row.names(count)<-row.names(st_count)
      colnames(count)<-colnames(st_count)
      cell_count_list[[k]]<-count
    }
    parallel::stopCluster(cl)

  } else {
    cell_count_list0<-list()
    for (i in (1:ncol(st_count))) {
      new_ref_celltype_gep=generate_spot_exp(alpha = alpha,sample = colnames(st_count)[i],marker_ref_celltype_gep = ref_celltype_gep,
                                             neighbor_prop = mat_neighbor_abundance,sig_beta_list = sig_beta_data)
      res<-calculate_cell_exp(spot_count=st_count[,i],spot_ref=new_ref_celltype_gep,spot_nnls_output = nnls_output[i,])
      cell_count_list0[[colnames(st_count)[i]]]<-res
    }
    cell_count_list<-list()
    for (k in colnames(ref_celltype_gep)) {
      count<-lapply(cell_count_list0,function(t) t[,k])
      count<-t(as.matrix(do.call(rbind, count)))
      row.names(count)<-row.names(st_count)
      colnames(count)<-colnames(st_count)
      cell_count_list[[k]]<-count
    }
  }
  return(cell_count_list)
}

calculate_localmoranI_bivariat<-function(input_x,input_y,weight_mtx){
  z_x<-(input_x-mean(input_x))/sd(input_x)
  z_y<-(input_y-mean(input_y))/sd(input_y)
  local_I<-z_x*(z_y%*%weight_mtx)
  x<-intersect(which(z_x<0),which(z_y%*%weight_mtx<0))
  if (length(x)>0){
    local_I[x]<-local_I[x]*(-1)
  }
  return(local_I)
}

get_correlated_LR<-function(filepath=getwd(),st_coords,cell_count_list,LR_list,nicheDE_res,keep_moranI=T,target_celltype,neighbor_celltype,p_cut=1e-5,cor_cut=1e-3,count_cut=10,max_cores=8){
  count_L<-cell_count_list[[neighbor_celltype]][intersect(unique(LR_list$L_gene),row.names(cell_count_list[[1]])),]
  count_R<-cell_count_list[[target_celltype]][intersect(unique(LR_list$R_gene),row.names(cell_count_list[[1]])),]
  count_L<-count_L[which(rowSums(count_L)>=count_cut),]
  count_R<-count_R[which(rowSums(count_R)>=count_cut),]
  LR_list<-LR_list[which(LR_list$L_gene%in%row.names(count_L)),]
  LR_list<-LR_list[which(LR_list$R_gene%in%row.names(count_R)),]
  print(paste0("Calculate correlated LR for niche effects from ",neighbor_celltype," to ",target_celltype))
  nicheDE_res<-nicheDE_res[which(nicheDE_res$interaction==paste0(target_celltype,"-",neighbor_celltype)),]
  sig_nicheDE_res<-nicheDE_res[which(nicheDE_res$padj<p_cut),]
  unsig_nicheDE_res<-nicheDE_res[which(nicheDE_res$padj>=p_cut),]

  dist_mat<-generate_dist_mtx(coord_df=st_coords,use_rank=T)
  dist_mat<-dist_mat[colnames(cell_count_list[[1]]),colnames(cell_count_list[[1]])]
  weight_mtx<-exp(-dist_mat*0.05)

  if (max_cores > 1){
    numCores = parallel::detectCores();
    if(parallel::detectCores() > max_cores)
      numCores <- max_cores
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="")
    doParallel::registerDoParallel(cl)
    environ = c('calculate_localmoranI_bivariat','LR_list','count_L','count_R','weight_mtx')
    local_moranI_list <- foreach::foreach(i = (1:nrow(LR_list)), .export = environ) %dopar% { #
      L=count_L[LR_list$L_gene[i],]
      R=count_R[LR_list$R_gene[i],]
      L=L[intersect(which(is.na(L)==F),which(is.na(R)==F))]
      R=R[intersect(which(is.na(L)==F),which(is.na(R)==F))]
      weight_mtx0<-weight_mtx[names(L),names(L)]
      local_moranI=calculate_localmoranI_bivariat(input_x = L,input_y = R,weight_mtx = weight_mtx0)
      as.numeric(local_moranI)
    }
    names(local_moranI_list)<-LR_list$LR_pair

    count_target<-cell_count_list[[target_celltype]][intersect(unique(nicheDE_res$geneName),row.names(cell_count_list[[target_celltype]])),]
    environ2 = c('local_moranI_list','count_target','nicheDE_res','sig_nicheDE_res','unsig_nicheDE_res','cor_cut')
    results <- foreach::foreach(i = (1:nrow(LR_list)), .export = environ) %dopar% { #
      local_moranI<-as.numeric(local_moranI_list[[i]])
      cor_test<-apply(count_target,1,function(t) cor.test(t,local_moranI)[c("p.value","estimate")])
      cor_test<-unlist(cor_test)
      cor_coef<-as.numeric(cor_test[paste0(nicheDE_res$geneName,".estimate.cor")])
      cor_pval<-as.numeric(cor_test[paste0(nicheDE_res$geneName,".p.value")])
      names(cor_coef)<-nicheDE_res$geneName
      names(cor_pval)<-nicheDE_res$geneName

      sign_check<-sign(cor_coef*nicheDE_res$coefficient)
      y<-names(cor_pval)[intersect(which(cor_pval<cor_cut),which(sign_check>0))]
      z<-setdiff(names(cor_pval),y)
      res<-fisher.test(matrix(c(length(intersect(y,sig_nicheDE_res$geneName)),length(intersect(y,unsig_nicheDE_res$geneName)),length(intersect(z,sig_nicheDE_res$geneName)),length(intersect(z,unsig_nicheDE_res$geneName))),nrow = 2,byrow = T))
      final_res<-c(LR_list$LR_pair[i],res$p.value,res$estimate,paste0(y,collapse=","))
      final_res
    }
    cor_LR_moranI<-as.data.frame(do.call(rbind, results))
    colnames(cor_LR_moranI)<-c("LR_pair","enrich_pval","enrich_OR","correlated_nicheDE")
    cor_LR_moranI$interaction<-paste0(target_celltype,"-",neighbor_celltype)
    parallel::stopCluster(cl)
  } else {
    local_moranI_list<-list()
    cor_LR_moranI<-list()
    for (i in (1:nrow(LR_list))) {
      L=count_L[LR_list$L_gene[i],]
      R=count_R[LR_list$R_gene[i],]
      L=L[intersect(which(is.na(L)==F),which(is.na(R)==F))]
      R=R[intersect(which(is.na(L)==F),which(is.na(R)==F))]
      weight_mtx0<-weight_mtx[names(L),names(L)]
      local_moranI=calculate_localmoranI_bivariat(input_x = L,input_y = R,weight_mtx = weight_mtx0)
      local_moranI<-as.numeric(local_moranI)
      local_moranI_list[[LR_list$LR_pair[i]]]<-local_moranI

      cor_test<-apply(count_target,1,function(t) cor.test(t,local_moranI)[c("p.value","estimate")])
      cor_test<-unlist(cor_test)
      cor_coef<-as.numeric(cor_test[paste0(nicheDE_res$geneName,".estimate.cor")])
      cor_pval<-as.numeric(cor_test[paste0(nicheDE_res$geneName,".p.value")])
      names(cor_coef)<-nicheDE_res$geneName
      names(cor_pval)<-nicheDE_res$geneName

      sign_check<-sign(cor_coef*nicheDE_res$coefficient)
      y<-names(cor_pval)[intersect(which(cor_pval<cor_cut),which(sign_check>0))]
      z<-setdiff(names(cor_pval),y)
      res<-fisher.test(matrix(c(length(intersect(y,sig_nicheDE_res$geneName)),length(intersect(y,unsig_nicheDE_res$geneName)),length(intersect(z,sig_nicheDE_res$geneName)),length(intersect(z,unsig_nicheDE_res$geneName))),nrow = 2,byrow = T))
      final_res<-c(LR_list$LR_pair[i],res$p.value,res$estimate,paste0(y,collapse=","))

      cor_LR_moranI[[i]]<-final_res
    }
    cor_LR_moranI<-as.data.frame(do.call(rbind, cor_LR_moranI))
    colnames(cor_LR_moranI)<-c("LR_pair","enrich_pval","enrich_OR","correlated_nicheDE")
    cor_LR_moranI$interaction<-paste0(target_celltype,"-",neighbor_celltype)
  }

  if (keep_moranI==T){
    saveRDS(local_moranI_list,file = paste0(filepath,"/nicheDE_result/LR_local_moranI_",neighbor_celltype,"_to_",target_celltype,".rds"))
  }
  return(cor_LR_moranI)
}
