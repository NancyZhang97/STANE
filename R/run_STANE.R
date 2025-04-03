run_STANE<-function(filepath=getwd(),st_count,st_coords,ref_celltype_gep,marker_ref_celltype_gep,cell_cut_list=NULL,loop=1,c=0.05,consider_self=T,alpha=0.2,p_cut=1e-5,prop_cut=0.1,max_cores=8){
  common_gene<-intersect(row.names(st_count),row.names(ref_celltype_gep))
  common_marker_gene<-intersect(row.names(st_count),row.names(marker_ref_celltype_gep))
  st_count<-st_count[common_gene,]
  ref_celltype_gep<-ref_celltype_gep[common_gene,]
  marker_ref_celltype_gep<-marker_ref_celltype_gep[common_marker_gene,]
  type<-colnames(ref_celltype_gep)

  if (loop == 1){
    print("Loop 1: Start the process")
    dir.create(paste0(filepath,"/deconv_result"))
    dir.create(paste0(filepath,"/nicheDE_result"))
    print("Step1: Spot-level cell type-specific deconvolution")

    results<-list()
    for(i in (1:ncol(st_count))) {
      results[[i]] <- calculate_nnls_cut(spot_count=st_count[,i],marker_ref_celltype_gep=marker_ref_celltype_gep,cell_cut_list = cell_cut_list,cutoff=prop_cut)
    }
    nnls_output = as.data.frame(do.call(rbind, results))
    row.names(nnls_output)<-colnames(st_count)
    colnames(nnls_output)<-type
  }

  else if (loop > 1){
    print(paste0("Loop ",loop,": Start the process"))
    print("Step1: Spot-level cell type-specific deconvolution")
    pre_loop<-loop-1
    mat_neighbor_abundance0<-readRDS(paste0(filepath,"/deconv_result/neighbor_abundance_loop_",pre_loop,".rds"))
    pre_nicheDE<-read.csv(paste0(filepath,"/nicheDE_result/nicheDE_res_loop_",pre_loop,".csv"))
    sig_beta_data<-pre_nicheDE[which(pre_nicheDE$padj<p_cut),]

    if(max_cores > 1) {
      numCores = parallel::detectCores();
      if(parallel::detectCores() > max_cores)
        numCores <- max_cores
      cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="")
      doParallel::registerDoParallel(cl)
      environ = c('calculate_nnls_cut','generate_spot_exp','st_count','marker_ref_celltype_gep','cell_cut_list','alpha','mat_neighbor_abundance0','sig_beta_data','prop_cut')
      results <- foreach::foreach(i = (1:ncol(st_count)), .packages = c("nnls"), .export = environ) %dopar% { #
        new_marker_ref_celltype_gep=generate_spot_exp(alpha = alpha,sample = colnames(st_count)[i],marker_ref_celltype_gep = marker_ref_celltype_gep,
                                                      neighbor_prop = mat_neighbor_abundance0,sig_beta_list = sig_beta_data)
        res<-calculate_nnls_cut(spot_count=st_count[,i],marker_ref_celltype_gep=new_marker_ref_celltype_gep,cell_cut_list = cell_cut_list,cutoff=prop_cut)
        res
      }
      parallel::stopCluster(cl)
      nnls_output = as.data.frame(do.call(rbind, results))
      row.names(nnls_output)<-colnames(st_count)
      colnames(nnls_output)<-type

    } else {
      results<-list()
      for(i in 1:(ncol(st_count))) {
        results[[i]] <- calculate_nnls_cut(spot_count=st_count[,i],marker_ref_celltype_gep=marker_ref_celltype_gep,cell_cut_list = cell_cut_list,cutoff=prop_cut)
      }
      nnls_output = as.data.frame(do.call(rbind, results))
      row.names(nnls_output)<-colnames(st_count)
      colnames(nnls_output)<-type
    }

  } else {
    print("Please input a number as parameter loop")
  }
  x<-apply(nnls_output, 1, sum)
  if (length(which(x==0))>0){
    x<-x[-which(x==0)]
    nnls_output<-nnls_output[names(x),]
  }
  prop_estimate<-nnls_output
  for (i in type) {
    prop_estimate[,i]<-nnls_output[,i]/x
  }
  prop_estimate<-as.data.frame(prop_estimate)
  st_count<-st_count[,row.names(prop_estimate)]

  print("Step2: Neighbor-cell abundance calculation")
  dist_mat<-generate_dist_mtx(coord_df=st_coords,use_rank=T)
  dist_mat<-dist_mat[row.names(prop_estimate),row.names(prop_estimate)]
  results<-list()
  for(i in 1:(ncol(st_count))) {
    results[[i]] <- neighbor_abundance(feature_mat=prop_estimate,dist_mat=dist_mat,spotID=colnames(st_count)[i],c=c,consider_self=consider_self)
  }
  mat_neighbor_abundance<-as.data.frame(do.call(rbind, results))
  row.names(mat_neighbor_abundance)<-colnames(st_count)
  colnames(mat_neighbor_abundance)<-type

  saveRDS(nnls_output,paste0(filepath,"/deconv_result/cell_deconvolute_res_loop_",loop,".rds"))
  saveRDS(mat_neighbor_abundance,paste0(filepath,"/deconv_result/neighbor_abundance_loop_",loop,".rds"))

  print("Step3: Niche effect estimation")

  if(max_cores > 1) {
    numCores = parallel::detectCores();
    if(parallel::detectCores() > max_cores)
      numCores <- max_cores
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="")
    doParallel::registerDoParallel(cl)
    environ = c('run_lasso_reg','mat_neighbor_abundance','nnls_output','st_count','ref_celltype_gep')
    res_nicheDE <- foreach::foreach(i = (1:nrow(st_count)), .packages = c("glmnet","islasso"), .export = environ) %dopar% { #
      res<-run_lasso_reg(select_gene=row.names(st_count)[i],count = st_count,ref_celltype_gep = ref_celltype_gep,neighbor_abundance_mat = mat_neighbor_abundance,cell_deconv_res = nnls_output)
      res
    }
    parallel::stopCluster(cl)
    names(res_nicheDE)<-row.names(st_count)

  } else {
    res_nicheDE<-list()
    for(i in 1:(nrow(st_count))) {
      res_nicheDE[[row.names(st_count)[i]]] <- run_lasso_reg(select_gene=row.names(st_count)[i],ref_celltype_gep = ref_celltype_gep,neighbor_abundance_mat = mat_neighbor_abundance,cell_deconv_res = nnls_output)
    }
  }

  all_reg_result<-data.frame(geneName="x",interaction=NA,coefficient=NA,pval=NA)
  for (i in (1:length(res_nicheDE))) {
    gene<-names(res_nicheDE)[i]
    res<-res_nicheDE[[i]]
    if (is.null(res)==F){
      res<-res$coefficients
      row.names(res)<-gsub("`","",row.names(res))
      if (nrow(res)>1){
        res<-as.data.frame(res)
        res<-res[-1,]
        data<-data.frame(geneName=gene,interaction=row.names(res),coefficient=res$Estimate,pval=res$`Pr(>|z|)`)
        all_reg_result<-rbind(all_reg_result,data)
      }
    }
  }
  all_reg_result_padj_sum<-all_reg_result
  all_reg_result_padj_sum$padj<-p.adjust(all_reg_result_padj_sum$pval,method = "BH")
  all_reg_result_padj_sum<-all_reg_result_padj_sum[order(all_reg_result_padj_sum$padj),]
  write.csv(all_reg_result_padj_sum,file = paste0(filepath,"/nicheDE_result/nicheDE_res_loop_",loop,".csv"),row.names = F,quote = F)
  print(paste0("Loop ",loop,": Finish the process"))
}

run_STANE_loop<-function(filepath=getwd(),st_count,st_coords,ref_celltype_gep,marker_ref_celltype_gep,cell_cut_list=NULL,start_loop=1,c=0.05,consider_self=T,alpha=0.2,p_cut=1e-5,prop_cut=0.1,max_cores=8){
  if (start_loop==1){
    run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = start_loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
    loop = start_loop+1
    run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
    rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
    print(paste0("The converge level for loop ", loop, " is ", rmse))
    while (rmse > 0.05 | loop < 6) {
      loop = loop+1
      run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
      rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
      print(paste0("The converge level for loop ", loop, " is ", rmse))
    }
  } else if (start_loop == 2){
    loop = start_loop
    run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
    rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
    print(paste0("The converge level for loop ", loop, " is ", rmse))
    while (rmse > 0.05 | loop < 6) {
      loop = loop+1
      run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
      rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
      print(paste0("The converge level for loop ", loop, " is ", rmse))
    }
  } else if (start_loop > 2){
    loop = start_loop
    rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
    print(paste0("The converge level for loop ", loop, " is ", rmse))
    while (rmse > 0.05 | loop < 6) {
      loop = loop+1
      run_STANE(filepath=filepath,st_count = st_count,st_coords = st_coords,ref_celltype_gep = ref_celltype_gep,marker_ref_celltype_gep = marker_ref_celltype_gep,cell_cut_list=cell_cut_list,loop = loop,c=c,consider_self=consider_self,alpha=alpha,p_cut=p_cut,prop_cut=prop_cut,max_cores=max_cores)
      rmse <- check_deconvolution_coverge(filepath = filepath, loop = loop)
      print(paste0("The converge level for loop ", loop, " is ", rmse))
    }
  } else {
    print("Please input a number as parameter loop")
  }

  print("Finish the interaction process and produce the final results")
  nicheDE_result<-read.csv(paste0(filepath,"/nicheDE_result/nicheDE_res_loop_",loop,".csv"))
  nnls_output<-readRDS(paste0(filepath,"/deconv_result/cell_deconvolute_res_loop_",loop,".rds"))
  prop_estimate<-nnls_output
  for (i in colnames(nnls_output)) {
    prop_estimate[,i]<-nnls_output[,i]/x
  }
  prop_estimate<-as.data.frame(prop_estimate)

  res<-list()
  res$nicheDE<-nicheDE_result
  res$propEstimate<-prop_estimate

  return(res)
}

