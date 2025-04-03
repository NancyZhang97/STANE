check_deconvolution_coverge<-function(filepath=getwd(),loop){
  pre_loop<-(loop-1)
  output_dir<-paste0(filepath,"/deconv_result")
  pre_res<-readRDS(paste0(output_dir,"/cell_deconvolute_res_loop_",pre_loop,".rds"))
  new_res<-readRDS(paste0(output_dir,"/cell_deconvolute_res_loop_",loop,".rds"))
  type=colnames(pre_res)
  common_cell<-intersect(row.names(pre_res),row.names(new_res))
  pre_res<-pre_res[common_cell,]
  new_res<-new_res[common_cell,]
  error<-c()
  for (i in type) {
    error<-c(error,pre_res[,i]-new_res[,i])
  }
  rmse<-sqrt(sum((error)^2)/length(error))

  return(rmse)
}

run_lasso_reg<-function(select_gene,count,ref_celltype_gep,neighbor_abundance_mat,cell_deconv_res){
  type<-colnames(cell_deconv_res)
  prop_estimate<-as.data.frame(cell_deconv_res)
  tryCatch(expr = {
    u_st_ept<-cell_deconv_res
    for (i in type) {
      u_st_ept[,i]<-cell_deconv_res[,i]*ref_celltype_gep[select_gene,i]
    }
    prop_estimate$u_exp<-apply(u_st_ept,1,sum)
    prop_estimate$true_exp<-count[select_gene,]
    prop_estimate$diff<-prop_estimate$true_exp-prop_estimate$u_exp
    prop_estimate$diff[which(abs(prop_estimate$diff)<0.01)]<-0

    cell_neighbor<-type
    exp_celltype<-round(ref_celltype_gep[select_gene,type],digits = 2)
    cell_target<-names(exp_celltype)[which(exp_celltype>0)]

    data<-as.data.frame(u_st_ept[row.names(neighbor_abundance_mat),cell_target,drop=F])
    data$sum<-apply(data,1,sum)
    p_st_list<-data[,cell_target,drop=F]
    if (length(which(is.nan(p_st_list[,1])))>0){
      p_st_list<-p_st_list[-which(is.nan(p_st_list[,1])),,drop=F]
    }
    cell_comb<-c()
    for (i in cell_target) {
      for (j in setdiff(cell_neighbor,i)) {
        cell_comb<-c(cell_comb,paste0(i,"-",j))
      }
    }
    neighbor_abundance_mat<-neighbor_abundance_mat[row.names(p_st_list),]
    reg_input<-matrix(data=NA,nrow = nrow(p_st_list),ncol = length(cell_comb))
    row.names(reg_input)<-row.names(p_st_list)
    colnames(reg_input)<-cell_comb
    neighbor_abundance_mat_scale<-apply(neighbor_abundance_mat, 2, scale)
    row.names(neighbor_abundance_mat_scale)<-row.names(neighbor_abundance_mat)
    for (i in cell_target) {
      for (j in setdiff(cell_neighbor,i)) {
        reg_input[,paste0(i,"-",j)]<-p_st_list[,i]*neighbor_abundance_mat_scale[,j]
      }
    }
    reg_input<-as.data.frame(reg_input)
    reg_input$y<-as.numeric(prop_estimate[row.names(reg_input),"diff"])
    reg_input$y<-as.numeric(scale(reg_input$y))
    x<-as.matrix(reg_input[,c(1:(ncol(reg_input)-1))])
    y<-as.numeric(reg_input$y)
    cv.out <- cv.glmnet(x, y, alpha = 0, nfolds = 10)
    best_lambda <- cv.out$lambda.min
    out <- islasso(y ~ ., data = reg_input, lambda = best_lambda)
    #print(select_gene)
    return(summary(out))
  }, error = function(e){message(paste("An error occurred for item", select_gene,":\n"), e)})
}
