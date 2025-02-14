output_HVG_mtx=function(input_mtx,top_N,gene_filter=NULL){
  if(is.null(gene_filter)){gene_filter=rownames(input_mtx)}
  scaled_input_mtx=t(scale(t(input_mtx)))
  Mad_values=apply(scaled_input_mtx,1,mad)
  ordered_all_express_Genes=rownames(scaled_input_mtx)[order(Mad_values,decreasing=T)]
  ordered_all_express_query_Genes=ordered_all_express_Genes[ordered_all_express_Genes%in%gene_filter]
  output_mtx=input_mtx[ordered_all_express_query_Genes[1:top_N],]
  return(output_mtx)
}



run_GraphCluster=function(mtx,r=0.8,k=20){
  library(Seurat)
  d=nrow(mtx)
  singlecell_dt=as.data.frame(t(mtx))
  sim_embed=singlecell_dt
  singlecell_dt$CellID=rownames(sim_embed)
  rownames(sim_embed)=singlecell_dt$CellID
  colnames(sim_embed)=paste0("TOY_",1:ncol(sim_embed))
  sim_reduc<- CreateDimReducObject(
    embeddings = as.matrix(sim_embed),
    key = "TOY_",
    assay = "RNA")
  tmpcount=matrix(data = 1,nrow = 1000,ncol = nrow(singlecell_dt),
                  dimnames = list(paste0("G",1:1000),singlecell_dt$CellID))
  testdt=CreateSeuratObject(counts = tmpcount)
  testdt@reductions$toy=sim_reduc
  testdt <- FindNeighbors(testdt, reduction = "toy", dims = 1:d,k.param = k,verbose = F)
  Obj=FindClusters(testdt,resolution = r,verbose = F)
  return(Obj@active.ident)
}

ConsensusGraphCluster=function(input_matrix,r=0.8,k=20,reps=1000,pItem=0.9){
  library(reshape2)
  n=ncol(input_matrix)
  ClusterResults=lapply(1:reps,function(x){
    set.seed(x)
    subset_matrix = input_matrix[,sample(1:n,round(n*pItem),replace = F)]
    cl=run_GraphCluster(subset_matrix,r=r,k=k)
    cl=factor(cl,labels = paste0("C",levels(cl)))
    return(cl)
  })
  
  ##
  cv_cl_num=unlist(lapply(lapply(ClusterResults, unique),length))
  
  #idx=which(cv_cl_num==target_cl)
  
  valid_recluster=ClusterResults#[idx]
  message(reps," Cluster stat:")
  print(cv_cl_num)
  message("Target Cluster Test:")
  print(length(valid_recluster))
  if(length(valid_recluster)==0){
    CC_mtx=NA
    
  }else{
    ClusterResults_mtx=lapply(1:length(valid_recluster), function(n){
      subtype=ClusterResults[[n]]
      mtx=do.call(rbind,lapply(subtype, function(x){return(subtype==x)}))
      colnames(mtx)=names(subtype)
      mtx2=ifelse(mtx,1,0)
      mtx3=melt(mtx2)
      mtx3$rep=n
      return(mtx3)
    })
    
    rep_consensus_df=do.call(rbind,ClusterResults_mtx)
    rep_consensus_df$Pair=apply(rep_consensus_df[,1:2],1,function(x) return(paste0(x,collapse = "|")))
    count.list=split(rep_consensus_df,rep_consensus_df$Pair)
    cc_ratio=lapply(count.list, function(x){
      tmp=factor(as.character(x[,3]),levels = c("0","1"))
      out=prop.table(table(tmp))
      return(out)
    })
    cc_ratio_df=do.call(rbind,cc_ratio)
    cc_ratio_df0=as.data.frame(cc_ratio_df)
    cc_ratio_df1=data.frame(do.call(rbind,strsplit(rownames(cc_ratio_df),"\\|")),ratio=cc_ratio_df0$`1`) 
    CC_mtx=dcast(cc_ratio_df1,X1~X2) %>% tibble::column_to_rownames("X1")
  }
  return(list(ClusterResults,CC_mtx))
}

filterLowExprGenes=function(TPM,lowExpTh=2,lowExpRatio=0.95){
  EXP=TPM
  zero_count=apply(EXP, 1, function(x){ return(length(which(x<=lowExpTh)))})
  EXP=EXP[-which(zero_count>=(ncol(EXP)*lowExpRatio)),]
  message("filtered data:",paste(dim(EXP),collapse =  "*"))
  return(EXP)
  
}


cal_median_silhouette=function(m,clusterResult){
  n=length(unique(clusterResult))
  if(n==1){
    Median_silhouette=NA
  }else{
    all_silhouette <- cluster::silhouette(as.numeric(clusterResult), dist(t(m),method = "euclidean"))
    Median_silhouette=median(all_silhouette[,3])
  }
return(Median_silhouette)
}

cal_calinhara=function(m,clusterResult){
  n=length(unique(clusterResult))
  if(n==1){
    RNAsubtype_calinhara=NA
  }else{
    RNAsubtype_calinhara <- fpc::calinhara(t(m),as.numeric(clusterResult) )
  }
return(RNAsubtype_calinhara)
}

plot_cluster_metrics_score=function(x,xlab_title,ylab_title,plot_title){
  df=data.frame(score=as.numeric(x),index=as.numeric(names(x)))
  p=ggplot(df,aes(x=index,y=score))+geom_line()+theme_light(base_size = 10)+xlab(xlab_title)+ylab(ylab_title)+ggtitle(plot_title)
return(p)
}
