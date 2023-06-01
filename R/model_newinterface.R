
.drLouvain <- function(hZ, verbose=TRUE){
  ### Louvain cluster based on low-dimensional embeddings.

  require(Seurat)
  n <- nrow(hZ); q <- ncol(hZ)
  row.names(hZ) <- paste0("spot", 1:n)
  colnames(hZ) <- paste0("gene", 1:q)
  seu <- CreateSeuratObject(counts= t(hZ), assay='RNA')
  DefaultAssay(seu) <- "RNA"
  pca1 <- CreateDimReducObject(embeddings = hZ, key = "PC_", assay='RNA')
  seu@reductions$"pca" <- pca1
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:q, verbose=verbose)
  seu <- FindClusters(seu, verbose=verbose)
  return(seu$seurat_clusters)
}

getInitialK <- function(Z, verbose=TRUE){
 ## Use Louvain to obatin initail clusters number to speed up computation.
    clusters <- .drLouvain(hZ=Z, verbose=verbose)
    K <- length(unique(clusters))
    return(K)
}


#' Set parameters for SC-MEB model
#' @description  Prepare parameters setup for SC-MEB model fitting.
#'
#' @param beta_grid a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0,4,0.2).
#' @param coreNum an integer value to decide how many cores are used to run SC-MEB model in parallel.
#' @param maxIter_ICM the maximum iteration of ICM step, default as 6.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set.
#' @export
#' @examples
#' model_set_SCMEB(coreNum=5)
#'
model_set_SCMEB <- function(beta_grid=seq(0.2, 4, by=0.2),
                            coreNum=1, maxIter_ICM=6, maxIter=30,
                            verbose=TRUE, seed=1){

  para_settings <- list(
    coreNum = coreNum,
    beta_grid= beta_grid,
    maxIter_ICM=maxIter_ICM,maxIter= maxIter,
    verbose=verbose, seed=seed)
  return(para_settings)

}




#' Set parameters for DR-SC model
#' @description  Prepare parameters setup for DR-SC model fitting.
#'
#' @param error_heter a logical value, whether use the heterogenous error for DR-SC model, default as TRUE. If error.heter=FALSE, then the homogenuous error is used for probabilistic PCA model in DR-SC.
#' @param wpca_int an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as FALSE which means the ordinary PCA is used.
#' @param int.model an optional string, specify which Gaussian mixture model is used in evaluting the initial values for DR-SC, default as "EEE"; and see Mclust for more models' names.
#' @param approxPCA an optional logical value, whether use approximated PCA to speed up the computation for initial values.
#' @param beta_grid a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0.2,4,0.2).
#' @param coreNum an integer value to decide how many cores are used to run in parallel.
#' @param maxIter_ICM the maximum iteration of ICM step, default as 6.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 25
#' @param epsLogLik an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set.
#' @export
#' @examples
#' model_set_DRSC(approxPCA=TRUE, coreNum=5)
#'
model_set_DRSC <- function(error_heter= TRUE,
                           wpca_int=FALSE, int.model="EEE", approxPCA=FALSE,
                           beta_grid=seq(0.2, 4, by=0.2), coreNum = 5,
                           maxIter_ICM=6, maxIter=25, epsLogLik=1e-5,
                           verbose=FALSE, seed=1){

  para_settings <- list(
                        error_heter=error_heter, wpca_int=wpca_int,int.model=int.model,
                        coreNum = coreNum, approxPCA = approxPCA,
                        beta_grid= beta_grid,
                        maxIter_ICM=maxIter_ICM,maxIter= maxIter, epsLogLik=epsLogLik,
                        verbose=verbose, seed=seed)
  return(para_settings)

}

#' Set parameters for PRECAST model
#' @description  Prepare parameters setup for PRECAST model fitting.
#'
#' @param Sigma_equal a logical value, whether set the mixture covariance matrices equal, default as FALSE.
#' @param Sigma_diag a logical value, whether set the mixture covariance matrices diagonal, default as TRUE.
#' @param mix_prop_heter a logical value, whether set the smoothing parameter of each data batch to be different, default as TRUE.
#' @param error_heter a logical value, whether use the heterogenous error for PRECAST model, default as TRUE.
#' @param Sp2 a logical value, whether add the intrisical CAR component in the PRECAST model, default as TRUE.
#' @param wpca_int an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as FALSE which means the ordinary PCA is used.
#' @param int.model an optional string, specify which Gaussian mixture model is used in evaluting the initial values for DR-SC, default as "EEE"; and see Mclust for more models' names.
#' @param beta_grid a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0.2,4,0.2).
#' @param coreNum an integer value to decide how many cores are used to run  in parallel, default as 5.
#' @param coreNum_int an integer value to decide how many cores are used in parallel computing initial values.
#' @param maxIter_ICM the maximum iteration of ICM step, default as 6.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 25
#' @param epsLogLik an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set.
#' @export
#' @examples
#' model_set_PRECAST(Sigma_equal=TRUE, coreNum=5)
#'
model_set_PRECAST <- function(Sigma_equal=FALSE, Sigma_diag=TRUE,mix_prop_heter=TRUE,
                              error_heter=TRUE, Sp2=TRUE, wpca_int=FALSE,int.model='EEE',
                              coreNum = 5, coreNum_int=coreNum,
                              beta_grid=seq(0.2,4, by=0.2),
                              maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE, seed=1){
  para_settings <- list(Sigma_equal=Sigma_equal, Sigma_diag=Sigma_diag,
                        mix_prop_heter=mix_prop_heter,
                        error_heter=error_heter, Sp2=Sp2, wpca_int=wpca_int,int.model=int.model,
                        coreNum = coreNum, coreNum_int=coreNum_int,
                        beta_grid= beta_grid,
                        maxIter_ICM=maxIter_ICM,maxIter= maxIter, epsLogLik=epsLogLik,
                        verbose=verbose, seed=seed)
  return(para_settings)

}

#' Set parameters for iSC.MEB model
#' @description  Prepare parameters setup for iSC.MEB model fitting.
#'
#' @param Sigma_equal a logical value, whether set the mixture covariance matrices equal, default as FALSE.
#' @param Sigma_diag a logical value, whether set the mixture covariance matrices diagonal, default as TRUE.
#' @param int.model an optional string, specify which Gaussian mixture model is used in evaluting the initial values for DR-SC, default as "EEE"; and see Mclust for more models' names.
#' @param beta_grid a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0.2,4,0.2).
#' @param coreNum an integer value to decide how many cores are used to run  in parallel, default as 5.
#' @param maxIter_ICM the maximum iteration of ICM step, default as 6.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 25
#' @param epsLogLik an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set.
#' @param init.start a positive integer, the number of times to compute initial values.
#' @param c_penalty a positive real, the penalty constant used in modified BIC to determine the number of clusters.
#' @export
#' @examples
#' model_set_iSCMEB(Sigma_equal=TRUE, coreNum=5)
#'
model_set_iSCMEB <- function(Sigma_equal = FALSE,
                             Sigma_diag = TRUE,
                             int.model = "EEE",
                             beta_grid = seq(0, 5, by = 0.2),
                             coreNum = 1,
                             maxIter_ICM = 6,
                             maxIter = 25,
                             epsLogLik = 1e-05,
                             verbose = TRUE,
                             seed = 1,
                             init.start = 1,
                             c_penalty = 1){

  para_settings <- list(Sigma_equal = Sigma_equal,
                        Sigma_diag = Sigma_diag,
                        beta_grid = beta_grid,
                        maxIter_ICM = maxIter_ICM,
                        maxIter = maxIter,
                        epsLogLik = epsLogLik,
                        int.model = int.model,
                        init.start = init.start,
                        coreNum = coreNum,
                        verbose = verbose,
                        seed = seed,
                        c_penalty = c_penalty)
  return(para_settings)
}

SC.MEB_structrue <- function(Z, K, Adj, parameterList=NULL){


  if(is.null(parameterList)){
    parameterList <- model_set_SCMEB()
  }
  ### initialize variables
  beta_grid <- maxIter_ICM<- maxIter<-  verbose<-  NULL
  seed<- coreNum <-  NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }
  set.seed(seed)
  prallel <- TRUE
  if(length(K)==1){
    prallel <- FALSE
  }


  resList <- SC.MEB(
    Z,
    Adj_sp=Adj,
    beta_grid = beta_grid,
    K_set = K,
    parallel = prallel,
    num_core = coreNum,
    PX = TRUE,
    maxIter_ICM = maxIter_ICM,
    maxIter = maxIter
  )

  return(resList)
}
DR.SC_structure <- function(XList,  K, AdjList, q=15, parameterList=NULL){

  if(length(XList)>1 || length(AdjList)>1)
    stop("DR.SC_structure: There are more than one data batch in XList or AdjList, in which DR.SC model can not be fitted\n")
  if(is.null(parameterList)){
    parameterList <- model_set_DRSC()
  }
  ### initialize variables
  beta_grid <- maxIter_ICM<- maxIter<- epsLogLik<- verbose<-  NULL
  error_heter<- wpca_int<-int.model<- approxPCA<-seed<- coreNum <-  NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }
  set.seed(seed)
  resList <- DR.SC::DR.SC_fit(X=XList[[1]], q=q, K=K, Adj_sp=AdjList[[1]],
                    beta_grid=beta_grid,
                    maxIter_ICM=maxIter_ICM,maxIter=maxIter,
                    epsLogLik=epsLogLik, verbose=verbose,
                    error.heter=error_heter, approxPCA=approxPCA,
                    wpca.int=wpca_int,int.model=int.model,
                    coreNum = coreNum)

  return(resList)
}

iSCMEB_structure <- function(VList, AdjList, K, parameterList=NULL){

  #require(iSC.MEB)
  if(is.null(parameterList)){
    parameterList <- model_set_iSCMEB()
  }
  ### initialize variables
  beta_grid <- maxIter_ICM<- maxIter<-  epsLogLik <- verbose<-  NULL
  int.model <- init.start<- Sigma_equal<-  Sigma_diag<- seed<- coreNum <-c_penalty <-  NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }
  set.seed(seed)
  if(length(K)==1){
    coreNum <- 1
  }


  resList <- fit.iscmeb(
    VList,
    AdjList,
    K,
    Sigma_equal = Sigma_equal,
    Sigma_diag = Sigma_diag,
    beta_grid = beta_grid,
    maxIter_ICM = maxIter_ICM,
    maxIter = maxIter,
    epsLogLik = epsLogLik,
    int.model = int.model,
    init.start = init.start,
    coreNum = coreNum,
    verbose = verbose,
    seed = seed,
    c_penalty = c_penalty,
    criteria = "MBIC"
  )

  return(resList)
}
# Addmodel <- function(SRTProj, model=c("PRECAST", "DR-SC", "SC-MEB", "iSC.MEB"), ...){
#
#   prefix <- getSRToutputPrefix()
#   verbose <- getSRTVerbose()
#   model <- match.arg(model)
#
#   if(verbose){
#     message(prefix, "Add ", model, "  model!")
#   }
#
#   NumSamples <- nrow(SRTProj@sampleMetaData)
#   if(toupper(model) == "PRECAST"){
#
#     list_set <- model_set_PRECAST(...)
#     SRTProj@models[["PRECAST"]][["modelSettings"]] <- list_set
#
#   }else if(toupper(model) == "DR-SC"){
#     if(NumSamples>1)
#       stop("Addmodel: there must be only one data bathes in SRTProj to add DR-SC model!\n You can use PRECAST or iSC.MEB model")
#     list_set <- model_set_DRSC(...)
#     SRTProj@models[["DR-SC"]][["modelSettings"]] <- list_set
#   }else if(toupper(model) == 'SC-MEB'){
#     if(NumSamples>1)
#       stop("Addmodel: there must be only one data bathes in SRTProj to add SC-MEB model! \n You can use PRECAST or iSC.MEB model")
#     list_set <- model_set_SCMEB(...)
#     SRTProj@models[["SC-MEB"]][["modelSettings"]] <- list_set
#   }
#   dm <- attr(SRTProj@models, "DefaultModel")
#   if(is.null(dm)){ ##  if there is no default model, set default model to model
#     DefaultModel(SRTProj) <- model
#     if(verbose)
#       message(prefix, "Set default model to ", model, "  model; users can change it by using DefaultModel(obj)<- 'xx'")
#   }
#
#   return(SRTProj)
# }
#
#
# DefaultModel <- function(object){
#
#    dm <- attr(object@models, "DefaultModel")
#    return(dm)
# }
#
#
# `DefaultModel<-` <-function (object, value) {
#
#   #object <- SRTProj
#   attr(object@models, "DefaultModel") <- value
#   return(object)
# }

#' Joint embedding, clustering and alignment for SRT data by fitting a PRECAST model
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param K An optional integer or integer vector, specify the candidates of number of clusters. if K=NULL, it will be set to 4~12.
#' @param q An optional integer, specify the number of low-dimensional embeddings to extract in PRECAST.
#' @param ... other parameters pass to \link{model_set_PRECAST}.
#' @references \href{https://www.nature.com/articles/s41467-023-35947-w}{Wei Liu, Liao, X., Luo, Z. et al, Jin Liu* (2023). Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST. Nature Communications, 14, 296}
#' @export
#' @importFrom PRECAST ICM.EM_structure selectModel
#'
#'
Integrate_PRECAST <- function(SRTProj,  K=NULL, q= 15, ...){

  ## Get the global settings
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()


  tstart <- Sys.time()
  ## Get variable features:
  var.features <- row.names(SRTProj@geneMetaData)[SRTProj@geneMetaData[, 1]]
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')

  XList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[var.features,]
    Matrix::t(data)
  })

  AdjList <- .mylapply(seq_along(samplenames), function(id){
    .getSparseMatrixFromH5file(hfile, groupname = paste0('AdjMat/AdjMat_', id))
  })
  o <- h5closeAll()





  XList <- lapply(XList, scale, scale=FALSE)
  modelSettings <-  model_set_PRECAST(...)
  if(length(K)==1){
    modelSettings$coreNum <- 1
    modelSettings$coreNum_int <- 1
  }
  .logDiffTime(sprintf(paste0("%s Read normalized data matrices and adjacency matrices from ", SRTProj@projectMetadata$projectName,
                              ".h5 file"), prefix), t1 = tstart, verbose = verbose)

  .logTime("Start fitting PRECAST model...", prefix, verbose)
  tstart <- Sys.time()
  resList <- PRECAST::ICM.EM_structure(XList, K=K, q=q, AdjList = AdjList,
                                       parameterList = modelSettings)

  ## clear object to reduce memory burden
  rm(XList, AdjList)
  ## Write the results to path rds file
  outputPath <- SRTProj@projectMetadata$outputPath
  filename <- SRTProj@projectMetadata$filename
  file_output <- paste0(outputPath, "/", filename)
  if(!dir.exists(file_output)){
    dir.create(file_output)
  }
  .logDiffTime(sprintf(paste0("%s Finish model fitting"), prefix), t1 = tstart, verbose = verbose)
  .logTime("Start writing results of PRECAST to disk...", prefix, verbose)
  tstart <- Sys.time()
  save(resList, file= paste0(file_output, "/resPRECAST.rds"))
  .logDiffTime(sprintf(paste0("%s Finish results writing"), prefix), t1 = tstart, verbose = verbose)

  ## select the best Model
  reslist <- PRECAST::SelectModel(resList)


  SRTProj@reductions[["microEnv.PRECAST"]] <-  matlist2mat(reslist$hV)
  SRTProj@reductions[["aligned.PRECAST"]] <- matlist2mat(reslist$hZ)
  SRTProj@clusters <- factor(unlist(reslist$cluster))
  SRTProj@models$PRECAST[["Rf"]] <- reslist$Rf


  return(SRTProj)
}

#' Joint clustering and alignment for SRT data by fitting an iSC.MEB model
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param K an optional integer or integer vector, specify the candidates of number of clusters. if K=NULL, it will be set to 4~12.
#' @param reduction  a string, specify the reduction slot name in SRTProj for fitting iSC.MEB model.
#' @param ... other parameters pass to \link{model_set_iSCMEB}.
#' @references \href{https://doi.org/10.1093/bioadv/vbad019}{ iSC.MEB: an R package for multi-sample spatial clustering analysis of spatial transcriptomics data, Bioinformatics Advances, 2023, vbad019}
#' @export
#' @importFrom iSC.MEB fit.iscmeb
#'
#'
Integrate_iSCMEB <- function(SRTProj, K= NULL, reduction='PCA', ...){

  if(!require("iSC.MEB", quietly = TRUE)){
    if(require("remotes", quietly = TRUE)){
      install_github("XiaoZhangryy/iSC.MEB")
    }else{
      stop("Integrate_iSCMEB: there is no iSC.MEB package installed!")
    }
  }

  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()

  if(is.null(K)){
    if(verbose)
      message("Integrate_iSCMEB: when K is NULL, Louvain is used to search a range for candidates of number of clusters!")
    hK <- getInitialK(SRTProj@reductions[[reduction]], verbose)
    K <- seq(max(4, hK-3), min(hK+3, 20), by = 1)

  }
  if(verbose){
    message(paste0("Integrate_iSCMEB: The candidates of number of clusters are: ", paste(K, collapse=' ')))
  }

  modelSettings <- model_set_iSCMEB(...)
  if(length(K)==1){
    modelSettings$coreNum <- 1
  }
  # Read adjacence matrix
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')
  AdjList <- .mylapply(seq_along(samplenames), function(id){
    .getSparseMatrixFromH5file(hfile, groupname = paste0('AdjMat/AdjMat_', id))
  })
  o <- h5closeAll()

  .logTime("Start fitting iSC.MEB model...", prefix, verbose)
  tstart <- Sys.time()
  Z <- SRTProj@reductions[[reduction]]
  VList <- mat2list(Z, nvec=SRTProj@sampleMetaData$NumOfSpots)
  rm(Z)
  resList <- iSCMEB_structure(VList,  AdjList = AdjList,K=K,
                              parameterList = modelSettings)

  .logDiffTime(sprintf(paste0("%s Finish model fitting"), prefix), t1 = tstart, verbose = verbose)
  .logTime("Start writing results of iSC.MEB to disk...", prefix, verbose)
  ## Write the results to path rds file
  outputPath <- SRTProj@projectMetadata$outputPath
  filename <- SRTProj@projectMetadata$filename
  file_output <- paste0(outputPath, "/", filename)
  if(!dir.exists(file_output)){
    dir.create(file_output)
  }
  save(resList, file= paste0(file_output, "/resiSCMEB.rds"))
  .logDiffTime(sprintf(paste0("%s Finish results writing"), prefix), t1 = tstart, verbose = verbose)

  opt <- resList@paramList$opt

  SRTProj@reductions[["microEnv.iSC.MEB"]] <-  matlist2mat(resList@fitList[[opt]]$hV)
  SRTProj@reductions[["aligned.iSC.MEB"]] <-  matlist2mat(resList@reduction$iSCMEB)
  SRTProj@clusters <- as.factor(unlist(resList@idents))
  SRTProj@models$iSC.MEB[["Rf"]] <- resList@fitList[[opt]]$Rf ## biological information


  return(SRTProj)

}


#' Joint embedding and clustering for SRT data by fitting an iSC.MEB model a DR-SC model
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param K An optional integer or integer vector, specify the candidates of number of clusters. if K=NULL, it will be set to 4~12.
#' @param q An optional integer, specify the number of low-dimensional embeddings to extract in DR-SC model.
#' @param ... other parameters pass to \link{model_set_DRSC}.
#' @references \href{https://academic.oup.com/nar/article/50/12/e72/6555431}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi and Jin Liu* (2022). Joint dimension reduction and clustering analysis of single-cell RNA-seq and spatial transcriptomics data, Nucleic Acids Research, gkac219}
#' @export
#' @importFrom DR.SC selectModel DR.SC_fit
#'
Cluster_DRSC <- function(SRTProj,  K=NULL, q= 15, ...){


  verbose <- getSRTVerbose()
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')
  if(length(samplenames)>1) stop("Cluster_DRSC: can not applied to multiple data batches!")
  ## Get variable features:
  var.features <- row.names(SRTProj@geneMetaData)[SRTProj@geneMetaData[, 1]]
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  XList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[var.features,]
    Matrix::t(data)
  })
  AdjList <- .mylapply(seq_along(samplenames), function(id){
    .getSparseMatrixFromH5file(hfile, groupname = paste0('AdjMat/AdjMat_', id))
  })

  o <- h5closeAll()



  XList <- lapply(XList, scale, scale=FALSE)

  modelSettings <- model_set_DRSC(verbose=verbose,...)
  if(length(K)==1){
    modelSettings$coreNum <- 1
  }

  resList <- DR.SC_structure(XList, K=K, q=q, AdjList = AdjList,
                                       parameterList = modelSettings)
  ## Write the results to path rds file
  outputPath <- SRTProj@projectMetadata$outputPath
  filename <- SRTProj@projectMetadata$filename
  file_output <- paste0(outputPath, "/", filename)
  if(!dir.exists(file_output)){
    dir.create(file_output)
  }
  save(resList, file= paste0(file_output, "/resDR-SC.rds"))

  ## select the best Model
  reslist <- DR.SC::selectModel(resList)



  SRTProj@reductions[["DR.SC"]] <- reslist$hZ
  SRTProj@clusters <- factor(reslist$cluster)
  # SRTProj@models$PRECAST[["R"]] <- reslist$R

  return(SRTProj)
}

#' Spatial clustering for SRT data by fitting a SC-MEB model
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param K an optional integer or integer vector, specify the candidates of number of clusters. if K=NULL, it will be set to 4~12.
#' @param reduction  a string, specify the reduction slot name in SRTProj for fitting SC-MEB model.
#' @param ... other parameters pass to \link{model_set_iSCMEB}.
#' @references \href{https://academic.oup.com/bib/article/23/1/bbab466/6440124}{Yi Yang, Xingjie Shi, Wei Liu, et.al and Jin Liu* (2021). SC-MEB: spatial clustering with hidden Markov random field using empirical Bayes. Briefings in Bioinformatics}
#' @export
#' @importFrom SC.MEB selectK SC.MEB
#'
Cluster_SCMEB <- function(SRTProj,  K=NULL, reduction='PCA', ...){


  # require(SC.MEB)
  verbose <- getSRTVerbose()
  if(is.null(K)){
    if(verbose)
     message("SRTfitSCMEB: when K is NULL, Louvain is used to search a range for candidates of number of clusters!")
    hK <- getInitialK(SRTProj@reductions[[reduction]], verbose)
    K <- seq(max(4, hK-3), min(hK+3, 20), by = 1)

  }
  if(verbose){
    message(paste0("SRTfitSCMEB: The candidates of number of clusters are: ", paste(K, collapse=' ')))
  }

  modelSettings <- model_set_SCMEB(...)
  if(length(K)==1){
    modelSettings$coreNum <- 1
  }
  # Read adjacence matrix
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')
  if(length(samplenames)>1 )
    stop("SRTfitSCMEB: There are more than one data batch in SRTProj, in which SC.MEB model can not be fitted\n")

  AdjList <- .mylapply(seq_along(samplenames), function(id){
    .getSparseMatrixFromH5file(hfile, groupname = paste0('AdjMat/AdjMat_', id))
  })
  o <- h5closeAll()
  Adj <- AdjList[[1]]



  Z <- SRTProj@reductions[[reduction]]

  resList <- SC.MEB_structrue(Z, K=K,  Adj = Adj,
                             parameterList = modelSettings)
  ## Write the results to path rds file
  outputPath <- SRTProj@projectMetadata$outputPath
  filename <- SRTProj@projectMetadata$filename
  file_output <- paste0(outputPath, "/", filename)
  if(!dir.exists(file_output)){
    dir.create(file_output)
  }
  save(resList, file= paste0(file_output, "/resSC-MEB.rds"))

  ## select the best Model
  reslist <- SC.MEB::selectK(resList, K_set = K)
  SRTProj@clusters <- as.factor(reslist$best_K_label)


  return(SRTProj)
}

.get_correct_mean_exp <- function(XList,  hVList, hW, features){

  r_max <- length(XList)
  X0 <- XList[[1]][,features]
  hV0 <- hVList[[1]]
  if(r_max>1){
    for(r in 2:r_max){
      X0 <- rbind(X0, XList[[r]][,features])
      hV0 <- rbind(hV0, hVList[[r]])
    }
  }

  as.matrix(X0 - hV0%*% base::t(hW))

}
#' Set parameters for Louvain model
#' @description  Prepare parameters setup for Louvain clustering.
#'
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param method Method for running leiden (defaults to matrix which is fast for small datasets). Enable method = "igraph" to avoid casting large data to a dense matrix.
#' @param n.start Number of random starts.
#' @param n.iter Maximal number of iterations per random start.
#' @param seed Seed of the random number generator.
#' @param group.singletons Group singletons into nearest cluster. If FALSE, assign all singletons to a "singleton" group.
#' @export
#' @examples
#' model_set_Louvain(n.start=5)
#'
model_set_Louvain <- function(algorithm=1, method='matrix', n.start=10, n.iter=10, seed=1, group.singletons=TRUE){

  para_settings <- list(
    algorithm=algorithm, method=method,n.start=n.start,
    n.iter = n.iter, seed = seed,
    group.singletons= group.singletons)
  return(para_settings)
}

#' Louvain clustering for a SRT data
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param resolution an positive real, specify the resolution.
#' @param reduction  a string, specify the reduction slot name in SRTProj for fitting Louvain model.
#' @param ... other parameters pass to \link{model_set_Louvain}.
#' @export
#' @importFrom Seurat FindNeighbors FindClusters
#'
Cluster_Louvain <- function(SRTProj, resolution=0.8, reduction="PCA", ...){
  ### Louvain cluster based on low-dimensional embeddings.

  # require(Seurat)
  verbose <- getSRTVerbose()

  parameterList <- model_set_Louvain(...)
  ### initialize variables
  algorithm <- method<- n.start<-n.iter<- seed<-group.singletons <- NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }

  set.seed(seed)


  hZ <- SRTProj@reductions[[reduction]]
  n <- nrow(hZ); q <- ncol(hZ)
  row.names(hZ) <- paste0("spot", 1:n)
  colnames(hZ) <- paste0("gene", 1:q)
  seu <- CreateSeuratObject(counts= t(hZ), assay='RNA')

  DefaultAssay(seu) <- "RNA"
  pca1 <- CreateDimReducObject(embeddings = hZ, key = "PC_", assay='RNA')
  rm(hZ)
  seu@reductions$"pca" <- pca1
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:q, verbose=verbose)
  seu <- FindClusters(seu, verbose=verbose, resolution =resolution,
                      algorithm=algorithm, method=method,n.start=n.start,
                      n.iter = n.iter, seed = seed,
                      group.singletons= group.singletons)
  SRTProj@clusters <- factor(as.numeric(seu$seurat_clusters))
  return(SRTProj)
}



#' Set parameters for SpatialAnno model
#' @description  Prepare parameters setup for SpatialAnno model fitting.
#'
#' @param error_heter a logical value, whether use the heterogenous error for SpatialAnno model, default as TRUE. If error.heter=FALSE, then the homogenuous error is used for probabilistic PCA model in DR-SC.
#' @param wpca_int an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as FALSE which means the ordinary PCA is used.
#' @param int.model an optional string, specify model for computing the initial values for SpatialAnno (one of "SCINA" and "scSorter"), default as "SCINA".
#' @param approxPCA an optional logical value, whether use approximated PCA to speed up the computation for initial values.
#' @param beta_grid a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0.2,4,0.2).
#' @param maxIter_ICM the maximum iteration of ICM step, default as 6.
#' @param maxIter the maximum iteration of ICM-EM algorithm. The default is 25
#' @param epsLogLik an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed a postive integer, the random seed to be set.
#' @param add_Unknown a logical value, whether allow "Unknown" in the identified cell types, default as TRUE.
#' @export
#' @examples
#' model_set_SpatialAnno(int.model='scSorter')
#'
model_set_SpatialAnno <- function(error_heter = TRUE, wpca_int=TRUE, int.model='SCINA',
                                  beta_grid = seq(0.1, 2.5, by = 0.2),
                                  maxIter_ICM=6, maxIter=30, epsLogLik=1e-5,
                                   Sigma_diag=FALSE,
                                  add_Unknown = TRUE, verbose=TRUE, seed = 1){

  para_settings <- list(
    error_heter = error_heter, add_Unknown = add_Unknown, beta_grid = beta_grid,
    maxIter_ICM=maxIter_ICM, maxIter= maxIter, epsLogLik=epsLogLik,
    wpca_int=wpca_int, int.model=int.model, Sigma_diag=Sigma_diag,
     seed = seed, verbose=verbose)
  return(para_settings)

}

SpatialAnno_structure <- function(XList, markers, AdjList, q=15, parameterList){

  # require(SpatialAnno)

  #icmem <-  getFromNamespace("icmem", ns ='SpatialAnno')
  if(length(XList)>1 || length(AdjList)>1)
    stop("SpatialAnno_structure: There are more than one data batch in XList or AdjList, in which SpatialAnno model can not be fitted\n")
  X <- XList[[1]]; Adj_sp <- AdjList[[1]]
  ### initialize variables
  error_heter<-add_Unknown <- beta_grid <- maxIter_ICM<- maxIter<- epsLogLik<-  NULL
  wpca_int<-int.model<- Sigma_diag<-seed<-  verbose <- NULL
  n_par <- length(parameterList)
  for(i in 1:n_par){
    assign(names(parameterList)[i], parameterList[[i]])
  }
  set.seed(seed)

  rho <- marker_list_to_mat(markers, add_Unknown)
  df_all = as.data.frame(matrix(0, 0, 2))
  colnames(df_all) = c("Marker", "Type")
    for (i in 1:(dim(rho)[2] - 1)) {
      # print(rownames(rho)[rho[, i] == 1])
      df = as.data.frame(rownames(rho)[rho[, i] == 1])
      colnames(df) = "Marker"
      df$Type = colnames(rho)[i]
      df_all = rbind(df_all, df)
    }
    anno = df_all[, c(2, 1)]
    anno_processed = scSorter:::design_matrix_builder(anno, weight = 2)
    dat <- scSorter:::data_preprocess(t(X), anno_processed)
    if (add_Unknown == TRUE)
      dat$designmat$Unknown = 0
    rho <- as.matrix(dat$designmat)
    m <- nrow(rho)
    X_m <- t(dat$dat[1:m, ])
    X_u <- t(dat$dat[-(1:m), ])
    n <- nrow(X_m)
    K <- ncol(rho)
    if (int.model == "SCINA") {
      results = adjSCINA(t(X), markers, max_iter = 100, convergence_n = 10,
                         rm_overlap = 0, allow_unknown = TRUE)
      bet_min <- numeric(length(results$theta))
      bet_max <- numeric(length(results$theta))
      alpha_int <- numeric(m)
      bet_int <- matrix(0, nrow = m, ncol = K)
      for (i in 1:length(results$theta)) {
        bet_max[i] <- max(results$theta[[i]]$mean[, 1] -
                            results$theta[[i]]$mean[, 2])
        bet_min[i] <- min(results$theta[[i]]$mean[, 1] -
                            results$theta[[i]]$mean[, 2])
        bet_int[rho[, i] != 0, i] <- results$theta[[i]]$mean[,
                                                             1] - results$theta[[i]]$mean[, 2]
        alpha_int[rho[, i] != 0] <- results$theta[[i]]$mean[,
                                                            2]
      }
      lfc <- median(bet_min)
      colnames(rho)[K] = "Unknown"
      y_hat <- results$cell_labels
      y_hat[y_hat == "unknown"] = "Unknown"
      if (length(unique(y_hat)) < K) {
        idx = match(names(markers), unique(y_hat))
        idx2 = which(is.na(idx) == TRUE)
        if (length(idx2) > 0) {
          y_hat[1:length(names(markers)[idx2])] = names(markers)[idx2]
        }
        if (is.na(match("Unknown", unique(y_hat)))) {
          y_hat[length(names(markers)[idx2]) + 1] = "Unknown"
        }
      }
      y_int <- match(y_hat, colnames(rho))
      R_int <- matrix(0, nrow = n, ncol = K)
      for (i in 1:nrow(X_m)) R_int[i, match(y_hat, colnames(rho))[i]] <- 1
      Pi_u_int <- colMeans(R_int)
      sigma_int <- update_sigma(X_m, rho, R_int, alpha_int,
                                bet_int)
      princ <- wpca(X_u, q = q, weighted = wpca_int)
      Lam_u_int <- princ$Lam_vec
      W_u_int <- princ$loadings
      hZ <- princ$PCs
      n_c <- colSums(R_int)
      Mu_u_int <- t(sapply(1:K, function(k) 1/n_c[k] * colSums(R_int[,
                                                                     k] * hZ)))
      q <- ncol(Mu_u_int)
      Sgm_u_int <- init_Sgm(R_int, hZ, matrix(0, ncol = q,
                                              nrow = q), Mu_u_int, FALSE)
    }
    else if (int.model == "scSorter") {
      library(scSorter)
      dat <- t(X)
      rts <- scSorter(dat, anno)
      lfc = 0
      y_hat <- rts$Pred_Type
      if (length(unique(y_hat)) < K) {
        idx = match(names(markers), unique(y_hat))
        idx2 = which(is.na(idx) == TRUE)
        if (length(idx2) > 0) {
          y_hat[1:length(names(markers)[idx2])] = names(markers)[idx2]
        }
        if (is.na(match("Unknown", unique(y_hat)))) {
          y_hat[length(names(markers)[idx2]) + 1] = "Unknown"
        }
      }
      y_int <- match(y_hat, colnames(rho))
      R_int <- matrix(0, nrow = n, ncol = K)
      for (i in 1:n) R_int[i, match(y_hat, colnames(rho))[i]] <- 1
      idx_colsum_R_int <- which(colSums(R_int) == 0)
      if (length(idx_colsum_R_int) != 0) {
        for (i in 1:length(idx_colsum_R_int)) {
          R_int[i, idx_colsum_R_int[i]] = 1
        }
      }
      Pi_u_int <- colMeans(R_int)
      bet_int <- rts$Pred_param[1:m, 1:(K - 1)] * rho[, (K -
                                                           1)]
      bet_int <- cbind(bet_int, matrix(0, m, 1))
      alpha_int <- apply((rts$Pred_param[1:m, 1:(K - 1)] -
                            bet_int[, 1:(K - 1)]), 1, max)
      sigma_int <- update_sigma(X_m, rho, R_int, alpha_int,
                                bet_int)
      princ <- wpca(X_u, q = 15, weighted = TRUE)
      Lam_u_int <- princ$Lam_vec
      W_u_int <- princ$loadings
      hZ <- princ$PCs
      n_c <- colSums(R_int)
      Mu_u_int <- t(sapply(1:K, function(k) 1/n_c[k] * colSums(R_int[,
                                                                     k] * hZ)))
      q <- ncol(Mu_u_int)
      Sgm_u_int <- init_Sgm(R_int, hZ, matrix(0, ncol = q,
                                              nrow = q), Mu_u_int, FALSE)
    }
    else {
      print("You specify the wrong way of generating initial value")
      break
    }
    beta_int <- 1.5
    fit <- icmem(X_m, X_u, Adj_sp, rho, lfc, y_int, Pi_u_int *
                   0, beta_int, beta_grid, alpha_int, bet_int, sigma_int, Mu_u_int,
                 W_u_int, Sgm_u_int, Lam_u_int, maxIter, maxIter_ICM,
                 epsLogLik, verbose, !error_heter, Sigma_diag)

    fit$type <- colnames(rho)[fit$type]
    return(fit)
}

#' Make annotations for SRT data by fitting a SpatialAnno model
#'
#' @param SRTProj an object named "SRTProject". The object SRTProject is created by \link{CreateSRTProject}.
#' @param markers a named string vector, the marker genes used for annotation.
#' @param q An optional integer, specify the number of low-dimensional embeddings to extract in DR-SC model.
#' @param ... other parameters pass to \link{model_set_SpatialAnno}.
#' \href{https://www.biorxiv.org/content/10.1101/2023.02.08.527590v1}{Shi, X., Yang, Y., Ma, X., Zhou, Y., Guo, Z., Wang, C., & Liu, J. (2023). Probabilistic cell/domain-type assignment of spatial transcriptomics data with SpatialAnno. bioRxiv, 2023-02.}
#' @export
#' @importFrom SpatialAnno icmem
#'
Annotation_SpatialAnno <- function(SRTProj, markers, q= 15, ...){


  verbose <- getSRTVerbose()
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')
  if(length(samplenames)>1) stop("Cluster_DRSC: can not applied to multiple data batches!")
  ## Get variable features:
  var.features <- row.names(SRTProj@geneMetaData)[SRTProj@geneMetaData[, 1]]
  use_features <- c(var.features, unlist(markers))
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  XList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[use_features,]
    Matrix::t(data)
  })
  AdjList <- .mylapply(seq_along(samplenames), function(id){
    .getSparseMatrixFromH5file(hfile, groupname = paste0('AdjMat/AdjMat_', id))
  })

  o <- h5closeAll()


  XList <- lapply(XList, scale, scale=TRUE)

  modelSettings <- model_set_SpatialAnno(verbose=verbose, ...)

  resList <- SpatialAnno_structure(XList, markers = markers, AdjList = AdjList, q=q,
                             parameterList = modelSettings)
  ## Write the results to path rds file
  outputPath <- SRTProj@projectMetadata$outputPath
  filename <- SRTProj@projectMetadata$filename
  file_output <- paste0(outputPath, "/", filename)
  if(!dir.exists(file_output)){
    dir.create(file_output)
  }
  save(resList, file= paste0(file_output, "/resSpatialAnno.rds"))





  SRTProj@reductions[["SpatialAnno"]] <- resList$Ez_u
  SRTProj@clusters <- factor(resList$type)
  # SRTProj@models$PRECAST[["R"]] <- reslist$R

  return(SRTProj)
}



