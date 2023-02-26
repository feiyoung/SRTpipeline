
#' Add tSNE embeddings for a SRTProject object
#'
#' @param SRTProj a \code{SRTProject} object.
#' @param n_comp an optional positive integer, specify the number of features to be extracted.
#' @param reduction an optional string, means which dimensional reduction (e.g. PRECAST, PCA) to use for the tSNE, default as the last one in the SRTProject object.
#' @param seed  a non-negative integer, the random seed for reproducibility.
#
#' @return Return a revised SRTProject object by adding tSNE embeddings
#' metadata
#'
#' @export
#'
AddTSNE <- function(SRTProj,  n_comp=2, reduction = NULL, seed=1) {
  UseMethod(generic = 'AddTSNE', object = SRTProj)
}

#' @rdname AddTSNE
#' @method AddTSNE SRTProject
#' @export
AddTSNE.SRTProject <- function(SRTProj, n_comp=2, reduction = NULL, seed=1){

  if(is.null(reduction)){
    reduction <- tail(names(SRTProj@reductions),1)
  }

  set.seed(seed)
  hZ_tsne <- scater::calculateTSNE(t(SRTProj@reductions[[reduction]]), ncomponents=n_comp)
  if(n_comp==3){
    name <- paste0("tSNE3")

  }else if(n_comp==2){
    name <- paste0("tSNE")
  }
  SRTProj@plotEmbeddings[[name]] <- hZ_tsne
  return(SRTProj)
}

#' Add UMAP embeddings for a SRTProject object
#'
#' @param SRTProj a \code{SRTProject} object.
#' @param n_comp an optional positive integer, specify the number of features to be extracted.
#' @param reduction an optional string, means which dimensional reduction (e.g. PRECAST, PCA) to use for the UMAP, default as the last one in the SRTProject object.
#' @param seed  a non-negative integer, the random seed for reproducibility.
#
#' @return Return a revised SRTProject object by adding UMAP embeddings
#'
#' @export
AddUMAP <- function(SRTProj,  n_comp=2, reduction= NULL, seed=1) {
  UseMethod(generic = 'AddUMAP', object = SRTProj)
}

#' @rdname AddUMAP
#' @method AddUMAP SRTProject
#' @export
AddUMAP.SRTProject <- function(SRTProj, n_comp=2, reduction = NULL, seed=1){

  if(is.null(reduction)){
    reduction <- tail(names(SRTProj@reductions),1)
  }

  set.seed(seed)
  hZ_umap <- scater::calculateUMAP(t(SRTProj@reductions[[reduction]]), ncomponents=n_comp)

  if(n_comp==3){
    name <- paste0("UMAP3")

  }else if(n_comp==2){
    name <- paste0("UMAP")
  }


  SRTProj@plotEmbeddings[[name]] <- hZ_umap

  return(SRTProj)
}

.approxPCA <- function(X, q){ ## speed the computation
  require(irlba)
  n <- nrow(X)
  svdX  <- irlba(A =X, nv = q)
  PCs <- svdX$u %*% diag(svdX$d[1:q])
  loadings <- svdX$v
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(PCs = PCs, loadings = loadings, Lam_vec = Lam_vec))
}

#' Add PCA embeddings for a SRTProject object
#'
#' @param obj a \code{SRTProject} object.
#' @param n_comp an optional positive integer, specify the number of features to be extracted.
#' @param Method a string, the method to compute PCA scores. There are three versions to provide: approximated PCA (APCA), PCA, and weighted PCA (WPCA).
#' @param seed  a non-negative integer, the random seed for reproducibility.
#' @return Return a revised SRTProject object by adding PCA embeddings
#' @export
#'
AddPCA <- function(obj, n_comp=15,Method = c('APCA', 'PCA', 'WPCA'), seed=1) UseMethod("AddPCA")

#' @rdname AddPCA
#' @method AddPCA SRTProject
#' @export
AddPCA.SRTProject <- function(obj, n_comp=15, Method = c('APCA', 'PCA', 'WPCA'), seed=1){


  if(!require(DR.SC, quietly = TRUE) && toupper(Method) %in% c("PCA", "WPCA")){
     stop("AddPCA: DR.SC package is required to run PCA or WPCA!")
  }

  Method <- match.arg(Method)
  if(toupper(Method) %in% c("PCA", "WPCA")){
    wpca <- getFromNamespace("wpca", ns ='DR.SC')
  }else if(toupper(Method)!="APCA"){
    stop(paste0("AddPCA: Method =", Method, " is not supported!") )
  }

  ## read log-normalized data
  ## Get variable features:
  var.features <- row.names(obj@geneMetaData)[obj@geneMetaData[, 1]]
  o <- h5closeAll()
  ## Get normalized var.features*spots data
  hfile <- obj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')
  XList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[var.features,]
    Matrix::t(data)
  })
  o <- h5closeAll()


  XList <- lapply(XList, scale, scale=FALSE)

  if(toupper(Method) == "PCA"){
    princ <- wpca(matlist2mat(XList), q=n_comp, weighted = FALSE)
  }else if(toupper(Method) == "WPCA"){
    princ <- wpca(matlist2mat(XList), q=n_comp, weighted = FALSE)
  }else{
    princ <- .approxPCA(matlist2mat(XList), q=n_comp)
  }
  obj@reductions[['PCA']] <- princ$PCs
  rm(princ)

  return(obj)
}
