

# Integration data --------------------------------------------------------

.get_correct_exp <- function(XList, RfList, M0){

  require(bigmemory)
  require(bigalgebra)
  require(S4Vectors)

  Rf <- matlist2mat(RfList)
  hK <- ncol(Rf)
  rm(RfList)
  X0 <- as.big.matrix(as.matrix(matlist2mat(XList)))
  rm(XList) #

  Xreg <- as.big.matrix(cbind(Rf, M0))
  Xregt <- as.big.matrix(rbind(t(Rf), t(M0)))
  tmpMat <- Xregt %*% Xreg
  beta.hat <- Xregt %*% X0; # ,
  rm(Xreg, Xregt, Rf)

  # message("good!")
  coefmat <- (qr.solve(as.matrix.Vector(tmpMat)) %*% as.matrix.Vector(beta.hat))[-c(1:hK),]
  # message("good2!")
  tmp <- as.big.matrix(M0) %*% coefmat
  rm(M0)
  hX <- X0 - tmp
  hX2 <- as.matrix.Vector(hX)

  return(hX2)
}


#' Get the batch-corrected integrated gene expresssions
#'
#' @param SRTProj a SRTProject object after finishing integration analysis using \code{\link{Integrate_PRECAST}} or \code{\link{Integrate_iSCMEB}}.
#' @param species a string, one of "Human", "Mouse" and "Unknown", default as "Human".
#' @param Method a string, one of "iSC.MEB" and "PRECAST", specify which method is used to remove the batch effects in the gene expression levels.
#' @param custom_housekeep a string vector, a group of user-specified housekeeping genes used for removing batch.
#' @param only.var.features a logical value, whether only return the variable genes.
#' @param add_empty_counts a logical value, whether add a empty count matrix. This argument makes it easy to transfer the returned SpatialExperiment object to a Seurat object using \code{\link{as.Seurat}}.
#' @return return a SpatialExperiment object with integrated normalized gene expression matrix and other information from the SRTProject object.
#' @export
#'
getIntegratedData <- function(SRTProj, species="Human", Method=c("iSC.MEB", "PRECAST"),
                              custom_housekeep=NULL,
                              only.var.features=TRUE, add_empty_counts=FALSE){

  # species="Human"; Method="iSC.MEB"; custom_housekeep=NULL;only.var.features=TRUE
  # suppressMessages(require(Matrix))
  # suppressMessages(require(Seurat))
  require(SpatialExperiment)

  if(!inherits(SRTProj, "SRTProject"))
    stop("getIntegrateSpaData: Check the argument: SRTProj!  SRTProj must be a SRTProject object.")
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')

  if(length(samplenames)==1)
    stop(paste0("getIntegrateSpaData: To run IntegrateSpaData(), the SRTProject must has more than one batches!"))

  if(!tolower(species) %in% c("human", "mouse", "unknown"))
    stop("getIntegrateSpaData: Check the argument: species! it must be one of 'Human', 'Mouse' and 'Unknown'!")

  Method <- match.arg(Method)

  ## Get variable features:
  var.features <- row.names(SRTProj@geneMetaData)[SRTProj@geneMetaData[, 1]]
  ## Integrating all genes
  XList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))
    if(only.var.features){
      Matrix::t(data[var.features,])
    }else{
      Matrix::t(data)
    }

  })


  n_r <- length(XList)
  barcodes_all <- lapply(XList, row.names)
  if(any(duplicated(unlist(barcodes_all)))){
    for(r in 1:n_r){
      row.names(XList[[r]]) <- paste0(row.names(XList[[r]]), r)
    }
  }
  barcodes_all <- lapply(XList, row.names)
  genelist <- tolower(colnames(XList[[1]]))
  lower_species <- tolower(species)
  houseKeep <- switch (lower_species,
                       human = {
                         data(Human_HK_genes, package="PRECAST")
                         intersect((genelist), tolower(Human_HK_genes$Gene))
                       },
                       mouse={
                         data(Mouse_HK_genes, package="PRECAST")
                         intersect((genelist), tolower(Mouse_HK_genes$Gene))
                       },
                       unknown={
                         character()
                       }
  )
  houseKeep <- tolower(c(houseKeep, custom_housekeep))

  if(length(houseKeep) < 5){
    message("Using  fitting results to obtain the batch corrected gene expressions.")
    Mname <- paste0("microEnv.", Method)
    hX <- .get_correct_exp(XList, RfList=SRTProj@models[[Method]]$Rf, M0=SRTProj@reductions[[Mname]])
  }else{
    message("Using both housekeeping gene and fitting results to obtain the batch corrected gene expressions.")
    wpca <- getFromNamespace("wpca", "DR.SC")
    idx <- which((genelist) %in% houseKeep)

    MList <- pbapply::pblapply(XList, function(x) wpca(as.matrix(x[,idx]), q=min(10, length(houseKeep)), F)$PCs)
    M0 <- matlist2mat(MList)
    rm(MList)
    hX <- .get_correct_exp(XList, SRTProj@models[[Method]]$Rf, M0=M0)
  }
  hX <- t(hX)
  row.names(hX) <- colnames(XList[[1]])
  colnames(hX) <- unlist(barcodes_all)

  cellData <- SRTProj@cellMetaData
  cellData$clusters <- SRTProj@clusters
  if(only.var.features){
    rowData <- SRTProj@geneMetaData[var.features,]
  }else{
    rowData <- SRTProj@geneMetaData
  }

  if(!is.null(cellData$sample_id)){ # remove the existing sample_id
    cellData$sample_id <- NULL
  }

  obj <- SpatialExperiment(
      assays = list(logcounts=unname(hX)),
      colData = cellData, rowData = rowData,
      spatialData= SRTProj@spatialCoords,
      reducedDims=c(SRTProj@reductions,SRTProj@plotEmbeddings),
      spatialCoordsNames=colnames(SRTProj@spatialCoords)[1:2],
      sample_id=cellData$batch)

  reducedDim(obj, "Coord") = as.matrix(SRTProj@spatialCoords)

  if(add_empty_counts){ ## Add the empty count matrix so that it can be easily transferred to Seurat object.
    assay(obj, "counts", withDimnames = F) <- sparseMatrix(i=1,j=1, dims=dim(obj))
  }

  return(obj)


}

#' Get the expression matrix from a SRTProject object
#'
#' @param SRTProj a SRTProject object.
#' @return return a a SpatialExperiment object with  gene expression matrix and other information from the SRTProject object.
#' @export
getGeneSpotData <- function(SRTProj){
  ## Obtain the expression data from h5file for downstream analyses
  require(SpatialExperiment)
  if(!inherits(SRTProj, "SRTProject"))
    stop("getIntegrateSpaData: Check the argument: SRTProj!  SRTProj must be a SRTProject object.")
  hfile <- SRTProj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')


  if(length(samplenames)==1){

    count <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('count/count_', 1))
    cellData <- SRTProj@cellMetaData
    cellData$clusters <- SRTProj@clusters
    if("data" %in% h5ls(hfile, recursive = 1)$name){
      data <- .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', 1))
    }
  }else{
    count <- .mylapply(seq_along(samplenames), function(id){
      count <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('count/count_', id))
      return(count)
    })
    count <- Reduce(cbind, count)
    if("data" %in% h5ls(hfile, recursive = 1)$name){
      data <- .mylapply(seq_along(samplenames), function(id){
        data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))
        return(data)
      })
      data <- Reduce(cbind, data)
    }
  }


  cellData <- SRTProj@cellMetaData
  if(length(SRTProj@clusters)>0){
      cellData$clusters <- SRTProj@clusters
  }

  if(length(samplenames)==1){
    reduceD <- c(SRTProj@reductions,SRTProj@plotEmbeddings)
  }else{
    reduceD <- SRTProj@reductions
  }


  geneNames <- row.names(count)
  row.names(count) <- row.names(data) <- paste0("gene", 1:nrow(count))
  gMetaData <- as.data.frame(SRTProj@geneMetaData)[geneNames,]
  if(ncol(SRTProj@geneMetaData)==1){

    gMetaData <- data.frame('y'=gMetaData)
    colnames(gMetaData) <- colnames(SRTProj@geneMetaData)
  }
  row.names(gMetaData) <- row.names(count)
  gMetaData <- DataFrame(gMetaData)


  obj <- SpatialExperiment(
    assays = list(counts = count, logcounts=data), #
    colData = cellData, rowData = gMetaData,
    spatialData= SRTProj@spatialCoords,
    reducedDims=reduceD,
    spatialCoordsNames=colnames(SRTProj@spatialCoords)[1:2],
    sample_id=cellData$batch)
  row.names(obj) <- geneNames

  reducedDim(obj, "Coord") = as.matrix(SRTProj@spatialCoords)

  return(obj)

}


# ## DEG analysis ---------------------------------------------------------


#' Find differential expression genes
#' Find differential expression genes by comparing two clusters.
#' @description
#'
#' @export
#' @param spe spe A SpatialExperiment object with cluster labels for each cell/spot in colData.
#' @param use_cluster a string, the name in colData denoting the cluster labels.
#' @param cluster1 An optional positive integer, specify the number of features to be extracted.
#' @param cluster2 A random seed to be used.
#' @param features a string vector, specify the gene list to do differential gene expression analysis.
#' @param logfc.threshold a real, the threshold of log fold change to determine which genes are returned.
#' @param test.use a string, the test to use.
#' @param ... other arguments pass to \code{\link{FindAllMarkers}}.
#' @return return a dataframe with the DEG analysis results.
#'
#'
#' @importFrom Seurat FindMarkers
#' @seealso None
#'
FindDEGs <- function(spe, use_cluster="clusters", cluster1 = NULL,
                     cluster2 = NULL, features = NULL,
                     logfc.threshold = 0.25,
                     test.use = "wilcox",...){
  require(Seurat)
  require(Matrix)
  require(SpatialExperiment)

  ## Get the global settings
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()

    count <- logcounts(spe)
    seu <- CreateSeuratObject(counts=count)

    Idents(seu) <- factor(colData(spe)[,use_cluster])

    rm(spe)
    tstart <- Sys.time()
    .logTime("Start finding marker genes...", prefix, verbose)
    dat_deg <- FindMarkers(seu,slot = "count", ident.1=cluster1, ident.2=cluster2, features = features,
                           logfc.threshold = logfc.threshold,
                           test.use = test.use,... )
    # dat_deg <- FindMarkers(seu, slot = "count", ident.1='2', ident.2='3')
    .logDiffTime(sprintf(paste0("%s Finish finding marker genes"), prefix), t1 = tstart, verbose = verbose)
    return(DataFrame(dat_deg))

}
#' Find differential expression genes for all clusters by comparing one cluster to the remaining clusters.
#' @description
#'
#' @export
#' @param spe A SpatialExperiment object with cluster labels for each cell/spot in colData.
#' @param use_cluster a string, the name in colData denoting the cluster labels.
#' @param features a string vector, specify the gene list to do differential gene expression analysis.
#' @param logfc.threshold a real, the threshold of log fold change to determine which genes are returned.
#' @param only.pos a logical value, whether only return the genes with positive log fold change.
#' @param test.use a string, the test to use.
#' @param ... other arguments pass to \code{\link{FindAllMarkers}}.
#' @return return a dataframe with the DEG analysis results.
#'
FindAllDEGs <- function(spe, use_cluster='clusters',features = NULL,
                        logfc.threshold = 0.25,only.pos = FALSE,
                        test.use = "wilcox",...){

  require(Seurat)
  require(Matrix)

  ## Get the global settings
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()



  count <- logcounts(spe) # use logcount to conduct the testing
  seu <- CreateSeuratObject(counts=count)

  Idents(seu) <- factor(colData(spe)[,use_cluster])

  tstart <- Sys.time()
  .logTime("Start finding marker genes...", prefix, verbose)
  dat_deg <- FindAllMarkers(seu,slot = "count",  features = features,
                              logfc.threshold = logfc.threshold, only.pos = only.pos,
                              test.use = test.use,... )
  .logDiffTime(sprintf(paste0("%s Finish finding marker genes"), prefix), t1 = tstart, verbose = verbose)
  return(DataFrame(dat_deg))

}

#' Get top pathways from enrichment analysis results
#'
#' @param df a data frame
#' @param ntop a positive integer, how many top terms for each source are returned
#' @param source_set a string vector, specify the sources.
#' @return return a dataframe object including the top pathways.
#' @export
#'
get_top_pathway <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }

  return(df_sub)
}



# Trajectory inference ----------------------------------------------------



#' Implement trajectory inference
#'
#' @description
#' Implement trajectory inference based on low-dimensional embeddings and clusters using Slingshot or other methods.
#' @export
#' @param spe A SpatialExperiment object.
#' @param use_cluster a string, the name of clusters in the colData of the  SpatialExperiment object used for trajectory inference.
#' @param reduction a string, the name of reduction slot used in trajectory inference.
#' @param name a string, the name of inferrend pseudotime, default as "PT".
#' @param traject.method a string, the method used to conduct trajectory inference, default as 'Slingshot'.
#' @param target.clus a string, the cluster to be used as end point.
#' @param start.clus a string, the cluster to be used as the start point.
#' @param seed a positive integer, the random seed.
#' @param return.all.branch a logical value, whether return all inferred braches, default as \code{FALSE} implies only the merged pseudotime is returned.
#' @param ... other arguments pass to other method.
#'
#' @return return a revised SpatialExperiment object. The colData of the SpatialExperiment object adds a column named name (PT).
#'
#' @importFrom slingshot slingshot slingPseudotime
AddTrajectory <- function(spe,  use_cluster='clusters', reduction= NULL,  name = NULL, traject.method='Slingshot',
                           target.clus = NULL, start.clus=NULL, seed = 1,
                          return.all.branch=FALSE, ...){

  suppressPackageStartupMessages(require(slingshot))

  if(!inherits(spe, "SpatialExperiment")){
    stop("AddTrajectory: spe must be a SpatialExperiment object!")
  }
  if(is.null(name)){
    name <- 'PT'
  }

  if(is.null(reduction)){
    reduction <- tail(reducedDimNames(spe),1)
  }
  embeds <- reducedDim(spe, reduction)
  nrow_raw <- nrow(embeds)
  common.pt <- rep(NA, nrow_raw)


  clusters <- colData(spe)[,use_cluster]

  if(!is.null(target.clus)){
    idx_flag <- clusters %in% target.clus
    subclusters <- clusters[idx_flag]
    embeds <- embeds[idx_flag, ]
  }else{
    subclusters <- clusters
  }
  if(tolower(traject.method) == 'slingshot'){
    set.seed(seed)
    sds <- slingshot(embeds, subclusters, start.clus = start.clus, ...)
    pt <- slingPseudotime(sds)
    common.pseudo <- rowMeans(pt, na.rm=TRUE)
    if(!is.null(target.clus)){
      common.pt[idx_flag] <- common.pseudo
    }else{
      common.pt <- common.pseudo
    }
    colData(spe)[[name]] <- common.pt
    if(return.all.branch){
      if(!is.null(target.clus)){
        all.pt <- matrix(NA, nrow=nrow_raw, ncol=ncol(pt))
        colnames(all.pt) <- colnames(pt)
        all.pt[idx_flag,] <- pt
      }else{
        all.pt <- pt
      }
      for(j in 1:ncol(pt)){
        colData(spe)[[paste0(name,".",colnames(all.pt)[j])]] <- all.pt[,j]
      }
    }
  }else{

  }

  return(spe)

}


