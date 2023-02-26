
#' SRTProject Class to store SRT datasets
#' @description
#' SRTProject object is a class with a number of slots designed for storing the SRT datasets and the analysis results.
#'  Key slots to access are listed below.
#' @export
#' @param projectMetadata a SimpleList object to record the project-related info, such as outputPath, projectName
#' @param projectSummary a SimpleList object to record the project summary data, such as how many spots for each data batches.
#' @param sampleMetaData a DataFrame object to record the info. about each data batch, such as gender, platform and so on.
#' @param cellMetaData a DataFrame object to record the meta data of each spot.
#' @param geneMetaData a SimpleList object to record the meta data of genes, such as SVGs data.frame, HVGs data.frame.
#' @param spatialCoords a DataFrame object to record the spatial coordinates.
#' @param reductions a SimpleList to record the inferred dimension reductions (DRs) by dimension reduction methods.
#' @param plotEmbeddings a SimpleList to record the embeddings based on DRs, such as tSNE, UMAP.
#' @param clusters a factor object to record the inferred clusters by clustering methods.
#' @param models a SimpleList to record the model parameters for specified model.
#'
#' @return Returns a UMAP representation.
#'
#' @importFrom scater calculateUMAP
#' @seealso \code{\link{CreateSRTProject}}
#'


setClass("SRTProject",
         representation(
           projectMetadata = "SimpleList", # record the project-related info, such as outputPath, projectName
           projectSummary = "SimpleList", # record the project summary data
           sampleMetaData = "DataFrame", # record the info. about each data batch
           cellMetaData = "DataFrame",  # record the meta data of each spot
           geneMetaData = "SimpleList", # record the meta data of genes, such as SVGs data.frame, HVGs data.frame;
           spatialCoords = "DataFrame", # record the spatial coordinates
           reductions = "SimpleList",   # record the dimension reductions
           plotEmbeddings = "SimpleList", # record the embeddings based on DR, such as tSNE, UMAP.
           clusters = "factor",    # record the clusters by clustering methods.
           models = "SimpleList"      # Record the model parameters for specified model
         )
)
.SRTpiplneLogo <- function (messageLogo = TRUE)
{

  logo <-c( " .oooooo..o ooooooooo.   ooooooooooooo             o8o                       oooo   o8o             \n",
            "d8P'    `Y8 `888   `Y88. 8'   888   `8             `''                       `888   `''                      \n",
            "Y88bo.       888   .d88'      888      oo.ooooo.  oooo  oo.ooooo.   .ooooo.   888  oooo  ooo. .oo.    .ooooo.\n",
            "`'Y8888o.    888ooo88P'       888       888' `88b `888   888' `88b d88' `88b  888  `888  `888P'Y88b  d88' `88b\n",
            "``Y88b888`   88b.             888       888   888  888   888   888 888ooo888  888   888   888   888  888ooo888    \n",
            "oo     .d8P  888 `88b.        888       888   888  888   888   888 888    .o  888   888   888   888  888    .o\n",
            "8''88888P'  o888o  o888o     o888o      888bod8P' o888o  888bod8P' `Y8bod8P' o888o o888o o888o o888o `Y8bod8P' \n",
            "                                        888              888                                                   \n",
            "                                       o888o            o888o                                                  \n")

  SRTlogo <- paste0(logo)
  if (messageLogo) {
    message(SRTlogo)
  }
  else {
    cat(SRTlogo)
  }
}

setMethod("show", "SRTProject",
          function(object) {
            require(S4Vectors)
            scat <- function(fmt, vals=character(), exdent=2, n = 5, ...){
              vals <- ifelse(nzchar(vals), vals, "''")
              lbls <- paste(S4Vectors:::selectSome(vals, maxToShow = n), collapse=" ")
              txt <- sprintf(fmt, length(vals), lbls)
              cat(strwrap(txt, exdent=exdent, ...), sep="\n")
            }
            .SRTpiplneLogo()
            cat("class:", class(object), "\n")
            cat("outputPath:", object@projectMetadata$outputPath, "\n")
            cat("h5filePath:", object@projectMetadata$h5filePath, "\n")
            o <- tryCatch({
              object@cellMetaData
            }, error = function(x){
              stop(paste0("\nError accessing sample info from ArchRProject.",
                          "\nThis is most likely the issue with saving the ArchRProject as an RDS",
                          "\nand not with save/loadArchRProject. This bug has mostly been attributed",
                          "\nto bioconductors DataFrame saving cross-compatability. We added a fix to this.",
                          "\nPlease Try:",
                          "\n\trecoverArchRProject(ArchRProj)",
                          "\n\nIf that does not work please report to Github: https://github.com/GreenleafLab/ArchR/issues"
              ))
            })
            cat("---------Datasets basic information-----------------\n")
            scat("samples(%d): %s\n", rownames(object@sampleMetaData))
            scat("sampleColData names(%d): %s\n", names(object@sampleMetaData))
            scat("cellMetaData names(%d): %s\n", names(object@cellMetaData))
            scat("numberOfSpots(%d): %s\n", table(object@cellMetaData$batch))

            cat("---------Downstream analyses information-----------------\n")
            if(is.element("isVGs#combine", colnames(object@geneMetaData))){
              cat("Variable features: ", sum(object@geneMetaData[,"isVGs#combine"]), "\n")
            }
            # scat("Prepared models(%d): %s\n", names(object@models))
            scat("Low-dimensional embeddings(%d): %s\n", names(object@reductions))
            if(length(object@clusters)){
              cat("Inferred cluster labels: Yes\n")
            }else{
              cat("Inferred cluster labels: No\n")
            }

            scat("Embedding for plotting(%d): %s\n", names(object@plotEmbeddings))
          }

)

#' Add the summary data for a SRTProject object
#'
#' @param SRTProj a SRTProject object.
#' @param name a string, set a name for the new summary data to be added.
#' @param summary any object, the new summary data to be added.
#'
#' @return return a revised SRTProject object.
#' @export
#'
addProjectSummary <- function (SRTProj = NULL, name = NULL, summary = NULL)
{

  pS <- SRTProj@projectSummary
  name <- paste0(length(pS) + 1, "_", name)
  pS <- append(pS, SimpleList(summary))
  names(pS)[length(pS)] <- name
  SRTProj@projectSummary <- pS
  SRTProj
}

createSRTFiles <- function (cntList, coordList, filename, cellMetaDataList = NULL,
                            sampleMetadata = NULL, min.spots = 0, min.genes = 0, filter.plot= FALSE,
                            prefix="****",
                            force=FALSE, verbose=TRUE){

  require(S4Vectors)
  # require(dplyr)
  require(rhdf5)

  if(!inherits(cntList, "SimpleList") && !inherits(cntList, "list")) stop("cntList must be a SimpleList or list!")
  if(length(names(cntList))!=length(cntList)) stop("cntList must have a name for each component!")
  if(length(unique(names(cntList))) !=length(cntList))
    stop("The componnets of cntList must have different names!")
  for(r in seq_along(cntList)){
    if(!all(colnames(cntList[[r]]) == row.names(coordList[[r]]))){
      stop("cntList colnames must be equal to row names of coordList")
    }
  }

  if(!is.null(metadataList)){
    for(r in seq_along(cntList)){
      if(!all(colnames(cntList[[r]]) == row.names(metadataList[[r]]))){
        stop("cntList colnames must be equal to row names of metadataList")
      }
    }
  }

  # if (grepl(".", filename)) {
  #   stop("filename cannot have a dot in the name! Name : ",
  #        filename)
  # }

  if(verbose)
    message(prefix,"Remove the unshared genes for all batches...")


  ## Start record time
  tstart <- Sys.time()

  samplenames <- names(cntList)
  #### Merge count matrices into a matrix
  ## filter genes
  rowsum_nCount_cntList <- lapply(cntList, function(x) Matrix::rowSums(x>0))
  if(filter.plot){
    p1 <- volin.filter.plot(rowsum_nCount_cntList)
  }
  # For each data batch, each gene has at least min.spots spots reads.
  keep_gene_list <- lapply(rowsum_nCount_cntList, function(x) names(x)[x>= min.spots])
  if(filter.plot){
    p2 <- volin.filter.plot(rowsum_nCount_cntList, keep_gene_list)
  }

  if(verbose){
    if(length(cntList)>1){
      message(prefix, "Filter out ", paste0(sapply(cntList, nrow)-sapply(keep_gene_list, length),
                                            collapse = ', '), " genes for ",
              length(cntList), " data batches, respectively.")
    }else{
      message(prefix, "Filter out ", paste0(sapply(cntList, nrow)-sapply(keep_gene_list, length),
                                            collapse = ', '), " genes for ",
              length(cntList), " this data batch.")
    }
  }

  shared_genes <- Reduce(intersect, keep_gene_list)
  make_unique <- function(x, verbose){
    if(anyDuplicated(x)){
      if(verbose)
        message("Revise the duplicated spots names or gene names")
      x <- make.unique(x)
    }
    return(x)
  }

  if(verbose)
    message(prefix, "To differ the spots across data batch,\n",
            prefix, "revise the name of each spot in each data batch by adding the name of each data batch\n")
  for(r in seq_along(cntList)){
    colnames(cntList[[r]]) <- paste0(names(cntList)[r],"#",colnames(cntList[[r]]))
  }
  ## obtain gene filtered data
  cntList <- lapply(cntList, function(x) x[shared_genes, ])

  ## Filter Spots
  colsum_nCount_cntList <- lapply(cntList, function(x) Matrix::colSums(x>0))
  if(filter.plot){
    boxplot(rowsum_nCount_cntList)
    p1 <- volin.filter.plot(rowsum_nCount_cntList)

  }
  raw_number <- sapply(cntList, ncol)
  keep_spot_list <- lapply(colsum_nCount_cntList, function(x) which(x>= min.genes))

  cntList <- lapply(seq_along(cntList), function(r) cntList[[r]][, keep_spot_list[[r]]])
  if(verbose)
    message(prefix, "Filter out ", paste0(raw_number -sapply(keep_spot_list, length),
                                          collapse = ', '), " spots for ",
            length(cntList), " data batches, respectively")

  # counts_all <- Reduce(cbind2, cntList)
  # row.names(counts_all) <- make_unique(row.names(counts_all), verbose)
  # colnames(counts_all) <- make_unique(colnames(counts_all), verbose)
  for(r in seq_along(cntList)){
    row.names(cntList[[r]]) <- make_unique(row.names(cntList[[r]]), verbose)
    colnames(cntList[[r]]) <- make_unique(colnames(cntList[[r]]), verbose)
  }

  #### Merge spatial coordinates into a matrix
  ## filter spots
  coordList <- lapply(seq_along(cntList), function(r) coordList[[r]][keep_spot_list[[r]], ])
  ## merge
  coord_all <- as.matrix(Reduce(rbind, coordList))
  row.names(coord_all) <- NULL


  #### Merge meta data
  if(!is.null(cellMetaDataList)){
    ## filter spots
    cellMetaDataList <- lapply(seq_along(cntList), function(r) cellMetaDataList[[r]][keep_spot_list[[r]], ])

    ## merge
    meta_data_all <- Reduce(rbind, cellMetaDataList)
    row.names(meta_data_all) <- NULL
  }

  .logDiffTime(sprintf(paste0("%s Data check and preprocessing"),
                       prefix), t1 = tstart, verbose = verbose)




  ####### write to h5 file #######
  ## Restart record time
  tstart <- Sys.time()
  hfile <- paste0(filename,".h5")
  o <- h5closeAll()
  if(!file.exists(hfile)){
    o <- h5createFile(hfile)
  }
  ### write counts to file
  # o <- .writeData2H5file(hfile, groupname='counts', object=counts_all,
  #                        objectclass= "SparseMatrix", force=force)

  ## require to revise cntList writer for each data batch for convenience's analysis.
  ## write countList to file
  if(!is.element("count", h5ls(hfile, recursive=1)$name)){
    o <- h5createGroup(hfile, 'count')
  }else{
    if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
      h5delete(hfile, 'count')
      o <- h5createGroup(hfile, 'count')
    }else{
      stop(paste0('count', " Group already exists! Set force = TRUE to overwrite!"))
    }
  }
  o <- .mysapply(seq_along(cntList), function(id){

    ### write counts to file
    groupname <- paste0("count/count_", id)
    Group <- paste0('count_', id)
    if(!is.element(Group, h5ls(hfile, recursive=2)$name)){# If no Group, create it and write data
      SparseMatWrite(cntList[[id]], groupname = groupname, hfile, write_names=TRUE)
    }else{
      if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
        h5delete(hfile, groupname)
        SparseMatWrite(cntList[[id]], groupname = groupname, hfile, write_names=TRUE)
      }else{
        stop(paste0(groupname, " Group already exists! Set force = TRUE to overwrite!"))
      }
    }
    return(invisible(0))
  })

  ### write spatial coordinates to hfile
  o <- .writeData2H5file(hfile, groupname='SpatialCoords', object=coord_all,
                         force=force)

  ### write cell meta data to hfile

  if(!is.null(cellMetaDataList)){
    # o <- h5write(obj = meta_data_all, file = hfile, name = "cellMetaData")
    o <- .writeData2H5file(hfile, groupname='cellMetaData', object=meta_data_all,
                           force=force)
  }

  ### write sample meta data to hfile
  if(is.null(sampleMetadata)){
    sampleMetadata <- DataFrame(row.names = names(cntList))
  }
  sampleMetadata[['NumOfSpots']] <- sapply(cntList, ncol)
  # o <- h5write(obj = sampleMetadata, file = hfile, name = "sampleMetadata")
  o <- .writeData2H5file(hfile, groupname='sampleMetadata', object=sampleMetadata,
                         force=force)
  ## write samplenames to hfile
  o <- .writeData2H5file(hfile, groupname='samplenames', object=samplenames,
                         force=force)

  o <- h5closeAll()
  SRTFile <- hfile
  .logDiffTime(sprintf(paste0("%s Write count and meta information into ", SRTFile, " file"),
                       prefix), t1 = tstart, verbose = verbose)

  return(list(SRTFile=SRTFile, shared_genes=shared_genes))
}

#' Determine whether output the information for all functions in SRTpipeline package.
#'
#' @param verbose a data frame
#' @export
#' @examples
#' setSRTverbose(verbose=FALSE)
#'
#'
setSRTverbose <- function(verbose=TRUE){
  options(SRTpipeline.verbose = verbose)
  return(invisible(0))
}

#' Get the output status
#'
#' @export
#' @examples
#' getSRTVerbose()
#'
#'
getSRTVerbose <- function (){
  SRTVerbose <- options()[["SRTpipeline.verbose"]]
  if (!is.logical(SRTVerbose)) {
    options(SRTpipeline.verbose = TRUE)
    return(TRUE)
  }
  SRTVerbose
}

#' Set the Prefix when output informtion
#'
#' @param SRTprefix a string, default as \code{"******"}.
#' @export
#' @examples
#' setSRToutputPrefix()
#'
#'
setSRToutputPrefix <-  function(SRTprefix = "******"){
  options(SRTprefix = SRTprefix)
  return(invisible(0))
}

#' Get the Prefix of output information
#'
#' @export
#' @examples
#' getSRToutputPrefix()
#'
#'
getSRToutputPrefix <- function (){
  prefix <- options()[["SRTprefix"]]
  if (!is.character(prefix)) {
    prefix <-  "******"
    options(SRTprefix = prefix)
    return(prefix)
  }
  prefix
}

#' Create a SRTProject object
#' @description Create a SRTProject object based on at least count matrix list and spatial coordinates list.
#' @param cntList a list with each component as a count matrix, has names for all components.
#' @param coordList a list with spatial coordinates matrix has the same names as that of cntList.
#' @param projectName a string, the name of this project, default as "SRTProject".
#' @param cellMetaDataList a list or NULL, provide the cell meta data for each data batch.
#' @param sampleMetadata a dataframe object with the number of rows same as the length of cntList, provide the sample meta data.
#' @param min.spots a non-negative integer, the least number of nonzero spots for a retained gene used to filter spots in the QC.
#' @param min.genes a non-negative integer, the least number of nonzero genes for a retained spot used to filter spots in the QC.
#' @param filter.plot a logical value, whether plot the QC plots in filtering step.
#' @param prefix a string, specify the prefix for the output information.
#' @param force a logical value, whether delete then rewrite the data to h5 file if the project exists.
#' @param verbose a logical value, whether output the information when creating a SRTProject object.
#' @export
#'
CreateSRTProject <- function(cntList, coordList, projectName="SRTProject", cellMetaDataList = NULL,
                             sampleMetadata = NULL, min.spots = 0, min.genes = 0, filter.plot= FALSE,
                             prefix="****",
                             force=FALSE, verbose=TRUE){
  ## projectName="SRTProject";
  require(rhdf5)
  require(S4Vectors)
  suppressPackageStartupMessages(require(ff)) ## to move SRT h5 file to other path
  nsamples <- length(cntList)
  file_genenameList <- createSRTFiles(cntList, coordList, filename=projectName, cellMetaDataList,
                                sampleMetadata, min.spots, min.genes, filter.plot,
                                prefix, force, verbose)
  SRTFile <- file_genenameList$SRTFile
  geneNames <- file_genenameList$shared_genes
  if(ncol(cellMetaDataList[[1]]) < 2) cm_colnames <- colnames(cellMetaDataList[[1]])

  ## remove the unused object to save memory
  rm(cntList, coordList, cellMetaDataList)
  tstart <- Sys.time()
  if (grepl(" ", projectName)) {
    stop("projectName cannot have a space in the project name! Name : ",
         projectName)
  }
  dir.create(projectName, showWarnings = FALSE)
  outputPath <- normalizePath(projectName)

  # Project meta data
  h5filePath_old <- normalizePath(SRTFile)
  h5filePath_new <- paste0(outputPath,'/' ,SRTFile)
  file.move(h5filePath_old, h5filePath_new)

  #Sample meta Information
  if(verbose)
    message("Getting Sample Metadata...")
  sampleMetadata <- .getDataFromH5file(h5filePath_new, groupname = "sampleMetadata")
  cellNamesList <- lapply(1:nsamples, function(id){
    .getDataFromH5file(h5filePath_new, groupname = paste0("count/count_",id,"/cellNames"))
  })

  #Cell meta Information
  if(verbose)
    message("Getting Cell Metadata...")
  cellMetaData <- .getDataFromH5file(h5filePath_new, groupname = "cellMetaData")
  if(length(dim(cellMetaData)) <2 || is.null(dim(cellMetaData))){
    cellMetaData <- data.frame(cellMetaData)
    colnames(cellMetaData) <- cm_colnames
  }

  samplenames <- .getDataFromH5file(h5filePath_new, groupname = "samplenames")
  cellMetaData$batch <- unlist(lapply(seq_along(samplenames),
                                      function(i) rep(samplenames[i],sampleMetadata$NumOfSpots[i])))
  row.names(cellMetaData) <- unlist(cellNamesList)
  # Spatial coordinates
  if(verbose)
    message("Getting Spatial Coordinates...")
  SpatialCoords <- .getDataFromH5file(h5filePath_new, groupname = "SpatialCoords")

  .logDiffTime(sprintf(paste0("%s Read meta information from ", SRTFile, " file"),
                       prefix), t1 = tstart, verbose = verbose)

  if(verbose)
    message("Initializing SRTProject...")
  row.names(sampleMetadata) <- samplenames
  Proj <- new("SRTProject",
              projectMetadata = SimpleList(projectName=projectName, outputPath = outputPath, h5filePath=h5filePath_new),
              projectSummary = SimpleList(),
              sampleMetaData = DataFrame(sampleMetadata),
              cellMetaData = DataFrame(cellMetaData),
              geneMetaData = DataFrame(row.names = geneNames),
              spatialCoords = DataFrame(SpatialCoords),
              reductions = SimpleList(),
              plotEmbeddings = SimpleList(),
              clusters = factor(),
              models = SimpleList()
  )

  Proj <- addProjectSummary(Proj, name = "DateOfCreation", summary = c("Date" = Sys.time()))

  Proj
}



#' Normalize the expression count data
#'
#' @param SRTProj a SRTProject object.
#' @param normalization.method a string, specify the method to do normalization.
#' @param force a logical value, whether overwrite the data objects in h5file.
#' @param ... other arguments pass to \code{\link{NormalizeData}}.
#' @export
#' @importFrom Seurat CreateSeuratObject NormalizeData
normalizeSRT <- function(SRTProj, normalization.method='LogNormalize', force=FALSE,  ...){


  require(Seurat)


  ## Get the global settings
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()

  hfile <- SRTProj@projectMetadata$h5filePath
  sampleNames <- row.names(SRTProj@sampleMetaData)


  ## write data to file
  o <- h5closeAll()
  group_new <- 'data'
  if(!is.element(group_new, h5ls(hfile, recursive=1)$name)){
    o <- h5createGroup(hfile, group_new)
  }else{
    if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
      h5delete(hfile, group_new)
      o <- h5createGroup(hfile, group_new)
    }else{
      stop(paste0(group_new, " Group already exists! Set force = TRUE to overwrite!"))
    }
  }
  o <- h5closeAll()
  # Create a seurate object for each unique sample_id
  dataFun <- function(id){
    # browser()

    tstart <- Sys.time()

    count <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('count/count_', id))
    ret_seurat <- CreateSeuratObject(counts= count)
    rm(count)
    ret_seurat <- NormalizeData(ret_seurat, normalization.method ==normalization.method, verbose=verbose, assay = "RNA" ,... )
    # ret_seurat <- NormalizeData(ret_seurat, normalization.method =normalization.method, verbose=verbose)
    data <- ret_seurat[['RNA']]@data

    ### write data to file
    groupname <- paste0(group_new, "/", group_new,"_", id)
    Group <- paste0(group_new, "_", id)
    if(!is.element(Group, h5ls(hfile, recursive=2)$name)){# If no Group, create it and write data
      SparseMatWrite(data, groupname = groupname, hfile, write_names=TRUE)
    }else{
      if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
        h5delete(hfile, groupname)
        SparseMatWrite(data, groupname = groupname, hfile, write_names=TRUE)
      }else{
        stop(paste0(groupname, " Group already exists! Set force = TRUE to overwrite!"))
      }
    }
    rm(data)
    .logDiffTime(sprintf(paste0("%s Finish normalization for sample ", id), prefix), t1 = tstart, verbose = verbose)


    return(invisible(0))
  }
  o <- .mysapply(seq_along(sampleNames), dataFun)
  o <- h5closeAll()

  return(SRTProj)
}
.selectIntFeatures <- function(countList, featureList, IntFeatures=2000, verbose=TRUE){
  ## This function is used for selecting common informative features
  if(length(countList) != length(featureList)) stop("The length of suelist and featureList must be equal!")


  if(any(sapply(featureList, length)< IntFeatures))
    stop("Feature list exists number of features less than IntFeatures!")

  geneUnion <- unique(unlist(lapply(featureList, function(x) x[1:IntFeatures])))

  ## ensure each seuobject has the genes in geneUnion
  gene_delete <- unique(unlist(lapply(countList, function(x) setdiff(geneUnion, row.names(x)))))

  geneUnion <- setdiff(geneUnion, gene_delete)
  # Remove zero-variance genes
  genes_zeroVar <- unique(unlist(.mylapply(countList, function(x)
    geneUnion[Matrix::rowSums(x[geneUnion,])==0])))

  if(length(genes_zeroVar)>0 && verbose){
    message("Remove the zero-vairance genes in some data batch!")
  }

  gene_Var <- setdiff(geneUnion, genes_zeroVar)

  # sort by number of datasets that identified this gene as variable features
  nsample <- length(countList)
  numVec <- rep(0, length(gene_Var))
  rankMat <- matrix(NA,length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for(i in 1:length(gene_Var)){
    for(j in 1:nsample){
      if(is.element(gene_Var[i], featureList[[j]])){
        numVec[i] <- numVec[i] +1
        rank1 <- which(featureList[[j]]==gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }

  }

  cutNum <- sort(numVec, decreasing = T)[min(IntFeatures, length(numVec))]
  if(max(numVec)> cutNum){
    genelist1 <- gene_Var[numVec>cutNum]
  }else{
    genelist1 <- NULL
  }
  num_rest_genes <- min(IntFeatures, length(numVec)) - length(genelist1)

  gene2 <- gene_Var[numVec==cutNum]

  ### select top 2000 genes that rank
  rankMat2 <- rankMat[gene2, ]
  rowMedian <- function(xmat, na.rm=TRUE){
    apply(xmat, 1, median, na.rm=na.rm)
  }
  genes1rank <- gene2[order(rowMedian(rankMat2, na.rm=T))[1:num_rest_genes]]
  genelist <- c(genelist1, genes1rank)

  return(genelist)
}

#' Select variable genes
#'
#' @param SRTProj a SRTProject object.
#' @param nfeatures a positive integer, the number of variable features to be chosen.
#' @param type a string, one of "HVGs" and "SVGs", represent the type of variable genes to be chosen.
#' @param method a string, the method to use, including 'vst' if \code{type='HVGs'}, "SPARK-X" if \code{type='SVGs'}.
#' @param use_custom_features a string vector, a group of custom genes to use. If this argument is not \code{NULL}, it will use this group of custom genes as variable genes.
#' @param exclude_features a string vector, a group of genes to be excluded when finding variable genes.
#' @return a revised SRTProject object.
#' @export
selectVariableFeatures <- function(SRTProj, nfeatures = 2000, type = c("HVGs", "SVGs"), method='vst',
                                   use_custom_features=NULL, exclude_features=NULL){

  require(Seurat)
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()

  type <- match.arg(type)
  hfile <- SRTProj@projectMetadata$h5filePath

  geneMetaData <- SRTProj@geneMetaData

  uniq_sample_id <- row.names(SRTProj@sampleMetaData)
  nTotalSpots <- nrow(SRTProj@cellMetaData)
  geneMetaData[[paste0("isVGs#combine")]] <- rep(NA, nrow(geneMetaData))
  if(!is.null(use_custom_features)){
    geneMetaData[[paste0("isVGs#combine")]] <- rep(FALSE, nrow(geneMetaData))
    geneMetaData[use_custom_features,"isVGs#combine"] <- TRUE
    SRTProj@geneMetaData <- geneMetaData
    return(SRTProj)
  }


  if(type=="HVGs"){
    for(id in uniq_sample_id){
      geneMetaData[[paste0("isHVGs#", id)]] <- rep(NA, nrow(geneMetaData))
      geneMetaData[[paste0("rankHVGs#", id)]] <- rep(NA, nrow(geneMetaData))
    }
  }else if(type=="SVGs"){
    for(id in uniq_sample_id){
      geneMetaData[[paste0("isSVGs#", id)]] <- rep(NA, nrow(geneMetaData))
      geneMetaData[[paste0("rankSVGs#", id)]] <- rep(NA, nrow(geneMetaData))
    }
  }



  # Create a seurate object for each unique sample_id
  tstart <- Sys.time()
  countList <- list()
  featureList <- list()
  for(id in seq_along(uniq_sample_id)){
    # browser()
    if(verbose){
      message("Find variable genes for sample ", uniq_sample_id[id])
    }
    meta_data <- SRTProj@spatialCoords[SRTProj@cellMetaData$batch == uniq_sample_id[id], 1:2]
    colnames(meta_data) <- c("row", "col")
    count <- .getSparseMatrixFromH5file(hfile, groupname = paste0('count/count_', id))
    row.names(meta_data) <- colnames(count)
    if(!is.null(exclude_features)){  ## Exclude these features when find the highly variable genes
      idx_remove <- which(row.names(count) %in% exclude_features)
      count <- count[-idx_remove,]
    }
    ret_seurat <- CreateSeuratObject(counts= count,
                                     meta.data = as.data.frame(meta_data))
    countList[[id]] <- count
    if(type=='HVGs'){
      #ret_seurat <- FindVariableFeatures(ret_seurat, nfeatures=nfeatures,selection.method=method, verbose=verbose, ...)

      ret_seurat <- FindVariableFeatures(ret_seurat, nfeatures=nfeatures,
                                         selection.method=method, verbose=verbose)
      assay <- DefaultAssay(ret_seurat)
      geneMetaData[, paste0("isHVGs#", uniq_sample_id[id])] <- FALSE
      geneMetaData[ret_seurat[[assay]]@var.features, paste0("isHVGs#", uniq_sample_id[id])] <- TRUE
    }else if(type=='SVGs'){
      if(tolower(method)=='spark-x'){
        # ret_seurat <- DR.SC::FindSVGs(ret_seurat, nfeatures=nfeatures, verbose=verbose, ...)
        ret_seurat <- DR.SC::FindSVGs(ret_seurat, nfeatures=nfeatures, verbose=verbose)
        assay <- DefaultAssay(ret_seurat)
        geneMetaData[, paste0("isSVGs#", uniq_sample_id[id])] <- FALSE
        geneMetaData[ret_seurat[[assay]]@var.features, paste0("isSVGs#", uniq_sample_id[id])] <- TRUE
        geneMetaData[, paste0("rankSVGs#", uniq_sample_id[id])] <- ret_seurat[[assay]]@meta.features$order.SVGs
      }
    }else{
      stop("selectVariableFeatures: check the argument: Method! It only supports 'SPARK-X' and 'vst' to select genes now. You can provide self-selected genes using customGenelist argument.")
    }
    var.features <- ret_seurat[[assay]]@var.features
    featureList[[id]] <- var.features
    rm(ret_seurat, count)
  }
  ## reduce memory
  countList <- .mylapply(countList, function(x) x[unique(unlist(featureList)),])

  .logDiffTime(sprintf(paste0("%s Select ", type,  "(variable features) for each data batch using  ",method ,
                              " method"), prefix), t1 = tstart, verbose = verbose)

  tstart <- Sys.time()
  if(length(uniq_sample_id)>1){
    genes_combined <- .selectIntFeatures(countList, featureList, nfeatures, verbose=verbose)
    geneMetaData[, "isVGs#combine"] <- FALSE
    geneMetaData[genes_combined, "isVGs#combine"] <- TRUE

  }else{
    geneMetaData[, "isVGs#combine"] <- geneMetaData[,2]
  }
  SRTProj@geneMetaData <- geneMetaData

  .logDiffTime(sprintf(paste0("%s Prioritize ", type,  " based on the number of times they were selected in all data batches"), prefix), t1 = tstart, verbose = verbose)


  return(SRTProj)
}

#' Get the top variables genes
#'
#' @param SRTProj a SRTProject object aftering running \code{\link{selectVariableFeatures}}.
#' @param ntop a positive integer, the number of top variable genes.
#' @param batch a string, number or NULL, specify the batch name to obtain the variable genes, default as \code{NULL}.
#' @param type a string, one of "HVGs" and "SVGs", represent the type of variable genes to be chosen.
#' @return return a string vector, the top variable genes.
#' @export

topVGs <- function(SRTProj, ntop=5, batch=NULL, type=c("HVGs", "SVGs")){

  if(nrow(SRTProj@geneMetaData)==0) stop("topVG: please run `selectVariableFeatures` first, then run topVG!")

  type <- match.arg(type)
  samplenames <- row.names(SRTProj@sampleMetaData)
  if(is.null(batch)){
    batch <- samplenames[1]
  }
  col_rank <- paste0("rank", type, "#", batch)
  col_vg <- paste0("is", type, "#", batch)
  VGs <- row.names(SRTProj@geneMetaData)[SRTProj@geneMetaData[[col_vg]]]
  order_features <- SRTProj@geneMetaData[[col_rank]][SRTProj@geneMetaData[[col_vg]]]


  idx <- order(order_features)[1:ntop]
  return(VGs[idx])

}


#' Add adjacency matrix list for a SRTProject object
#'
#' @param SRTProj a SRTProject object
#' @param platform a string, specify the data collection platform, one of "Visium", "ST" and "Other", default as "Visium". The platform helps to calculate the adjacency matrix by defining the neighborhoods when type="fixed_distance" is chosen.
#' @param type an optional string, specify which type of neighbors' definition. Here we provide two definition: one is "fixed_distance", the other is "fixed_number".
#' @param force a logical value, whether delete and rewrite the adjacency matrix list to h5file if there exists the adjacency matrix.
#' @param ... other arguments pass to \code{\link{getAdj_auto}}.
#' @return  return a SRTProject.
#' @export
#'
AddAdj <- function(SRTProj,platform=c("Visium", "ST", "Other"),  type="fixed_distance",  force=FALSE,...){

  if(!inherits(SRTProj, "SRTProject"))
    stop("AddAdj: Check the argument: SRTProj!  SRTProj must be a SRTProject object.")

  ## calculate Adj and write them to H5file
  verbose <- getSRTVerbose()
  prefix <- getSRToutputPrefix()

  platform <- match.arg(platform)
  tstart <- Sys.time()
  uniq_sample_id <- unique(SRTProj@cellMetaData$batch)
  # Create a seurate object for each unique sample_id
  calAdj <- function(id){
    # browser()
    pos <-  as.matrix(SRTProj@spatialCoords[SRTProj@cellMetaData$batch == id, 1:2])
    if(tolower(type)=='fixed_distance'){
      if(tolower(platform) %in% c("st", "visium")){
        Adj <- PRECAST::getAdj_reg(as.matrix(pos), platform=platform)
      }else{
        Adj <- DR.SC::getAdj_auto(pos, ...)
      }
    }else if (tolower(type) == "fixed_number") {
      Adj <- PRECAST::getAdj_fixedNumber(pos, ...)
    } else {
      stop("AddAdj: Unsupported adjacency  type \"", type, "\".")
    }

    return(Adj)
  }
  o <- h5closeAll()
  hfile <- SRTProj@projectMetadata$h5filePath
  if(!is.element("AdjMat", h5ls(hfile, recursive=1)$name)){
    o <- h5createGroup(hfile, 'AdjMat')
  }else{
    if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
      h5delete(hfile, 'AdjMat')
      o <- h5createGroup(hfile, 'AdjMat')
    }else{
      stop(paste0('AdjMat', " Group already exists! Set force = TRUE to overwrite!"))
    }
  }

  o <- .mysapply(seq_along(uniq_sample_id), function(id){

    Adj <- calAdj(uniq_sample_id[id])
    ### write counts to file
    datasets <- c('indices' = 'i', 'indptr' = 'p', 'Dim'='Dim', 'data' = 'x')
    groupname <- paste0("AdjMat/AdjMat_", id)
    Group <- paste0('AdjMat_', id)
    if(!is.element(Group, h5ls(hfile, recursive=2)$name)){# If no Group, create it and write data
      SparseMatWrite(Adj, groupname = groupname, hfile, write_names=FALSE)
    }else{
      if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
        h5delete(hfile, groupname)
        SparseMatWrite(Adj, groupname = groupname, hfile, write_names=FALSE)
      }else{
        stop(paste0(groupname, " Group already exists! Set force = TRUE to overwrite!"))
      }
    }
    return(invisible(0))
  })
  o <- h5closeAll()
  .logDiffTime(sprintf(paste0("%s Write adjacency matrices to ",SRTProj@projectMetadata$projectName,
                              ".h5 file"), prefix), t1 = tstart, verbose = verbose)



  return(SRTProj)
}



#' Add inferred clusters to the \code{cellMetaData}
#'
#' @param SRTProj a SRTProject object after fitting clustering method.
#' @param name a string, the name created in colData to put the cluster labels.
#' @return  return a revised SRTProject.
#' @export
#'
AddClusters2MetaData <- function(SRTProj, name='Clusters'){

  if(!inherits(SRTProj, "SRTProject"))
    stop("AddClusters2MetaData: Check the argument: SRTProj!  SRTProj must be a SRTProject object.")

  if(legnth(SRTProj@clusters)==0) stop("AddClusters2MetaData: there is no inferred clusters in the SRTProject object!")
  SRTProj@cellMetaData[[name]] <- SRTProj@clusters

  return(SRTProj)
}
