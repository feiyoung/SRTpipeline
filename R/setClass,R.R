
setClass("SRTProject",
         representation(
           projectMetadata = "SimpleList", # record the project-related info, such as outputPath, projectName
           projectSummary = "SimpleList", # record the project summary data
           batchMetaData = "DataFrame", # record the info. about each data batch
           cellMetaData = "DataFrame",  # record the meta data of each spot
           geneMetaData = "DataFrame", # record the meta data of each gene
           spatialCoords = "DataFrame", # record the spatial coordinates
           reductions = "SimpleList",   # record the dimension reductions
           plotEmbeddings = "SimpleList", # record the embeddings based on DR, such as tSNE, UMAP.
           clusters = "DataFrame"    # record the clusters by clustering methods.
         )
)

setMethod("show", "SRTProject",

          function(object) {
            require(S4Vectors)
            scat <- function(fmt, vals=character(), exdent=2, n = 5, ...){
              vals <- ifelse(nzchar(vals), vals, "''")
              lbls <- paste(S4Vectors:::selectSome(vals, maxToShow = n), collapse=" ")
              txt <- sprintf(fmt, length(vals), lbls)
              cat(strwrap(txt, exdent=exdent, ...), sep="\n")
            }
            #.ArchRLogo(ascii = "Package")
            cat("class:", class(object), "\n")
            cat("outputDirectory:", object@projectMetadata$outputDirectory, "\n")

            o <- tryCatch({
              object@cellColData$Sample
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

            scat("samples(%d): %s\n", rownames(object@sampleColData))
            scat("sampleColData names(%d): %s\n", names(object@sampleColData))
            scat("cellMetaData names(%d): %s\n", names(object@cellMetaData))
            scat("numberOfSpots(%d): %s\n", nrow(object@cellMetaData))
            scat("medianTSS(%d): %s\n", median(object@cellMetaData$TSSEnrichment))
            scat("medianFrags(%d): %s\n", median(object@cellMetaData$nFrags))

          }

)

.requirePackage <- function (x = NULL, load = TRUE, installInfo = NULL, source = NULL)
{
  if (x %in% rownames(installed.packages())) {
    if (load) {
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }
    else {
      return(0)
    }
  }
  else {
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0("install.packages(\"",
                              x, "\")")
      }
      else if (tolower(source) == "bioc") {
        installInfo <- paste0("BiocManager::install(\"",
                              x, "\")")
      }
      else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ",
                  installInfo))
    }
    else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

.safelapply <- function (..., threads = 1, preschedule = FALSE)
{
  if (tolower(.Platform$OS.type) == "windows") {
    threads <- 1
  }
  if (threads > 1) {
    .requirePackage("parallel", source = "cran")
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    errorMsg <- list()
    for (i in seq_along(o)) {
      if (inherits(o[[i]], "try-error")) {
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error",
                                capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x,
                                                           1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ",
                                                            i, " : "), capOut), "\n")
      }
    }
    if (length(errorMsg) != 0) {
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
    }
  }
  else {
    o <- lapply(...)
  }
  o
}

.summarizeArrowContent <- function (ArrowFile = NULL)
{ # How to read data from h5 file
  o <- h5closeAll()
  h5DF <- h5ls(ArrowFile)
  h5DF <- h5DF[-which(h5DF$group == "/"), ]
  groups <- stringr::str_split(h5DF$group, pattern = "/",
                               simplify = TRUE)[, 2]
  groupList <- split(h5DF, groups)
  groupList2 <- lapply(seq_along(groupList), function(x) {
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]),
                   "", groupDFx$group)
    if (all(groupx == "")) {
      groupDFx
    }
    else {
      subDF <- groupDFx[-which(groupx == ""), ]
      split(subDF, stringr::str_split(subDF$group, pattern = "/",
                                      simplify = TRUE)[, 3])
    }
  })
  names(groupList2) <- names(groupList)
  o <- h5closeAll()
  return(groupList2)
}
.getMetadata <- function (ArrowFile = NULL)
{
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- .summarizeArrowContent(ArrowFile)$Metadata
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name ==
                                                        "CellNames"]), ]
  md <- lapply(seq_len(nrow(arrowMD)), function(x) {
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x],
                                              "/", arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)
  md$CellNames <- paste0(sampleName, "#", md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md) == "CellNames")]
  md <- md[, order(colnames(md))]
  o <- h5closeAll()
  return(md)
}

getSRTVerbose <- function ()
{
  SRTVerbose <- options()[["SRTpiline.verbose"]]
  if (!is.logical(SRTVerbose)) {
    options(SRTpiline.verbose = TRUE)
    return(TRUE)
  }
  SRTVerbose
}
.logDiffTime <- function(main = "", t1 = NULL, verbose = TRUE, addHeader = FALSE,
          t2 = Sys.time(), units = "mins", header = "###########",
          tail = "elapsed.", precision = 3)
{

  # main = ""; t1 = NULL; verbose = TRUE; addHeader = FALSE;
  # t2 = Sys.time(); units = "mins"; header = "###########";
  # tail = "elapsed."; precision = 3
  if (verbose) {
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),
                      precision))
      if (addHeader) {
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s",
                       header, Sys.time(), main, dt, units, tail,
                       header)
      }
      else {
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(),
                       main, dt, units, tail)
      }
      if (getSRTVerbose())
        message(msg)
    }, error = function(x) {
      if (getSRTVerbose())
        message("Time Error : ", x)
    })
  }

  return(invisible(0))
}

writeSparseMat2h5 <- function(x){

}

CreateSRTProject <- function(seuList,projectName="SRTProject", outputDirectory = projectName){
  ## projectName="SRTProject"; outputDirectory = projectName

  seuList <- PRECAST:::gendata_seulist()
  if(grepl(" ", outputDirectory)){
    stop("outputDirectory cannot have a space in the path! Path : ", outputDirectory)
  }
  dir.create(outputDirectory,showWarnings=FALSE)
  if(grepl(" ", normalizePath(outputDirectory))){
    stop("outputDirectory cannot have a space in the full path! Full path : ", normalizePath(outputDirectory))
  }
  sampleDirectory <- file.path(normalizePath(outputDirectory), "SRTFiles")
  dir.create(sampleDirectory,showWarnings=FALSE)



  ## Write count information to h5file
  prefix <- '****'
  tstart <- Sys.time()
  .logDiffTime(sprintf("%s Tabix Bed To Temporary File",
                       prefix), t1 = tstart, verbose = verbose)
  require(rhdf5)
  o <- h5closeAll()
  tmpFile <- file.path(normalizePath(outputDirectory), "SRTFiles", "tmpFile")
  o <- h5createFile(tmpFile)
  o <- h5createGroup(tmpFile, paste0("counts"))
  o <- h5createGroup(tmpFile, paste0("logNormalizedData"))
  o <- h5createGroup(tmpFile, paste0("Metadata"))
  o <- h5write(obj = "Arrow", file = tmpFile, name = "Class")
  o <- h5write(obj = paste0(packageVersion("SRTpipeline")),
               file = tmpFile, name = "SRTpipelineVersion")
  cellMetaData <- seuList[[1]]@meta.data
  o <- h5write(obj = 'cellMetaData', file = tmpFile, name = "Metadata/cellMetaData")
  o <- h5write(obj = cellMetaData, file = tmpFile, name = "Metadata/cellMetaData2")
  ## Write count matrix
  count <- seuList[[1]][["RNA"]]@counts
  o <- h5write(obj = count, file = tmpFile, name = "counts/count")
  h5ls(tmpFile)


  #Sample Information
  sampleColData <- DataFrame(row.names = sampleNames, ArrowFiles = ArrowFiles)
  sampleMetadata <- SimpleList(lapply(sampleNames, function(x) SimpleList()))
  names(sampleMetadata) <- sampleNames

  #Cell Information
  message("Getting Cell Metadata...")
  metadataList <- .safelapply(seq_along(ArrowFiles), function(x){
    if(getArchRVerbose()) message(x, " ", appendLF = FALSE)
    .getMetadata(ArrowFiles[x])
  }, threads = threads)
  message("")
  message("Merging Cell Metadata...")
  allCols <- unique(c("Sample",rev(sort(unique(unlist(lapply(metadataList,colnames)))))))
  cellColData <- lapply(seq_along(metadataList), function(x){
    mdx <- metadataList[[x]]
    idx <- which(allCols %ni% colnames(mdx))
    if(length(idx) > 0){
      for(i in seq_along(idx)){
        mdx[,allCols[idx]] <- NA
      }
    }
    mdx[, allCols, drop = FALSE]
  }) %>% Reduce("rbind", .) %>% DataFrame


  message("Initializing SRTProject...")
  Proj <- new("SRTProject",
               projectMetadata = SimpleList(outputDirectory = normalizePath(outputDirectory)),
               projectSummary = SimpleList(),
               batchColData = sampleColData,
               sampleMetadata = sampleMetadata,
               cellMetaData = cellColData,
               geneMetaData = DataFrame(),
               spatialCoords = DataFrame(),
               reducedDims = SimpleList(),
               embeddings = SimpleList()
  )

  Proj <- addProjectSummary(AProj, name = "DateOfCreation", summary = c("Date" = Sys.time()))

  Proj
}
