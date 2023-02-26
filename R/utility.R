## read a matrix from h5 file
.MatrixRead <- function(hfile, groupname){

  # groupname <- 'Integrated_PRECAST'
  o <- h5closeAll()
  h5f = H5Fopen(hfile)
  groupName2 <- h5ls(h5f&groupname, recursive = 1)$name
  datasets <- c("x", "rownames", "colnames")
  if(all(datasets %in% groupName2)){
    X_id <- h5f&paste0(groupname, "/x")
    X <- X_id[]
    rowname_id <- h5f&paste0(groupname, "/rownames")
    colname_id <- h5f&paste0(groupname, "/colnames")
    row.names(X) <- rowname_id[]
    colnames(X) <- colname_id[]
  }else{
    X_id <- h5f&paste0(groupname, "/x")
    X <- X_id[]
  }
  o <- h5closeAll()
  return(X)
}

## Write a matrix to h5 file
.MatrixWrite <- function(x, groupname, hfile, force=FALSE, write_names=TRUE){
  require(hdf5r)
  o <- h5closeAll()

  if(is.element(groupname, h5ls(hfile, recursive=1)$name) && force){# If no Group, create it and write data
    h5delete(hfile, groupname)
  }
  o <- h5createGroup(hfile, groupname)
  datasets <- c("x", "rownames", "colnames")
  if(!write_names) datasets <- datasets[1]
  for(i in seq_along(datasets)){
    tmp_name <- paste0(groupname,"/" ,datasets[i])
    if(datasets[i]=='x'){
      o <- h5write(obj = x, file = hfile, name = tmp_name)
    }else if(datasets[i] == 'rownames'){
      o <- h5write(obj = row.names(x), file = hfile, name = tmp_name)
    }else{
      o <- h5write(obj = colnames(x), file = hfile, name = tmp_name)
    }

  }
  o <- h5closeAll()
  return(invisible(x = NULL))
}

SparseMatWrite <- function(x, groupname, hfile, write_names=TRUE) {


  # name = "Matrix"; x <- seuList[[1]][["RNA"]]@counts
  require(hdf5r)
  o <- h5closeAll()
  if(!file.exists(hfile)){
    o <- h5createFile(hfile)
  }
  o <- h5createGroup(hfile, groupname)
  datasets <- c('indices' = 'i', 'indptr' = 'p', "Dim"="Dim", 'data' = 'x', 'geneNames'='Dimnames',
                'cellNames'='Dimnames')
  if(!write_names){
    datasets <- datasets[-c(5:6)]
  }
  for (i in seq_along(along.with = datasets)) { # start to write datasets
    ds.i <- slot(object = x, name = datasets[i])
    tmp_name <- paste0(groupname,"/" ,names(datasets)[i])
    if(names(datasets)[i] == "geneNames"){
      geneNames <- ds.i[[1]]
      o <- h5write(obj = geneNames, file = hfile, name = tmp_name)
    }else if(names(datasets)[i] == "cellNames"){
      o <- h5write(obj = ds.i[[2]], file = hfile, name = tmp_name)
    }else{
      if(groupname=='counts'){
        o <- h5write(obj = as.integer(ds.i), file = hfile, name = tmp_name)
      }else{
        o <- h5write(obj = ds.i, file = hfile, name = tmp_name)
      }

    }

  }
  o <- h5closeAll()

  return(invisible(x = NULL))
}


.getDataFromH5file <- function(hfile, groupname){
  require(rhdf5)
  o <- h5closeAll()
  groupname_hirach <- strsplit(groupname, split='/', fixed=T)[[1]]
  flag_hirach <- sapply(seq_along(groupname_hirach), function(i){
    groupName1 <- h5ls(hfile, recursive=i)$name
    is.element(groupname_hirach[i], groupName1)

  } )
  if(!all(flag_hirach)) stop(paste0("groupname ", groupname, " is not found in the ", hfile, " H5file!"))


  h5f = H5Fopen(hfile)
  dat <- h5f&groupname
  dat_true <- dat[]
  o <- h5closeAll()

  return(dat_true)
}

.getSparseMatrixFromH5file <- function(hfile, groupname='counts'){

  require(rhdf5)
  require(Matrix)

  target <- groupname
  groupname1 <- groupname
  n <- 1
  if(grepl("/", groupname)){
    groupvec <- strsplit(groupname, split='/', fixed=T)[[1]]
    n <- length(groupvec)
    groupname <- groupvec[n]
    groupname1 <- groupvec[1]
  }
  o <- h5closeAll()
  groupName1 <- h5ls(hfile, recursive=n)$name

  if(!is.element(groupname,groupName1)){
    stop(paste0("hfile does not have a group named ", groupname) )
  }
  #h5read(hfile, name="counts/indices")
  h5f = H5Fopen(hfile)
  groupName2 <- h5ls(h5f&target, recursive = 1)$name
  indices <- h5f&paste0(target,"/indices")
  indptr <- h5f&paste0(target,"/indptr")
  dat <- h5f&paste0(target,"/data")
  Dim <- h5f&paste0(target,"/Dim")
  if(groupname1 %in% c("counts", "AdjMat") ){
    sparse.mat <- sparseMatrix(i = indices[] + 1,
                               p = indptr[], dims=Dim[],
                               x = as.integer(dat[]), repr = "T")
  }else{
    sparse.mat <- sparseMatrix(i = indices[] + 1,
                               p = indptr[],dims=Dim[],
                               x =as.double(dat[]), repr = "C")
  }


  if(all(c("geneNames", "cellNames") %in% groupName2)){
    geneNames <- h5f&paste0(target,"/geneNames")
    cellNames <- h5f&paste0(target,"/cellNames")
    row.names(sparse.mat) <- geneNames[]
    colnames(sparse.mat) <- cellNames[]
    # sparse.mat@Dimnames[[1]] <-  geneNames[]
    # sparse.mat@Dimnames[[2]] <- cellNames[]
  }


  o <- h5closeAll()

  return(sparse.mat)
}

.writeData2H5file <- function(hfile, groupname, object, objectclass=NULL,
                              write_names=TRUE, force=FALSE){

  o <- h5closeAll()
  if(is.null(objectclass)) objectclass <- "others"
  if(objectclass == 'SparseMatrix'){
    Group <- groupname
    if(!is.element(Group, h5ls(hfile, recursive=1)$name)){# If no Group, create it and write data
      SparseMatWrite(object, groupname = Group, hfile, write_names)
    }else{
      if(force){ # If Group exists, and enforce overwrite, then re-create it and write data.
        h5delete(hfile, Group)
        SparseMatWrite(object, groupname = Group, hfile, write_names)
      }else{
        stop(paste0(Group, " Group already exists! Set force = TRUE to overwrite!"))
      }
    }
  }else{
    Group <- groupname
    if(!is.element(Group, h5ls(hfile, recursive=1)$name)){ # If no Group, create it and write data
      o <- h5write(obj = object, file = hfile, name = Group)
    }else{
      if(force){
        h5delete(hfile, Group)
        o <- h5write(obj = object, file = hfile, name = Group)
      }else{
        stop(paste0(Group, " Group already exists! Set force = TRUE to overwrite!"))
      }
    }
  }
  o <- h5closeAll()
  return(invisible(0))
}

.mysapply <- function(...){
  if(getSRTVerbose()){
    pbapply::pbsapply(...)
  }else{
    sapply(...)
  }
}
.mylapply <- function(...){
  if(getSRTVerbose()){
    pbapply::pblapply(...)
  }else{
    lapply(...)
  }
}
.myapply <- function(...){
  if(getSRTVerbose()){
    pbapply::pbapply(...)
  }else{
    apply(...)
  }
}

.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
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


.logTime <- function(main='', prefix='*****', versoe=TRUE){

  if(versoe){
    message(paste0(Sys.time()," : ", prefix," ",  main))
  }


}
