# Homepage: https://feiyoung.github.io/SRTpipeline/index.html
# How to custom the website?
# https://www.cnblogs.com/payton/articles/7273790.html
# Add external .html files
# https://stackoverflow.com/questions/65893350/add-pkgdown-html-to-site-articles-generated-by-rmarkdownrender
# require(devtools)
# use_readme_rmd()
# use_news_md()
# use_vignette("test")  #substitute with the name of your package
# use_github_links()
# use_travis()
# use_cran_badge()
# library(rmarkdown)
# render("README.Rmd", md_document())
# library(pkgdown)
# build_site()
# build_home() #
# build_reference()
# build_article("Installation")
# build_article("integration_large_datasets")
# build_article("multibatch_brain12_tutorial")
# https://stackoverflow.com/questions/58744079/use-roxygen2-to-generate-namespace-a-small-example-or-template
# devtools::document()
#setwd(normalizePath('F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\SRTpipeline'))

## Run pkgdown on Linux
# setwd("/nas/LiuWei/Projects/SRTpipeline/SRTpipeline")
# Tools used in preprocessing ---------------------------------------------

firstup <- function(x) {
  ## First letter use upper capital
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


#' Transfer gene names from one fortmat to the other format
#' @description  Transfer gene names from one fortmat to the other format for two species: human and mouse.
#' @param genelist a string vector, the gene list to be transferred.
#' @param now_name a string, the current format of gene names, one of 'ensembl', 'symbol'.
#' @param to_name a string, the  format of gene names to transfer, one of 'ensembl', 'symbol'.
#' @param species a string, the species, one of 'Human' and 'Mouse'.
#' @param Method a string, the method to use, one of 'biomaRt' and 'eg.db'.
#' @return Return a string vector of transferred gene names. The gene names not matched in the database will not change.
#' @export
#' @examples
#' geneNames <- c("ENSG00000171885", "ENSG00000115756")
#' transferGeneNames(geneNames, now_name = "ensembl", to_name="symbol",species="Human", Method='eg.db')
#'
#'
transferGeneNames <- function(genelist, now_name = "ensembl", to_name="symbol",
                              species="Human", Method='biomaRt'){


  if(toupper(now_name) == toupper(to_name)) stop("now_name and to_name must be various!")
  if(! toupper(species) %in% c("HUMAN", "MOUSE")) stop("Check species: the current version only support Human and Mouse!")
  transferredNames <- switch (toupper(species),
    HUMAN = {
      if(tolower(Method)=='eg.db'){
        require(org.Hs.eg.db)
        mapIds(org.Hs.eg.db, keys = genelist,
               keytype = toupper(now_name), column=toupper(to_name))
      }else if(tolower(Method)=='biomart'){
        require(biomaRt)
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
        G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                        values=genelist,mart= mart)

        idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
        G_list <- G_list[idx_in_genelst,]
        idx_dup <- which(duplicated(G_list$ensembl_gene_id))
        G_list <- G_list[-idx_dup,]
        row.names(G_list) <- G_list$ensembl_gene_id
        symbol_list <- G_list[genelist,]$hgnc_symbol

        symbol_list[symbol_list==''] <- NA
        symbol_list

      }else{
        stop("Check Method: the current version only support biomaRt and eg.db!")
      }


    },
    MOUSE= {
      if(tolower(Method)=='eg.db'){
      require(org.Mm.eg.db)
      mapIds(org.Mm.eg.db, keys = genelist,
             keytype = toupper(now_name), column=toupper(to_name))
      }else if(tolower(Method)=='biomart'){
        require(biomaRt)
        mart <- useDataset(" mmusculus_gene_ensembl", useMart("ensembl"))
        G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                        values=genelist,mart= mart)

        idx_in_genelst <- which(G_list$ensembl_gene_id %in% genelist)
        G_list <- G_list[idx_in_genelst,]
        idx_dup <- which(duplicated(G_list$ensembl_gene_id))
        G_list <- G_list[-idx_dup,]
        row.names(G_list) <- G_list$ensembl_gene_id
        symbol_list <- G_list[genelist,]$hgnc_symbol

        symbol_list[symbol_list==''] <- NA
        symbol_list
      }else{
        stop("Check Method: the current version only support biomaRt and eg.db!")
      }

    }
  )

  if(toupper(to_name) == 'SYMBOL'){
    if(toupper(species) == "MOUSE"){
      transferredNames <- firstup(transferredNames)
    }else{
      transferredNames <- toupper(transferredNames)
    }
  }


  flag_na <- is.na(transferredNames)
  if(any(flag_na))
    transferredNames[flag_na] <- genelist[flag_na]

  return(transferredNames)
}




#' Filter genes for a Seurat object
#' @description Filter genes for a single-cell RNA sequencing or sptially-resolved transcriptomics data in the form of a Seurat object.
#' @param seu a Seurat object.
#' @param min_spots a non-negative integer, the least number of nonzero spots for a retained gene.
#' @param assay a string, the assay to perform filtering.
#' @return  Return a revised Seurat object.
#' @export
#'
#'
filter_gene <- function(seu, min_spots=20, assay= NULL){

  if(is.null(assay)) assay <- DefaultAssay(seu)
  if(sum(dim(seu[[assay]]@counts))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@counts>0)>min_spots
    return(seu[names(gene_flag), ])
  }else if(sum(dim(seu[[assay]]@data))!=0){
    gene_flag <- Matrix::rowSums(seu[[assay]]@data>0)>min_spots
    return(seu[names(gene_flag), ])
  }else{
    stop("filter_gene: Seuat object must provide slots count or data in assay!")
  }


}
#' Filter spots for a Seurat object
#' @description Filter spots for a single-cell RNA sequencing or sptially-resolved transcriptomics data in the form of a Seurat object.
#' @param seu a Seurat object.
#' @param min_feature a non-negative integer, the least number of nonzero genes for a retained spot.
#' @param assay a string, the assay to perform filtering.
#' @return  Return a revised Seurat object.
#' @export
#'
#'
filter_spot <- function(seu, min_feature=0, assay=NULL){ # each spots at least include 1 non-zero features

  if(is.null(assay)) assay <- DefaultAssay(seu)
  col_name <- paste0("nFeature_",assay)
  idx <- seu@meta.data[,col_name] > min_feature
  seu[, idx]
  # subset(seu, subset = nFeature_RNA > min_feature)
}


.validInput <- function (input = NULL, name = NULL, valid = NULL)
{
  valid <- unique(valid)
  if (is.character(valid)) {
    valid <- tolower(valid)
  }
  else {
    stop("Validator must be a character!")
  }
  if (!is.character(name)) {
    stop("name must be a character!")
  }
  if ("null" %in% tolower(valid)) {
    valid <- c("null", valid[which(tolower(valid) !=
                                     "null")])
  }
  av <- FALSE
  for (i in seq_along(valid)) {
    vi <- valid[i]
    if (vi == "integer" | vi == "wholenumber") {
      if (all(is.numeric(input))) {
        cv <- min(abs(c(input%%1, input%%1 - 1)), na.rm = TRUE) <
          .Machine$double.eps^0.5
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "null") {
      cv <- is.null(input)
    }
    else if (vi == "bool" | vi == "boolean" |
             vi == "logical") {
      cv <- is.logical(input)
    }
    else if (vi == "numeric") {
      cv <- is.numeric(input)
    }
    else if (vi == "vector") {
      cv <- is.vector(input)
    }
    else if (vi == "matrix") {
      cv <- is.matrix(input)
    }
    else if (vi == "sparsematrix") {
      cv <- is(input, "dgCMatrix")
    }
    else if (vi == "character") {
      cv <- is.character(input)
    }
    else if (vi == "factor") {
      cv <- is.factor(input)
    }
    else if (vi == "rlecharacter") {
      cv1 <- is(input, "Rle")
      if (cv1) {
        cv <- is(input@values, "factor") || is(input@values,
                                               "character")
      }
      else {
        cv <- FALSE
      }
    }
    else if (vi == "palette") {
      cv <- all(.isColor(input))
    }
    else if (vi == "timestamp") {
      cv <- is(input, "POSIXct")
    }
    else if (vi == "dataframe" | vi == "data.frame" |
             vi == "df") {
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
    }
    else if (vi == "fileexists") {
      cv <- all(file.exists(input))
    }
    else if (vi == "direxists") {
      cv <- all(dir.exists(input))
    }
    else if (vi == "granges" | vi == "gr") {
      cv <- is(input, "GRanges")
    }
    else if (vi == "grangeslist" | vi == "grlist") {
      cv <- .isGRList(input)
    }
    else if (vi == "list" | vi == "simplelist") {
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
    }
    else if (vi == "bsgenome") {
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text = input))
      }, error = function(e) {
        FALSE
      })
      cv <- any(cv1, cv2)
    }
    else if (vi == "se" | vi == "summarizedexperiment") {
      cv <- is(input, "SummarizedExperiment")
    }
    else if (vi == "seurat" | vi == "seuratobject") {
      cv <- is(input, "Seurat")
    }
    else if (vi == "txdb") {
      cv <- is(input, "TxDb")
    }
    else if (vi == "orgdb") {
      cv <- is(input, "OrgDb")
    }
    else if (vi == "bsgenome") {
      cv <- is(input, "BSgenome")
    }
    else if (vi == "parallelparam") {
      cv <- is(input, "BatchtoolsParam")
    }
    else if (vi == "archrproj" | vi == "archrproject") {
      cv <- is(input, "ArchRProject")
    }
    else {
      stop("Validator is not currently supported by ArchR!")
    }
    if (cv) {
      av <- TRUE
      break
    }
  }
  if (av) {
    return(invisible(TRUE))
  }
  else {
    stop("Input value for '", name, "' is not a ",
         paste(valid, collapse = ","), ", (",
         name, " = ", class(input), ") please supply valid input!")
  }
}




# Tools for visualization -------------------------------------------------


range01 <- function(x, ...){
  if(all(x==0)){
    return(x)
  }else{
    (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  }

}
#' Transfer a vector to a list
#'
#' @param y_int a vector with length n.
#' @param nvec an integer vector, specify the length of each component in the transferred list. \code{n=sum(nvec)}.
#' @export
#' @return return a list with each component vector
#' @examples
#' vec2list(c(rep(1,3), rep(2, 5)), nvec=c(3, 5))
#'
#'

vec2list <- function(y_int, nvec){
  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")

  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}

#' Transfer a matrix to a list
#'
#' @param z_int a vector with number of rows n.
#' @param nvec an integer vector, specify the length of each component in the transferred list. \code{n=sum(nvec)}.
#' @export
#' @return return a list with each componnet matrix.
#' @examples
#' mat2list(matrix(1, nrow=3+5, ncol=2), nvec=c(3,5))
#'
mat2list <- function(z_int, nvec){

  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

#' Transfer matrix list to a matrix
#'
#' @param XList a list with each component matrix
#' @return return a rbinded matrix
#' @export
#' @examples
#' matlist2mat(list(matrix(1, nrow=2, ncol=2), matrix(2, nrow=4, ncol=2)))
#'
#'
matlist2mat <- function(XList){
  # transfer a matrix list to a matrix stacked by rows.
  r_max <- length(XList)
  X0 <- XList[[1]]
  if(r_max>1){
    for(r in 2:r_max){
      X0 <- rbind(X0, XList[[r]])
    }
  }

  return(X0)
}


get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}


get_indexList <- function(alist){
  nsample <- length(alist)
  nr <- 0
  indexList <- list()
  for(i in 1:nsample){
    indexList[[i]] <- (nr+1):(nrow(alist[[i]] )+nr)
    nr <- nr + nrow(alist[[i]] )
  }
  return(indexList)
}

replace_cluster <- function(clusters, new_names){

  uni_clus <- unique(clusters)
  new_names <- new_names[uni_clus]
  n_clus <- length(unique(clusters))
  new_vec <- clusters
  for(i in 1:n_clus){
    new_vec[clusters==uni_clus[i]] <- new_names[i]
  }
  return(new_vec)
}


#' Prepare the reference panel data for the RCTD deconvolution analysis
#'
#' @param sc_obj a Seurat object, denote the single cell RNA reference panel data.
#' @param clust_vr a string, the cluster name in the meta data of sc_obj.
#' @return return a list with three components: meta_data, cell_type_dict and dge.
#' @export
#'
RCTD_structure <- function(sc_obj, clust_vr) {

  sc_obj[["Name"]] = sc_obj@meta.data[, clust_vr]

  # Cell type dictionary between cluster and cell-type
  ct <- unique(sc_obj@meta.data[, clust_vr])
  df_ct <- data.frame("Cluster" = 1:length(ct),
                      "Name" = ct)
  metadata <- sc_obj@meta.data %>%
    # Rownames to columns must be before left join since after it the rownames are erased
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(df_ct, by = c("Name" = "Name")) %>%
    # Change names to 鈥渂arcode鈥? 鈥渃luster鈥? 鈥渘UMI鈥?
    mutate(
      cluster = Cluster,
      nUMI = nCount_RNA
    ) %>%
    dplyr::select(barcode, cluster, nUMI)

  expr_mtrx <- sc_obj@assays$RNA@counts

  return(list("meta_data" = metadata,
              "cell_type_dict" = df_ct,
              "dge" = expr_mtrx))
}

#' Perform spatial deconvolution analysis using RCTD model
#'
#' @param sc_data a Seurat object, denote the single cell RNA reference panel data.
#' @param cell_type_columns a string, the cluster name in the meta data of sc_obj.
#' @param st_count a sparse matrix, denote the count matrix of the spatial transcriptomics data to be deconvoluted.
#' @param st_coord a matrix with two columns, the spatial coordinates.
#' @param output_dir a string, specify the output directory.
#' @param output_name a string, the file name used.
#' @return write the results as csv format in the output_dir.
#' @export
#'
#'
RCTD_run <- function(sc_data,cell_type_columns,st_count,st_coord,output_dir,output_name){
  #sc_data: seurat object; Rownames should be genes and colnames represent cells/barcodes names.
  #st_data: Rownames should be genes and colnames represent barcodes/pixel names.


  print('Prepare SC data')
  sc_ls <- RCTD_structure(sc_obj = sc_data,clust_vr = cell_type_columns)   #use which column as cell type
  meta_data <- sc_ls[[1]]
  sc_counts <- sc_ls[[3]]
  cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- as.numeric(meta_data$nUMI); names(nUMI) <- meta_data$barcode # create nUMI named list

  reference <- Reference(sc_counts, cell_types, nUMI)

  print('Prepare ST data')
  puck <- SpatialRNA(st_coord, st_count)  # rownames are barcodes/pixel names

  print('Run RCTD')
  # Creating and running RCTD
  myRCTD <- create.RCTD(spatialRNA = puck,
                        reference = reference,
                        max_cores = 1,
                        CELL_MIN = 18)
  myRCTD <- run.RCTD(RCTD = myRCTD,
                     doublet_mode = 'doublet')
  deconv_results <- myRCTD@results

  print('Save result')
  dir.create(file.path(output_dir), showWarnings = FALSE)
  saveRDS(object = myRCTD,paste(output_dir,"/deconv_results_",output_name,".RDS",sep=""))

  #normalize the cell type proportions to sum to 1.
  rctd_deconv <- sweep(deconv_results$weights, 1, rowSums(deconv_results$weights), '/')
  rctd_deconv = as.matrix(rctd_deconv)
  rctd_deconv[rctd_deconv<1e-3] = 0
  colnames(rctd_deconv) = sc_ls[[2]][,'Name']

  ### Save results
  write.csv(rctd_deconv, file = paste(output_dir,"/output_weights_",output_name,".csv",sep=""))
}


align_colors <- function(ClusterMat,base_cluster, cols_cluster){
  # j <- 1
  if(length(cols_cluster) < max(apply(ClusterMat, 2, max))) stop("Number of cols is not enough!")
  givecolors <- function(mapID, base_cols_cluster, other_cols){
    m <- ncol(mapID)
    color_ID <- unname(base_cols_cluster)
    row_clusterID <- as.numeric(names(mapID))
    row_clusterCols <- rep("", m)
    col_used_count <- rep(0, length(base_cols_cluster))
    tmp_base_cluster <- NULL
    for(i in 1: m){
      # i <- 1

      if((mapID[1,i] %in% color_ID) &&  (max(mapID[2, mapID[1,] == mapID[1,i]]) == mapID[2,i])){
        col_chose <- names(base_cols_cluster)[mapID[1, i]]
        row_clusterCols[i] <- col_chose
        tmp_base_cluster <- c(tmp_base_cluster, mapID[1, i])
        color_ID <- setdiff(color_ID, mapID[1, i])

      }else{
        row_clusterCols[i] <- other_cols[1]
        other_cols <- other_cols[-1]
      }

      col_used_count[mapID[1, i]] <- col_used_count[mapID[1, i]] + 1


    }
    if(any(is.na(row_clusterCols))) row_clusterCols[is.na(row_clusterCols)] <- names(base_cols_cluster)[color_ID]
    return(row_clusterCols)
  }
  base1 <- sort(unique(base_cluster))
  n_base <- length(base1)
  base_cols_cluster <- base1
  names(base_cols_cluster) <- cols_cluster[1:n_base]
  other_cols <- cols_cluster[-(1:n_base)]
  colorList <- list()
  for(j in 1:ncol(ClusterMat)){
    # j <- 2
    stab <- table(ClusterMat[,j], base_cluster)
    mapID <- rbind(apply(stab, 1, which.max), apply(stab, 1, max))
    colorList[[j]] <- givecolors(mapID, base_cols_cluster, other_cols)
  }
  return(colorList)
}




# ### Tools used in DR.SC -------------------------------------------------

# One sample function
#' SingleCellExperiment object to Seurat object
#' @description Transfer SingleCellExperiment object to a Seurat object for preparation for DR.SC model fitting; see our [DR.SC package website](https://feiyoung.github.io/DR.SC/index.html) for more usage of DR.SC.
#' @param sce a SingleCellExperiment object, at least including the raw gene count expression matrix.
#' @param verbose an optional logical value, whether output the information.
#' @return Return a Seurat object.
#' @examples
#' library(SingleCellExperiment)
#' dir <- system.file(
#'   file.path("extdata", "10xVisium", "section1"),
#'   package = "SpatialExperiment")
#'
#' # read in counts
#' fnm <- file.path(dir, "raw_feature_bc_matrix")
#' sce <- DropletUtils::read10xCounts(fnm)
#' colnames(sce) <- paste0("cell", 1:ncol(spe))
#'
#'
#' seu <- sce2seurat(sce)
#' head(seu)
#'
#' ## Fit DR-SC model
#' library(DR.SC)
#' library(Seurat)
#' seu <- NormalizeData(seu)
#' seu <- FindVariableFeatures(seu)
#' seu <- DR.SC(seu = seu, K=4, platform = "scRNAseq")
#' @export
sce2seurat <- function(sce, verbose= TRUE){

  ## Transfer SingleCellExperiment object to a Seurat object for preparation for DR.SC model fitting.
  if(verbose){
    message("Transfer SingleCellExperiment object to a Seurat object")
    message("preparation for  model fitting")
  }

  require(SingleCellExperiment)
  count <- counts(sce)
  meta_data <- as.data.frame(colData(sce))
  require(Seurat)
  seu <- CreateSeuratObject(counts=count, meta.data = meta_data)

  return(seu)
}

#' SpatialExperiment object to Seurat object
#' @description Transfer SpatialExperiment object to a Seurat object for preparation for DR.SC model fitting.
#' @param spe a SpatialExperiment object, at least including the raw gene count expression matrix ans sptial coordinates.
#' @param verbose an optional logical value, whether output the information.
#' @return Return a Seurat object, where the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col".
#' @examples
#' dir <- system.file(
#' file.path("extdata", "10xVisium", "section1"),
#' package = "SpatialExperiment")
#'
#' #' read in counts
#' fnm <- file.path(dir, "raw_feature_bc_matrix")
#' sce <- DropletUtils::read10xCounts(fnm)
#'
#' #' read in image data
#' img <- readImgData(
#'   path = file.path(dir, "spatial"),
#'   sample_id="foo")
#'
#' #' read in spatial coordinates
#' fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
#' xyz <- read.csv(fnm, header = FALSE,
#'                 col.names = c(
#'                   "barcode", "in_tissue", "array_row", "array_col",
#'                   "pxl_row_in_fullres", "pxl_col_in_fullres"))
#'
#' #' construct observation & feature metadata
#' rd <- S4Vectors::DataFrame(
#'   symbol = rowData(sce)$Symbol)
#'
#' #' construct 'SpatialExperiment'
#' (spe <- SpatialExperiment(
#'   assays = list(counts = assay(sce)),
#'   colData = colData(sce), rowData = rd, imgData = img,
#'   spatialData=DataFrame(xyz),
#'   spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"),
#'   sample_id="foo"))
#'
#' colnames(spe) <- paste0("spot", 1:ncol(spe))
#'
#'
#' seu <- spe2seurat(spe)
#' head(seu)
#'
#' #'#' Fit DR-SC model
#' library(DR.SC)
#' library(Seurat)
#' seu <- NormalizeData(seu)
#' seu <- FindVariableFeatures(seu)
#' seu <- DR.SC(seu = seu, K=4)
#' @export
spe2seurat <- function(spe, verbose= TRUE){

  ## Transfer SpatialExperiment object to a Seurat object for preparation for DR.SC model fitting.
  if(verbose){
    message("Transfer SpatialExperiment object to a Seurat object")
  }
  require(Seurat)
  suppressPackageStartupMessages(require(SpatialExperiment))
  if(is.null(colnames(spe))) stop("spe2seurat: Check argument spe! spe must have colnames!")
  if(is.null(row.names(spe))) stop("spe2seurat: Check argument spe! spe must have row.names!")

  col_meta_data <- as.data.frame(colData(spe))
  meta.data <- cbind.data.frame(col_meta_data,data.frame(row=spatialCoords(spe)[,1],
                                           col=spatialCoords(spe)[,2]))


  ret <- CreateSeuratObject(
    counts=assays(spe)$counts,
    meta.data=meta.data
  )
  return(ret)
}


## spe2seuratList(spe)
# Tools used in PRECAST ---------------------------------------------------

#' SpatialExperiment object to Seurat list object
#' @description Transfer SpatialExperiment object to a Seurat list object for preparation for PRECAST model fitting; see our [PRECAST package website](https://feiyoung.github.io/PRECAST/index.html) for more usage of PRECAST.
#' @param spe a SpatialExperiment object or a list consisting of SpatialExperiment objects. If spe is a SpatialExperiment object, it must at least contain the batch(sampel) id in the colData (i.e., sample_id), the raw gene count expression matrix and spatial coordinates. And the batch must be specified (i.e., batch='sample_id'). If spe is a list consisting of multiple SpatialExperiment objects, then each object represents a data batch.
#' @param batch a optional argument, NULL or a string. Only if spe is a list, batch can be NULL.
#' @param verbose an optional logical value, whether output the information.
#' @export
#' @examples
#' suppressPackageStartupMessages( library(SpatialExperiment))
#' dir <- system.file(
#'   file.path("extdata", "10xVisium", "section1"),
#'   package = "SpatialExperiment")
#'
#' # read in counts
#' fnm <- file.path(dir, "raw_feature_bc_matrix")
#' sce <- DropletUtils::read10xCounts(fnm)
#'
#' # read in image data
#' img <- readImgData(
#'   path = file.path(dir, "spatial"),
#'   sample_id="foo")
#'
#' # read in spatial coordinates
#' fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
#' xyz <- read.csv(fnm, header = FALSE,
#'                 col.names = c(
#'                   "barcode", "in_tissue", "array_row", "array_col",
#'                   "pxl_row_in_fullres", "pxl_col_in_fullres"))
#'
#' # construct observation & feature metadata
#' rd <- S4Vectors::DataFrame(
#'   symbol = rowData(sce)$Symbol)
#'
#' # construct 'SpatialExperiment'
#' (spe <- SpatialExperiment(
#'   assays = list(counts = assay(sce)),
#'   colData = colData(sce), rowData = rd, imgData = img,
#'   spatialData=DataFrame(xyz),
#'   spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"),
#'   sample_id="foo"))
#' colnames(spe) <- paste0("spot", 1:ncol(spe))
#'
#'
#' ## Transfer a SpatialExperiment to a seuList
#' colData(spe)$batch_id <- rep(c("a", "b"), each=25)
#' seuList1 <- spe2seuratList(spe, batch = 'batch_id')
#' seuList1
#' ## Transfer a list of SpatialExperiment to a seuList
#' seuList2 <- spe2seuratList(list(spe, spe))
#' seuList2
#' ## Create PRECAST object
#' library(PRECAST)
#' Pobject <- CreatePRECASTObject(seuList=seuList2, selectGenesMethod = "HVGs",
#'                                premin.spots = 0, premin.features = 0,
#'                                postmin.features = 0, postmin.spots = 0, verbose = F)
#'
#' @return Return a list consisting of Seurat objects, where the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col".
spe2seuratList <- function(spe, batch=NULL, verbose=TRUE){

  ## Transfer one SpatialExperiment object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.

  if(verbose){
    message("Transfer one SpatialExperiment object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.")
  }
   require(purrr)
  if(is.null(batch) && (!is.list(spe))) stop("If batch=NULL, spe must be a list with each component to be SpatialExperiment object!")

  if(!is.list(spe) && (!is.null(batch))){
    uniq_sample_id <- unique(colData(spe)[,batch])
    # Create a seurate object for each unique sample_id
    seuList <- map(uniq_sample_id,
        .f = function(smp_id, spe){
          # browser()
          ret_spe <- spe[, colData(spe)[,batch] == smp_id]
          ret_seurat <- spe2seurat(ret_spe)

          return(ret_seurat)
        },
        spe = spe)


  }

  if(is.list(spe)){
    seuList <- pbapply::pblapply(spe, spe2seurat, verbose=verbose)
  }

  return(seuList)
}

#' Seurat object to seuratList object
#' @description  Transfer Seurat object to a Seurat list object for preparation for PRECAST model fitting. see our [PRECAST package website](https://feiyoung.github.io/PRECAST/index.html) for more usage of PRECAST.
#' @param seu a Seurat object with multiple data bathes information. It must at least contain the batch(sampel) id in the meta.data (i.e., sample_id), the raw gene count expression matrix ans sptial coordinates. And the batch must be specified (i.e., batch='sample_id').
#' @param batch a string to specify the data batch field in the meta.data of seu.
#' @param verbose an optional logical value, whether output the information.
#' @return Return a list of multiple Seurat objects, where for each Seurat object, the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col".
#' @examples
#' suppressPackageStartupMessages( library(SpatialExperiment))
#' dir <- system.file(
#'   file.path("extdata", "10xVisium", "section1"),
#'   package = "SpatialExperiment")
#'
#' # read in counts
#' fnm <- file.path(dir, "raw_feature_bc_matrix")
#' sce <- DropletUtils::read10xCounts(fnm)
#'
#' # read in image data
#' img <- readImgData(
#'   path = file.path(dir, "spatial"),
#'   sample_id="foo")
#'
#' # read in spatial coordinates
#' fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
#' xyz <- read.csv(fnm, header = FALSE,
#'                 col.names = c(
#'                   "barcode", "in_tissue", "array_row", "array_col",
#'                   "pxl_row_in_fullres", "pxl_col_in_fullres"))
#'
#' # construct observation & feature metadata
#' rd <- S4Vectors::DataFrame(
#'   symbol = rowData(sce)$Symbol)
#'
#' # construct 'SpatialExperiment'
#' (spe <- SpatialExperiment(
#'   assays = list(counts = assay(sce)),
#'   colData = colData(sce), rowData = rd, imgData = img,
#'   spatialData=DataFrame(xyz),
#'   spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"),
#'   sample_id="foo"))
#'
#' colnames(spe) <- paste0("spot", 1:ncol(spe))
#'
#'
#' seu <- spe2seurat(spe)
#' head(seu)
#' seu2 <- seu
#' seu2$sample_id <- paste0(seu2$sample_id, "2")
#' library(SeuratObject)
#' seu2 <- RenameCells(seu2,  paste0(colnames(seu2), "2") )
#' seu_all <- merge(seu, seu2)
#' head(seu_all@meta.data)
#' table(seu_all$sample_id)
#' seuList <- seu2seuList(seu_all, batch='sample_id')
#' seuList
#'
#' ## Create PRECAST object
#' library(PRECAST)
#' Pobject <- CreatePRECASTObject(seuList=seuList, selectGenesMethod = "HVGs",
#'                                premin.spots = 0, premin.features = 0,
#'                                postmin.features = 0, postmin.spots = 0, verbose = F)

#' @export
seu2seuList <- function(seu, batch, verbose=TRUE){


  ## Transfer one Seurat object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.
  if(verbose){
    message("Transfer one Seurat object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.")
  }
    require(purrr)
    meta_data <- seu@meta.data
    uniq_sample_id <- unique(meta_data[,batch])
    # Create a seurate object for each unique sample_id
    map(uniq_sample_id,
        .f = function(smp_id, seu){
          # browser()
          ret_seurat <- seu[, meta_data[,batch] == smp_id]
          return(ret_seurat)
        },
        seu = seu)

}


AddCoord2Reduc <- function(seu, spatial_coords= c("row", "col"), embed_name= "position"){

   require(PRECAST)
   embed <- as.matrix(seu@meta.data[,spatial_coords])
   #embed_name <- "position"
   seu <- Add_embed(embed = embed, seu=seu, embed_name = embed_name, assay= DefaultAssay(seu))
   redu <- CreateDimReducObject(embeddings = embed,
                               key = paste0(toupper(embed_name),"_"), assay = assay)
   seu@reductions[[embed_name]] <- redu
   return(seu)
}

# Generate simulated data using splatter package-------------------------------------------

generate_count <- function(seu, annotated_label, NumSpatialDomain=7, NumBatches = 3,
                           sim_seed=1,J = 2000, batch_facLoc=0.1, batch_facScale = 0.1){

  library(SingleCellExperiment)
  library(splatter)
  ## read position and annotation label from real data dlpfc

  #pos = seu@meta.data[,spatial_coords]
  #y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  y <- as.numeric(seu@meta.data[,annotated_label])
  assay <- DefaultAssay(seu)
  cnts = as.matrix(seu[[assay]]@counts)
  n_spots <- ncol(cnts)
  init_params <- splatEstimate(cnts)


  #batch_facLoc = 0.1 ## Batch effects in location
  #batch_facScale = 0.1
  #C = length(unique(y)) # the number of spatial Domains

  C <- NumSpatialDomain
  I = NULL
  N = n_spots*2
  L = NumBatches ## number of NumBatches

  debug = FALSE

  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  group_prob <- as.vector(table(y)/length(y))


  params <- setParams(
    init_params,
    batchCells = rep(N, L), # 3N here represents a large number such that
    # we have sufficient cells of each type to be
    # allocated to the spatial transcriptomics data
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale,
    nGenes = J,
    group.prob = group_prob,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)
  return(sim_groups)

}

simu_seuList <- function(seu, annotated_label, spatial_coords= c("row", "col"),
                             NumBatches = 3, NumSpotsRatio=0.9, ngenes = 2000, sim_seed=1,
                             batch_facLoc=0.1, batch_facScale = 0.1){


  # annotated_label= "annotation"; spatial_coords= c("row", "col");
  # NumBatches = 3; NumSpotsRatio=0.2; ngenes = 20; sim_seed=1;
  # batch_facLoc=0.1; batch_facScale = 0.1

  if(NumSpotsRatio<0.1) stop("simu_seuList:: Check argument: NumSpotsRatio! this argument must be greater than 0.1!")
  require(Seurat)
  ### Spatial domain annotations
  y <- as.numeric(seu@meta.data[,annotated_label])
  NumSpatialDomain <- length(unique(y))
  sim_groups <- generate_count(seu,annotated_label, NumSpatialDomain ,sim_seed = sim_seed, J= ngenes,
                               batch_facLoc=batch_facLoc, batch_facScale = batch_facScale)
  ## reorder the sample to  be the same as y
  meta_data_all <- colData(sim_groups)


  Groupy = paste0("Group", y)
  nbatch <- NumBatches

  pos <- as.matrix(seu@meta.data[,spatial_coords])
  gen_pos_id <- function(iseed, pos){
    set.seed(iseed)
    quantx <- runif(1, NumSpotsRatio-0.1 ,NumSpotsRatio+0.1)
    quanty <- runif(1, NumSpotsRatio-0.1 ,NumSpotsRatio+0.1)
    pos_id <- which(pos[,1]< quantile(pos[,1], quantx) & pos[,2]< quantile(pos[,2], quanty)  )
    return(pos_id)
  }

  idxList_pos <- list()
  for(i_batch in 1: nbatch){

    # list(1:length(y), which(pos[,1]< quantile(pos[,1], 0.9)), which(pos[,2]< quantile(pos[,2], 0.9)))
    idxList_pos[[i_batch]] <- gen_pos_id(i_batch, pos)
  }
  posList <- lapply(idxList_pos, function(idx) pos[idx, ])


  yList <- lapply(idxList_pos, function(idx) y[idx])

  idxList_sim <- list()
  for(r in 1:nbatch){
    message("r = ", r)
    y_tmp <- yList[[r]]
    num_each_celltype = table(y_tmp)
    i_vec <- as.numeric(names(num_each_celltype))
    idx1 = rep(0, length(y_tmp))
    for (i in i_vec){
      idx1[y_tmp==i] = which(meta_data_all$Group == paste0("Group",i) & meta_data_all$Batch==paste0("Batch", r))[1:num_each_celltype[i]]
    }
    idxList_sim[[r]] <- idx1 ## align with posList

  }

  sceList_sim <- lapply(idxList_sim, function(idx) sim_groups[,idx])
  ## Add spatial coordinates
  sceList_sim <- lapply(1: nbatch, function(r){
    sce <- sceList_sim[[r]]
    colData(sce)$row <- posList[[r]][,1]
    colData(sce)$col <- posList[[r]][,2]
    return(sce)
  })

  seuList_use <- pbapply::pblapply(sceList_sim, sce2seurat)
  return(seuList_use)
}




# Metrics -----------------------------------------------------------------
#' Compute adjusted Rand index and normalized mutual information
#'
#' @param hy a vector, the predicted clusters
#' @param y a vector, the true clusters
#' @param type a string, the type of cluster metric, one of "ARI" and "NMI", default as "ARI".
#' @export
#' @examples
#' y <- rep(1:4, each=100)
#' hy <- rep(1:5, each=80)
#' cluster_metric(y, hy)
#'
#'
cluster_metric <- function(hy, y, type=c('ARI',"NMI")){

  require(mclust)
  require(aricode)
  type <- match.arg(type)
  switch(type,
         ARI= adjustedRandIndex(hy, y),
         NMI = NMI(as.vector(hy), y))
}

