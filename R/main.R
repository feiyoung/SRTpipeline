
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
# github_document("README.Rmd")
# library(pkgdown)
# build_site()
# build_home()
# build_reference()
# build_article("get_started")
# build_article("dlpfc_tutorial")



# Tools used in preprocessing ---------------------------------------------

firstup <- function(x) {
  ## First letter use upper capital
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
transferGeneNames <- function(genelist, now_name = "ensembl", to_name="symbol",
                              species="Human", Method='biomaRt'){

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

  if(toupper(to_name) == 'SYMBOL')
    transferredNames <- firstup(transferredNames)

  flag_na <- is.na(transferredNames)
  if(any(flag_na))
    transferredNames[flag_na] <- genelist[flag_na]

  return(transferredNames)
}





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
get_top_pathway <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }

  return(df_sub)
}


barPlot_enrich <- function(top_dat, source='Ont', term_name="Term", nlog10P='nlog10P',
                           bar_width=0.8, base_size=20, font_family='serif', cols= ggthemes::canva_pal()(4)){
  # source='Ont'; term_name="Term"; nlog10P='nlog10P'
  require(ggplot2) # y=term_name,
  order_idx <- order(top_dat[,nlog10P])
  top_dat <- top_dat[order_idx,]
  top_dat[, term_name] <- factor(top_dat[, term_name], levels=top_dat[order_idx,term_name])
  p1 <- ggplot(data=top_dat, aes_string(x=term_name,y=nlog10P, fill=source)) +
    scale_fill_manual(values=cols)+
    geom_bar(position = "dodge", stat="identity",width =bar_width)+ coord_flip() +
    theme_classic() + theme(text=element_text(size=base_size, family=font_family))
  return(p1)
}

SpaPlot <- function(seuInt, batch=NULL, item=NULL, point_size=2,text_size=12,
                    cols=NULL,font_family='', border_col="gray10",
                    fill_col='white', ncol=2, combine = TRUE, title_name="Sample"){

  ## Check arguments input
  require(ggplot2)
  if(!inherits(seuInt, "Seurat")) stop("SpaPlot: check argument: seuInt! it must be a Seurat Object.")
  if(!('batch' %in% colnames(seuInt@meta.data))) seuInt$batch <- 1
  seuInt@meta.data$ident <- Idents(seuInt)
  if(is.null(item)) item <- "ident"

  if(item %in% colnames(seuInt@meta.data)){
    if(!is.factor(seuInt@meta.data[, item])) seuInt@meta.data[, item] <- factor(seuInt@meta.data[, item])

    seuInt@meta.data$tmp_item_id <- as.numeric(seuInt@meta.data[, item])
  }else if(!(toupper(item) %in% c("RGB_UMAP", "RGB_TSNE"))){
    stop("SpaPlot: check the value of argument: item! It is not the colname of meta.data of seuInt!")
  }
  if(is.null(cols)&& toupper(item) != "RGB_UMAP" && toupper(item)!="RGB_TSNE"){
    # to determine the number of colors

    nn <- length(unique(seuInt@meta.data[,item]))
    cols <- PRECAST:::gg_color_hue(nn)
  }


  if(!is.vector(cols) && item != "RGB_UMAP" && item!="RGB_tSNE")
    stop("Check argument: cols! it must be a vector object.")


  ###Finish  Check  of arguments


  if(is.null(batch)){
    batch_vec <- unique(seuInt$batch)
  }else{
    batch_vec <- batch
  }


  if(length(batch_vec)<2) ncol <- 1




  pList <- list()
  k <- 1
  item_null_flag <- FALSE
  for(batchi in batch_vec){
    # batchi <- (batch_vec[2])
    seu <- subset(seuInt, batch==batchi)
    meta_data <- seu@meta.data
    meta_data$ident <- Idents(seu)

    embed_use <- seu@reductions$position@cell.embeddings
    if(item %in% colnames(meta_data)){

      sort_id <- sort(unique(meta_data[, 'tmp_item_id']))
      p1 <- plot_scatter(embed_use, meta_data, label_name=item,
                         point_size=point_size, cols =cols[sort_id])
    }else if(toupper(item)=="RGB_UMAP"){
      p1 <- plot_RGB(embed_use, seu@reductions$UMAP3@cell.embeddings, pointsize = point_size)
    }else if(toupper(item)=="RGB_TSNE"){
      p1 <- plot_RGB(embed_use, seu@reductions$tSNE3@cell.embeddings, pointsize = point_size)
    }
    p1 <- p1 + PRECAST:::mytheme_graybox(base_size = text_size, base_family = font_family, bg_fill = fill_col,
                                         border_color = border_col)
    if(!is.null(title_name)){
      p1 <- p1 + ggtitle(label=paste0(title_name, batchi))
    }

    pList[[k]] <- p1
    k <- k + 1
    if(item_null_flag){
      item <- NULL
    }
  }
  if(combine){

    pList <- patchwork::wrap_plots(pList, ncol=ncol)
  }
  return(pList)
}

featurePlot <- function(seu, feature, reduction="position",
                        assay=NULL, cols=c('#3AB370',"#FD1593"),
                        pt_size=1,title_size=16, quant = 0.9){
  require(ggplot2)
  if(is.null(assay)) assay <- DefaultAssay(seu)
  dat <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  x <- as.vector(seu[["RNA"]]@data[feature,])
  dat$Expression <- range01(x )
  med <- quantile(dat$Expression, quant)
  ggplot(data=dat, aes(x=POSITION_1 , y=POSITION_2, color=Expression)) + geom_point(size=pt_size) +
    scale_colour_gradient2(low = "#0571B099", mid = "white", high = "#CA0020", midpoint = med)+
    PRECAST:::mytheme_graybox() +
    ggtitle(feature) + theme(title =element_text(size=title_size, color=1, face='italic'),
                             legend.position = 'none')

}





range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


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



get_top_pathway1 <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }

  return(df_sub)
}

coordinate_rotate <- function(pos, theta=0){# counter-clock rotation
  pos_new <- pos
  pos_new[,1] <- pos[,1]*cos(theta) - pos[,2]*sin(theta)
  pos_new[,2] <- pos[,1]*sin(theta) + pos[,2]*cos(theta)
  return(pos_new)
}
rotate_angles <- c(90.1, -40.5)


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

spe2seurat <- function(spe, verbose= TRUE){

  ## Transfer SpatialExperiment object to a Seurat object for preparation for DR.SC model fitting.
  if(verbose){
    message("Transfer SpatialExperiment object to a Seurat object")
    message("preparation for  model fitting")
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


# Multiple sample

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


