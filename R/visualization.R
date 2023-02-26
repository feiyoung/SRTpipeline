

.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Choose color schema from a palette
#'
#' @param palettes_name a string, the palette name, one of "Nature 10", "Light 13", "Classic 20", "Blink 23" and "Hue n", default as 'Nature 10'.
#' @param n_colors a positive integer, the number of colors.
#' @param alpha a positive real, the transparency of the color.
#' @param plot_colors a logical value, whether plot the selected colors.
#' @export
#' @importFrom colorspace adjust_transparency
#' @examples
#' chooseColors()
#'
chooseColors <- function(palettes_name= c("Nature 10", "Light 13", "Classic 20", "Blink 23", "Hue n"), n_colors = 7,
                         alpha=1, plot_colors=FALSE){

  require(colorspace)
  palettes_name <- match.arg(palettes_name)
  colors <- if(palettes_name == "Classic 20"){
    require(ggthemes)
    # palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
    pal1 <- tableau_color_pal(palettes_name)
    pal1(n_colors)
  }else if(palettes_name == "Nature 10"){
     cols <- c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",
               "#70B5B0", "#DFE0EE" ,"#DFCDE4", "#FACB12", "#f9decf")
     cols[1:n_colors]
  }else if(palettes_name == "Blink 23"){
    cols <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#FE0092", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
                        "#D35400", "#00eefd", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
                        "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#007CC8")
    cols[1:n_colors]
  }else if(palettes_name == "Light 13"){
    cols <-c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
       "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
       "#FF9896","#91D1C2", "#C7E9C0" ,
       "#6B6ECF", "#7B4173" )
    cols[1:n_colors]
  }else if(palettes_name == "Hue n"){
    .gg_color_hue(n_colors)
  }else{
    stop(paste0("chooseColors: check palettes_name! Unsupported palettes_name: ", palettes_name))
  }
  require(colorspace)
  colors_new = adjust_transparency(colors,   alpha = alpha)
  if(plot_colors){
    barplot(rep(1, length(colors_new)), axes = FALSE, space = 0, col = colors_new)
  }

  return(colors_new)
}

.plot_scatter <- function (embed, labels, label_name, axis_names = c("tSNE1","tSNE2"),
                           cols = NULL, no_guides = FALSE, do_density = FALSE,
                           no_axis_name = FALSE, pt_size = 1, point_alpha = 1,
                           nrow.legend = NULL, base_size=10,border_col='gray',
                           legend.position = "bottom", no_axis = FALSE){

  require(ggplot2)
  plt_df <- as.data.frame(embed)
  colnames(plt_df) = c("x", "y")
  plt_df[label_name] = labels
  if (is.null(cols)) {
    cluster <- as.vector(plt_df[[label_name]])
    ngrp <- length(unique(cluster))
    cols <- .gg_color_hue(ngrp)
  }
  plt <- plt_df %>% ggplot(aes_string("x", "y",
         col = label_name, fill = label_name)) +
    geom_point(size = pt_size, alpha = point_alpha) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols)+
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.position = legend.position,
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col))

  if (!is.null(nrow.legend)){
    plt <- plt + guides(color = guide_legend(nrow = nrow.legend,override.aes = list(stroke = 1,
                          alpha = 1, shape = 16, size = 4)))
  }else{
    plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1,
                                                                 alpha = 1, shape = 16, size = 4)),
                        alpha = "none")
  }



  if (no_axis) {
    plt <- plt + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                       axis.title.y = element_blank(), axis.ticks.y = element_blank())
  }
  else if (!no_axis_name) {
    plt <- plt + labs(x = axis_names[1], y = axis_names[2])
  }
  if (do_density)
    plt <- plt + geom_density_2d()
  if (no_guides)
    plt <- plt + guides(col = "none", fill = "none",
                        alpha = "none")
  if (legend.position == "bottom") {
    plt <- plt + theme(legend.direction = "horizontal")
  }
  else if (legend.position == "right") {
    plt <- plt + theme(legend.direction = "vertical")
  }
  return(plt)
}


### Plot spatial heatmap
#' Plot spatial cluster heatmap for each data batch
#'
#' @description
#' The function \code{EachClusterSpaHeatMap} is used to visualize  clusters on spatial coordinates with the RGB colors for each data batch.
#' @export
#' @param obj A SRTProject object.
#' @param batch a string, number or NULL, specify the batch name to plot, default as \code{NULL} implying all batches plotted.
#' @param cols a string vector, colors used in the plot.
#' @param title_name a string, the title of plot.
#' @param combine a logical value,  whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component. TRUE by default.
#' @param layout.dim a 2 vector,  the dimension in the layout of plots when \code{combine = TRUE}.
#' @param no_axis a logical value, whether displys the axis.
#' @param ... Arguments passed to other methods.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom dplyr `%>%`
#'
EachClusterSpaHeatMap <- function(obj, ...){
  UseMethod(generic = 'EachClusterSpaHeatMap', object = obj)
}

#' @rdname EachClusterSpaHeatMap
#' @method EachClusterSpaHeatMap list
EachClusterSpaHeatMap.list <- function(obj, cols=NULL, title_name = NULL, combine = TRUE,
   layout.dim = c(round(length(obj[[2]])/2), 2), no_axis = FALSE,pt_size = 2, point_alpha = 1,
   nrow.legend = NULL, base_size=10, border_col='gray',
   axis_names = c("Coord x","Coord y"),no_guides = FALSE,
   no_axis_name = FALSE,
   legend.position = "bottom", common.legend = TRUE){

  require(dplyr)
  if (length(obj) != 2)
      stop("EachClusterSpaHeatMap: check argument: obj! The list must constain tow elements. The first is a list of two column matrices to specify the position information, and the second specify the features information!")
    if (!is.list(obj[[1]]))
      stop("EachClusterSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")
    if (any(sapply(obj[[1]], ncol) != 2))
      stop("EachClusterSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")



    nT = length(obj[[1]])
    sample_name = names(obj[[1]])
    if (is.null(sample_name))
      sample_name = paste0("Batch", c(1:nT))
    pList <- vector("list", nT)

      if (!is.list(obj[[2]]))
        stop("EachClusterSpaHeatMap: check argument: obj! For scatter plot, the second element of list must be a list of vector which specify the features information!")
      if (any(!sapply(obj[[2]], is.vector)) && any(!sapply(obj[[2]], is.factor)) )
        stop("EachClusterSpaHeatMap: check argument: obj! For scatter plot, the second element of list must be a list of vector which specify the features information!")
      cluster = obj[[2]]
      if(is.factor(unlist(cluster))){
        uni_cluster <- cluster %>% unlist %>% levels;
      }else{
        uni_cluster <- cluster %>% unlist %>% as.vector %>% unique %>% sort
      }


      ngrp <- length(uni_cluster)
      cols <- cols[1:ngrp]
      if (is.null(cols)) {
        cols <- .gg_color_hue(ngrp)

      }

      names(cols) <-  uni_cluster


      label_name = names(obj)[2]
      if (is.null(label_name))
        label_name = "cluster"
      for (i in 1:nT) {
        cluster_fac_i <- factor(cluster[[i]])
        uni_cluster_tmp <- levels(cluster_fac_i)
        p1 <- .plot_scatter(obj[[1]][[i]], cluster_fac_i, label_name,
                            cols=cols[uni_cluster_tmp], axis_names = axis_names,
                             no_guides = no_guides, do_density = FALSE,
                            no_axis_name = no_axis_name, pt_size = pt_size,
                            point_alpha = point_alpha,nrow.legend = nrow.legend,
                            base_size=base_size,border_col=border_col,
                            legend.position = legend.position, no_axis = no_axis)
        if (!is.null(title_name)) {
          p1 <- p1 + ggtitle(label = paste0(title_name,
                                            sample_name[i]))
        }
        pList[[i]] <- p1
      }

    if (combine) {
      pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                                 nrow = layout.dim[1], legend = legend.position,
                                 common.legend = common.legend)
    }
    return(pList)
}

#' @rdname EachClusterSpaHeatMap
#' @method EachClusterSpaHeatMap SRTProject
#' @export
EachClusterSpaHeatMap.SRTProject <- function(obj, batch=NULL, cols= NULL, layout.dim = NULL,
                                         combine=TRUE,  no_axis=TRUE, ...){
   #obj <- SRTProj

  ## batch can be a vector representing the number of batches or names of batches

   posList <- mat2list(obj@spatialCoords, obj@sampleMetaData$NumOfSpots)
   names(posList) <- row.names(obj@sampleMetaData)
   clusterList <- vec2list(obj@clusters, obj@sampleMetaData$NumOfSpots)
   names(clusterList) <- row.names(obj@sampleMetaData)
   if(!is.null(batch)){
     posList <- posList[batch]
     clusterList <- clusterList[batch]
   }
   list2 <- list(pos=posList, cluster=clusterList)
   pList <- EachClusterSpaHeatMap(list2, cols=cols,combine=combine,
                              no_axis=no_axis,layout.dim=layout.dim, ...)
   return(pList)
}


.plot_RGB <- function(position, embed_3d, no_axis=TRUE, no_axis_name = TRUE,axis_names=c("x", "y"),
                      pt_size =2, base_size=15, border_col='gray'){


  # no_axis=TRUE; no_axis_name = TRUE;axis_names=c("x", "y");
  # pt_size =2; base_size=15; border_col='gray'
  # suppressMessages(require(ggplot2))

  info = as.data.frame(position)
  colnames(info) = c("sdimx","sdimy")


  r = (embed_3d[,1]-min(embed_3d[,1]))/(max(embed_3d[,1])-min(embed_3d[,1]))
  g = (embed_3d[,2]-min(embed_3d[,2]))/(max(embed_3d[,2])-min(embed_3d[,2]))
  b = (embed_3d[,3]-min(embed_3d[,3]))/(max(embed_3d[,3])-min(embed_3d[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  plt <- ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
    geom_point(size=pt_size) +
    scale_color_identity()+
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col))
  if (no_axis) {
    plt <- plt + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                       axis.title.y = element_blank(), axis.ticks.y = element_blank())
  }
  else if (!no_axis_name) {
    plt <- plt + labs(x = axis_names[1], y = axis_names[2])
  }


  plt
}

#' Plot spatial RGB heatmaps for a group of SRT data batches
#'
#' @description
#' The function \code{EachRGBSpaHeatMap} is used to visualize the RGB colors of 3-dimensional embeddings on spatial coordinates.
#' @export
#' @param obj A SRTProject object.
#' @param batch a string vector or positive integer vector, the batch names or IDs to be plotted.
#' @param plot_type a string,  the type of plot, one of "UMAP" and "tSNE".
#' @param layout.dim a 2 integer vector, the dimension in the layout of plots when \code{combine = TRUE}.
#' @param combine a logical value, whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component. TRUE by default.
#' @param title_name a string, the title of plot.
#' @param common.legend a logical value,  whether combine the legend of multiple plots. TRUE by default.
#' @param no_axis a logical value, whether display the axis in plot.
#' @param pt_size a positive real, the point size in plot.
#' @param base_size a positive real, the baseline text size in the plot
#' @param border_col a string, the color of border in plot.
#' @param axis_names a 2 string vector, the axis names in plot.
#' @param no_axis_name a logical value, whether display the axis name in plot.
#' @param ... Arguments passed to \code{\link{ggarrange}}.
#'
#' @importFrom ggpubr ggarrange

EachRGBSpaHeatMap <- function(obj, ...){
  UseMethod(generic = 'EachRGBSpaHeatMap', object = obj)
}
#' @rdname EachRGBSpaHeatMap
#' @method EachRGBSpaHeatMap list
EachRGBSpaHeatMap.list <- function(obj, title_name = NULL, combine = TRUE,
                                   layout.dim = c(round(length(obj[[2]])/2), 2), no_axis = FALSE,
                                 pt_size = 2,  base_size=10, border_col='gray',
                                   axis_names = c("Coord x","Coord y"),
                                   no_axis_name = FALSE, ...){
  require(dplyr)
  if (length(obj) != 2)
    stop("EachRGBSpaHeatMap: check argument: obj! The list must constain tow elements. The first is a list of two column matrices to specify the position information, and the second specify the features information!")
  if (!is.list(obj[[1]]))
    stop("EachRGBSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")
  if (any(sapply(obj[[1]], ncol) != 2))
    stop("EachRGBSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")



  nT = length(obj[[1]])
  sample_name = names(obj[[1]])
  if (is.null(sample_name))
    sample_name = paste0("Batch", c(1:nT))
  pList <- vector("list", nT)

  if (!is.list(obj[[2]]))
    stop("EachRGBSpaHeatMap: check argument: obj! For RGB plot, the second element of list must be a list of matrix which specify the features information!")
  if (any(sapply(obj[[2]], ncol) != 3))
    stop("EachRGBSpaHeatMap: check argument: obj! For RGB plot, the second element of list must be a list of matrix with three columns!")

  for (i in 1:nT) {

    p1 <-.plot_RGB(obj[[1]][[i]], obj[[2]][[i]], no_axis=no_axis,no_axis_name=no_axis_name,
                   axis_names=axis_names,
                               pt_size =pt_size, base_size=base_size, border_col=border_col)
    if (!is.null(title_name)) {
      p1 <- p1 + ggtitle(label = paste0(title_name,
                                        sample_name[i]))
    }
    pList[[i]] <- p1
  }

  if (combine) {
    pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                               nrow = layout.dim[1], ...)
  }
  return(pList)
}

#' @rdname EachRGBSpaHeatMap
#' @method EachRGBSpaHeatMap SRTProject
#' @export
EachRGBSpaHeatMap.SRTProject <- function(obj, batch=NULL, plot_type = c("UMAP", "tSNE"), layout.dim = NULL,
                                     combine=TRUE,  title_name = NULL, no_axis=TRUE, pt_size = 2,  base_size=10, border_col='gray',
                                     axis_names = c("Coord x","Coord y"),
                                     no_axis_name = FALSE, ...){
  #obj <- SRTProj
  plot_type <- match.arg(plot_type)
  ## batch can be a vector representing the number of batches or names of batches

  names_embed3 <- paste0(plot_type, "3")
  posList <- mat2list(obj@spatialCoords, obj@sampleMetaData$NumOfSpots)
  names(posList) <- row.names(obj@sampleMetaData)
  if(!names_embed3 %in% names(obj@plotEmbeddings)){
    stop(paste0("EachRGBSpaHeatMap: Check plot_type! The ", names_embed3, " is not in the `plotEmbeddings` slot!"))
  }
  embed3List <- mat2list(obj@plotEmbeddings[[names_embed3]], obj@sampleMetaData$NumOfSpots)
  samplenames <- row.names(obj@sampleMetaData)
  names(embed3List) <- samplenames
  names(samplenames) <- samplenames
  if(is.null(batch)){
    batch <- samplenames
  }else{
    batch <- samplenames[batch]
  }

  posList <- posList[batch]
  embed3List <- embed3List[batch]

  if(is.null(layout.dim) && length(batch)>1)
    layout.dim <- c(round(length(batch)/2), 2)

  list2 <- list(pos=posList, embed3=embed3List)
  pList <- EachRGBSpaHeatMap(list2, combine=combine,title_name=title_name,
                        no_axis=no_axis,layout.dim=layout.dim,pt_size = pt_size,
                        base_size=base_size, border_col=border_col,
                        axis_names = axis_names,
                        no_axis_name = no_axis_name, ...)
  return(pList)
}


.plot_cont_scatter <- function(embed, feature, quant=0.9, pt_size=2,
                               cols = c("#0571B099",  "#CA0020"),base_size=14,border_col="gray",
                               no_axis=TRUE, no_axis_name=TRUE,legend.position="right",
                               axis_names=c("Coord x","Coord y")){
  require(ggplot2)
  if(!is.matrix(embed))  stop(".plot_cont_scatter: Check argument embed! it must be a matrix!")
  if(ncol(embed)!=2) stop(".plot_cont_scatter: Check argument embed! it must have two colums!")
  dat <- as.data.frame(embed)
  colnames(dat) <- c("row", "col")
  if(!is.vector(feature))  stop(".plot_cont_scatter: Check argument feature! it must be a vector!")
  if(length(cols)!=2) stop(".plot_cont_scatter: Check argument cols! it must have two color names!")

  x <- as.vector(feature)
  dat$Expression <- range01(x )
  med <- quantile(dat$Expression, quant)
  plt <- ggplot(data=dat, aes(x=row , y=col, color=Expression)) + geom_point(size=pt_size) +
    scale_colour_gradient2(low = cols[1], mid = "white", high = cols[2], midpoint = med)+
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col),
          legend.position=legend.position)
  if (no_axis) {
    plt <- plt + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                       axis.title.y = element_blank(), axis.ticks.y = element_blank())
  }
  else if (!no_axis_name) {
    plt <- plt + labs(x = axis_names[1], y = axis_names[2])
  }
  if (legend.position == "bottom") {
    plt <- plt + theme(legend.direction = "horizontal")
  }
  else if (legend.position == "right") {
    plt <- plt + theme(legend.direction = "vertical")
  }
 return(plt)

}
#' Plot spatial feature heatmaps for a group of SRT data batches
#'
#' @description
#' Visualize the gene expressions on spatial coordinates for a group of genes and SRT data batches.
#' @export
#' @param obj a SRTProject object.
#' @param features a string vector, a group of genes to be visualized.
#' @param batch a string or positive integer vector, the names or ID of batches to be plotted
#' @param quantVec a real or real vector between 0 and 1, set the midlle scale of continous color schema.
#' @param combine a logical value, whehter put all plots on one figure, default as TRUE.
#' @param scale a logical value, whether scale the gene expression for plot, default as \code{TRUE}.
#' @param layout.dim a 2-diemensional vector, specify the layout (numbers of rows and columns) of figure.
#' @param no_axis a logical value, whether display the axis in plot.
#' @param pt_size a positive real, the point size in plot.
#' @param  base_size a positive real, the baseline text size in plot.
#' @param border_col a string, the color of border in plot.
#' @param  common.legend a logcial value, whether use a common legend in figure.
#' @param axis_names a 2-dimensional string vector, the axis names in plot.
#' @param no_axis_name a logical value, whether display the axis name in plot.
#' @param legend.position a string, the position of legend.
#' @param ... other arguments pass to \code{\link{ggarrange}}.
#'
#' @return Returns a ggplot object or a list of ggplot objects.
#'
#' @seealso None.
#'
EachExprSpaHeatMap <- function(obj,features, batch=NULL, quantVec=NULL,
                           combine=TRUE, scale=TRUE, ...){
  UseMethod(generic = 'EachExprSpaHeatMap', object = obj)
}
#' @rdname EachExprSpaHeatMap
#' @method EachExprSpaHeatMap list
EachExprSpaHeatMap.list <- function(obj, quantVec = NULL, title_name = TRUE, combine = TRUE,
                                layout.dim = c(round(length(obj[[2]])), ncol(obj[[2]][[1]])), no_axis = TRUE,
                                pt_size = 1,  base_size=10, border_col='gray', common.legend=TRUE,
                                axis_names = c("Coord x","Coord y"),
                                no_axis_name = TRUE, legend.position="right",...){

  require(dplyr)
  if (length(obj) != 2)
    stop("EachExprSpaHeatMap: check argument: obj! The list must constain tow elements. The first is a list of two column matrices to specify the position information, and the second specify the features information!")
  if (!is.list(obj[[1]]))
    stop("EachExprSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")
  if (any(sapply(obj[[1]], ncol) != 2))
    stop("EachExprSpaHeatMap: check argument: obj! The first element of list must be a list of two column matrices which specify the position information!")



  nT = length(obj[[1]])
  sample_name = names(obj[[1]])
  if (is.null(sample_name))
    sample_name = paste0("Batch", c(1:nT))
  pList <- vector("list", nT)

  if (!is.list(obj[[2]]))
    stop("EachExprSpaHeatMap: check argument: obj! For EachExprSpaHeatMap plot, the second element of list must be a list of matrix which specify the features information!")
  nfeatures <- ncol(obj[[2]][[1]])
  featureNames <- colnames(obj[[2]][[1]])
  if (any(sapply(obj[[2]], ncol) != nfeatures))
    stop("EachExprSpaHeatMap: check argument: obj! For EachExprSpaHeatMap plot, the second element of list must be a list of matrix with same columns!")
  if(is.null(quantVec)){
    quantVec <- rep(0.5, nfeatures)
  }

  kk <- 1
  for (i in 1:nT) {
    for(j in 1:nfeatures){
      # p1 <- .plot_cont_scatter(as.matrix(posList[[1]]), feature=dataList[[1]][,1],pt_size=4)
      p1 <-.plot_cont_scatter(obj[[1]][[i]] %>%as.matrix(), feature=obj[[2]][[i]][,j]%>% as.vector(),
                              quant = quantVec[j], no_axis=no_axis,no_axis_name=no_axis_name,
                              axis_names=axis_names,legend.position=legend.position,
                              pt_size =pt_size, base_size=base_size, border_col=border_col, ...)

      if (title_name) {
        p1 <- p1 + ggtitle(label = paste0(featureNames[j])) +
          theme(title = element_text(size=base_size+3, color=1, face='italic'))
      }
      pList[[kk]] <- p1 + theme(legend.title = element_text(size=base_size+1, face="plain"))
      kk <- kk + 1
    }

  }

  if (combine) {
    pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                               nrow = layout.dim[1], common.legend = common.legend,
                               legend = legend.position, ...)
  }
  return(pList)

}

#' @rdname EachExprSpaHeatMap
#' @method EachExprSpaHeatMap SRTProject
#' @export
EachExprSpaHeatMap.SRTProject <- function(obj, features, batch=NULL, quantVec=NULL,
                                      combine=TRUE, scale=TRUE, layout.dim = NULL, no_axis = TRUE,
                                      pt_size = 1,  base_size=10, border_col='gray', common.legend=TRUE,
                                      axis_names = c("Coord x","Coord y"),
                                      no_axis_name = TRUE, legend.position="right",...){


  remain_features <- setdiff(features, row.names(obj@geneMetaData))
  if(length(remain_features)> 0){
    stop(paste0("There are undefined features: ", paste0(remain_features, collapse = ', ')))
  }
  ## Read normalized data from h5 file
  hfile <- obj@projectMetadata$h5filePath
  samplenames <- row.names(obj@sampleMetaData)

  if(is.null(layout.dim)){
    layout.dim <- c(length(samplenames), length(features))
  }


  dataList <- .mylapply(seq_along(samplenames), function(id){

    data <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[features,]
    Matrix::t(data)
  })
  if(length(features)==1){
    dataList <- .mylapply(dataList, function(x) matrix(x, ncol=1))
  }
  if(scale){
    dataList <- .mylapply(dataList, function(x) apply(x, 2, range01))
  }
  names(dataList) <- samplenames
  posList <- mat2list(obj@spatialCoords, obj@sampleMetaData$NumOfSpots)
  names(posList) <- row.names(obj@sampleMetaData)

  if(!is.null(batch)){
    posList <- posList[batch]
    dataList <- dataList[batch]
  }
  obj <- list(pos=posList, data=dataList)
  pList <- EachExprSpaHeatMap(obj, quantVec=quantVec, combine = combine,layout.dim = layout.dim, no_axis = no_axis,
                              pt_size = pt_size,  base_size=base_size, border_col=border_col, common.legend=common.legend,
                              axis_names = c("Coord x","Coord y"),
                              no_axis_name = no_axis_name, legend.position=legend.position,...)
  return(pList)

}


### Plot genes-by-cluster heatmap
doHeatmap <- function (seu, features = NULL,  cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
          grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
          legend.position='right', axis_text_y=TRUE, y_text_size= 5,
          ...){

  require(Seurat)
  if(do.scale){
    seu <- ScaleData(seu, verbose=FALSE)
    seu <- subset(seu, downsample = 400)
  }
  Idents(seu) <- factor(seu@meta.data$cluster)
  ngrp <- nlevels(Idents(seu))
  if (is.null(grp_color)) {
    grp_color <- .gg_color_hue(ngrp)
  }


  if(axis_text_y){
    axis.text.y <- element_text(size = y_text_size, face = "italic",
                                family = "serif")
  }else{
    axis.text.y <- element_blank()
  }

  p1 <- Seurat::DoHeatmap(object = seu, features = features, group.colors = grp_color[1:ngrp],
                    label = grp_label, ...) + guides(color = guide_legend(title = cell_label,
                                                                          ncol = ncol.legend,
                                                                          override.aes = list(stroke = 1, alpha = 1,
                                                                     shape = 16, size = pt_size, color = grp_color[1:ngrp])),
                                                     fill = guide_legend(title=fill_legend_title),
                                                     alpha = "none") + theme(legend.text = element_text(size = 10),
                                                                             legend.title = element_text(size = 13, face = "bold"),
                                                                             axis.text.y = axis.text.y)
  if (legend.position == "bottom") {
    p1 <- p1 + theme(legend.direction = "horizontal")
  }
  else if (legend.position == "right") {
    p1 <- p1 + theme(legend.direction = "vertical")
  }
  return(p1)

}

doHeatmap.spe <- function (spe, features,  cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                           grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                           legend.position='right', axis_text_y=TRUE, y_text_size= 5,
                           ...){

  require(Seurat)
  require(Matrix)
  require(SpatialExperiment)

  if(is.element("logcounts",assayNames(spe))){
    logcount <- assay(spe, "logcounts")[features,]
    count <- sparseMatrix(i=1,j=1, dims=dim(logcount))
    row.names(count) <- row.names(logcount)
    colnames(count) <- colnames(logcount)
    seu <- CreateSeuratObject(counts=count, meta.data= as.data.frame(colData(spe)))
    asy <- DefaultAssay(seu)
    seu[[asy]]@data <- logcount
  }else{
    count <- assay(spe, "counts")[features,]
    seu <- CreateSeuratObject(counts=count, meta.data= as.data.frame(colData(spe)))
    seu <- NormalizeData(seu, verbose = FALSE)
  }

  if(do.scale){
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- subset(seu, downsample = 400)
  }
  Idents(seu) <- factor(seu@meta.data$clusters)
  ngrp <- nlevels(Idents(seu))
  if (is.null(grp_color)) {
    grp_color <- .gg_color_hue(ngrp)
  }


  if(axis_text_y){
    axis.text.y <- element_text(size = y_text_size, face = "italic",
                                family = "serif")
  }else{
    axis.text.y <- element_blank()
  }

  p1 <- Seurat::DoHeatmap(object = seu, features = features, group.colors = grp_color[1:ngrp],
                          label = grp_label, ...) + guides(color = guide_legend(title = cell_label,
                                                                                ncol = ncol.legend,
                                                                                override.aes = list(stroke = 1, alpha = 1,
                                                                                                    shape = 16, size = pt_size, color = grp_color[1:ngrp])),
                                                           fill = guide_legend(title=fill_legend_title),
                                                           alpha = "none") + theme(legend.text = element_text(size = 10),
                                                                                   legend.title = element_text(size = 13, face = "bold"),
                                                                                   axis.text.y = axis.text.y)
  if (legend.position == "bottom") {
    p1 <- p1 + theme(legend.direction = "horizontal")
  }
  else if (legend.position == "right") {
    p1 <- p1 + theme(legend.direction = "vertical")
  }
  return(p1)

}
#' Plot gene-by-cell heatmaps for a group of SRT data batches
#'
#' @description
#' Visualize the expressions of the selected features (genes) for all cells/spots.
#' @export
#' @param obj a SRTProject object or a SpatialExperiment object.
#' @param features a string vector, a group of selected genes' names.
#' @param cell_label a string, specify name of clusters.
#' @param fill_legend_title a string, the legend title.
#' @param do.scale a logical value, whether scale the gene expression for plot, default as \code{TRUE}.
#' @param grp_label a logical value, whether display the cluster labels in plot.
#' @param pt_size a positive real, the point size in plot.
#' @param grp_color a string vector, the colors for each cluster in plot.
#' @param ncol.legend a positive number, the number of columns for legend.
#' @param legend.position a string, the postion of legend.
#' @param axis_text_y a logical value, whether display the tick text in y-axis.
#' @param y_text_size a positive real, the text size for y-axis.
#' @param combine a logical value, whether put all plots on one figure, default as \code{TRUE}.
#' @param ... other arguments pass to \code{ggarrange}.
#' @return Returns a ggplot object or a list of ggplot objects.
#'
#' @importFrom ggpubr ggarrange
#' @seealso None
#'
EachGCHeatMap <- function(obj,features, batch=NULL,  layout.dim= NULL,
                          common.legend=FALSE,  cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                          grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                          legend.position='right', axis_text_y=TRUE, y_text_size= 5,combine=TRUE, ...){
  UseMethod(generic = 'EachGCHeatMap', object = obj)
}

#' @rdname EachGCHeatMap
#' @method EachGCHeatMap SRTProject
#' @export
EachGCHeatMap.SRTProject <- function(obj, features,  batch=NULL,  layout.dim= NULL,
                                     common.legend=FALSE,  cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                                     grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                                     legend.position='right', axis_text_y=TRUE, y_text_size= 5,combine=TRUE, ...){

  # features <- paste0("gene",2:10)
  ## Get normalized features*spots data
  hfile <- obj@projectMetadata$h5filePath
  samplenames <- .getDataFromH5file(hfile, groupname = 'samplenames')


  names(samplenames) <- samplenames
  if(is.null(batch)){
    batch <- samplenames
  }
  if(length(obj@clusters)==0){
    stop("EachGCHeatMap: there is no cluster labels in `clusters` slot!")
  }

  if(!is.null(batch)){

    ids <- match(batch, samplenames) # select the samples' ids.

    if(is.null(layout.dim))
      layout.dim <- c(round(length(batch)/2), 2)

    seuList <- .mylapply(ids, function(id){
      dat1 <-  .getSparseMatrixFromH5file(hfile, groupname = paste0('data/data_', id))[features,]
      count <- sparseMatrix(i=1,j=1, x=0, dims=dim(dat1))
      row.names(count) <- row.names(dat1)
      colnames(count) <- colnames(dat1)
      meta_data_here <- data.frame(cluster=obj@clusters[obj@cellMetaData$batch==batch[id]])
      row.names(meta_data_here) <- colnames(count)
      seuInt <- CreateSeuratObject(counts = count, meta.data = meta_data_here)
      seuInt[['RNA']]@data <- as.sparse(dat1)
      Idents(seuInt) <- factor(meta_data_here$cluster)
      return(seuInt)
    })


   pList <- .mylapply(seuList, doHeatmap, features=features,  cell_label = cell_label, fill_legend_title=fill_legend_title, do.scale=do.scale,
                      grp_label = grp_label, pt_size = pt_size, grp_color = grp_color, ncol.legend = ncol.legend,
                      legend.position=legend.position, axis_text_y=axis_text_y, y_text_size=y_text_size,...)
   # pList <- .mylapply(seuList, doHeatmap, features=features)
   # doHeatmap(seuList[[1]], features=features)
   if(length(seuList)==1){
     return(pList[[1]])
   }
   if (combine && length(seuList)>1) {
     pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                                nrow = layout.dim[1], common.legend = common.legend,
                                legend = legend.position, ...)

   }
   return(pList)
  }


}


#' @rdname EachGCHeatMap
#' @method EachGCHeatMap SpatialExperiment
#' @export
EachGCHeatMap.SpatialExperiment <- function(obj, features,  batch=NULL,  layout.dim= NULL,
                                     common.legend=FALSE,  cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                                     grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                                     legend.position='right', axis_text_y=TRUE, y_text_size= 5,combine=TRUE, ...){

  # features <- paste0("gene",2:10)
  ## Get normalized features*spots data
  require(SpatialExperiment)
  samplenames <- unique(colData(obj)[,'sample_id'])
  names(samplenames) <- samplenames


  if(!is.null(batch)){
    samplenames <- samplenames[batch]
  }
  if(is.null(layout.dim) && length(samplenames)>1){
    layout.dim <- c(2, round(length(samplenames)/2) )
  }
  if(is.null(layout.dim) && length(samplenames) == 1){
    layout.dim <- c(1,1)
  }
  pList <- .mylapply(samplenames, function(id){

    ret_spe <- obj[, colData(obj)[,'sample_id'] == id]
    p <- doHeatmap.spe(ret_spe,features = features,  cell_label = cell_label, fill_legend_title=fill_legend_title, do.scale=do.scale,
                       grp_label = grp_label, pt_size = pt_size, grp_color = grp_color, ncol.legend = ncol.legend,
                       legend.position=legend.position, axis_text_y=axis_text_y, y_text_size=y_text_size, ...)

    return(p)
  })


  if (combine) {
    pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                               nrow = layout.dim[1], common.legend = common.legend,
                               legend = legend.position, ...)
  }
  return(pList)


}


#' Plot gene-by-cell heatmap
#'
#' @description
#' Visualize the expressions of the selected features (genes) for all cells/spots.
#' @export
#' @param obj a SpatialExperiment object.
#' @param features a string vector, a group of selected genes' names.
#' @param cell_label a string, specify name of clusters.
#' @param fill_legend_title a string, the legend title.
#' @param do.scale a logical value, whether scale the gene expression for plot, default as \code{TRUE}.
#' @param grp_label  a logical value, whether display the cluster labels in plot.
#' @param pt_size a positive real, the point size in plot.
#' @param grp_color a string vector, the colors for each cluster in plot.
#' @param ncol.legend a positive number, the number of columns for legend.
#' @param legend.position a string, the postion of legend.
#' @param axis_text_y  a logical value, whether display the tick text in y-axis.
#' @param y_text_size a positive real, the text size for y-axis.
#' @param ... other arguments pass to \code{DoHeatmap}.
#'
#' @return Returns a ggplot object.
#'
#' @seealso None
#'
GCHeatMap <- function(obj, features,cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                                        grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                                        legend.position='right', axis_text_y=TRUE, y_text_size= 5, ...){
  UseMethod(generic = 'GCHeatMap', object = obj)
}

#' @return Returns a ggplot object.
#' @rdname GCHeatMap
#' @method GCHeatMap SpatialExperiment
#' @export
GCHeatMap.SpatialExperiment <- function(obj, features,cell_label = "Domain", fill_legend_title='Expression', do.scale=TRUE,
                                            grp_label = FALSE, pt_size = 4, grp_color = NULL, ncol.legend = 1,
                                            legend.position='right', axis_text_y=TRUE, y_text_size= 5, ...){

  # features <- paste0("gene",2:10)
  ## Get normalized features*spots data
  require(SpatialExperiment)


  p <- doHeatmap.spe(obj, features = features,  cell_label = cell_label, fill_legend_title=fill_legend_title, do.scale=do.scale,
                       grp_label = grp_label, pt_size = pt_size, grp_color = grp_color, ncol.legend = ncol.legend,
                       legend.position=legend.position, axis_text_y=axis_text_y, y_text_size=y_text_size, ...)

  return(p)



}



###
#' Plot cluster-by-cluster correlation hetamap
#'
#' @description
#' Visualize the cluster-by-cluster correlation based on a specified low-dimensional embeddings.
#' @export
#' @param obj A SRTProject object with low dimension embeddings.
#' @param reduction a string, the name of reduction (low-dimensional embeddings).
#' @param legend_title a string, the title of legend.
#' @param grp_label a logical value, whether display the cluster labels in plot.
#' @param pt_size a positive real, the point size in plot.
#' @param grp_color  a string vector, the colors for each cluster in plot.
#' @param fill_legend_title a string, the legend title for filled colors.
#' @param seed a positive integer, a ranodm seed.
#' @param ... other arguments pass to \code{DoHeatmap}.
#'
#' @return Returns a ggplot object.
#'
#' @seealso None.
#'
CCHeatMap <- function(obj,reduction=NULL,legend_title='Domain', grp_label = FALSE,
                      pt_size=4, grp_color=NULL,fill_legend_title="Pearson's correlation", seed=1,  ...){
  UseMethod(generic = 'CCHeatMap', object = obj)
}
#' @rdname CCHeatMap
#' @method CCHeatMap SRTProject
#' @export
CCHeatMap.SRTProject <- function(obj, reduction=NULL,legend_title='Domain', grp_label = FALSE,
                                 pt_size=4, grp_color=NULL,fill_legend_title="Pearson's correlation", seed=1, ...){

  if(is.null(reduction)){
    reduction <- tail(names(obj@reductions), 1)
  }

  if(!reduction %in% names(obj@reductions)){
    stop("CCHeatMap: check argument `reduction`, there is no ", reduction, " in `reductions` slot!")
  }
  embeds <- SRTProj@reductions[[reduction]]
  clusters <- SRTProj@clusters
  set.seed(1)
  subsample <- function(hZ_idrsc, cluster_idrsc_vec, nn_sub, seed){
    set.seed(seed)
    nn <- length(cluster_idrsc_vec)
    idx <- sort(sample(nn, min(nn, nn_sub)))
    return(list(hZ=hZ_idrsc[idx, ], cluster=cluster_idrsc_vec[idx]))
  }
  sublist <- subsample(embeds, clusters, 400, seed=seed)
  order_idx <- order(sublist[[2]])
  cluster_orderd <- sublist[[2]][order_idx]
  corMat <- cor(t(sublist[[1]][order_idx, ]))
  row.names(corMat) <- paste0('spot', 1:nrow(corMat))
  colnames(corMat) <- paste0('spot', 1:nrow(corMat))
  p1 <- doHeatmap.matrix(corMat, cluster_orderd,legend_title='Domain', grp_label = grp_label,
                         pt_size=pt_size, grp_color=grp_color,fill_legend_title=fill_legend_title, ...)
  return(p1)

}

doHeatmap.matrix <- function(corMat, cluster_orderd,legend_title='Domain', grp_label = FALSE,
                             pt_size=4, grp_color=NULL,fill_legend_title="Pearson's correlation", ...){

  require(Seurat)
  meta_data <- data.frame(cluster=factor(cluster_orderd))
  row.names(meta_data) <- colnames(corMat)
  seu <- CreateSeuratObject(counts=corMat,
                            meta.data = meta_data)

  seu[["RNA"]]@scale.data <- corMat

  doHeatmap(seu, features = row.names(seu), do.scale = FALSE, cell_label=legend_title, pt_size=pt_size,
            grp_color=grp_color, grp_label=grp_label, axis_text_y=FALSE,fill_legend_title=fill_legend_title, ...)

}

##
#' Two-dimensional embeddings plot
#'
#' @description
#' Visualize the specified item on the two-dimensional embeddings such as tSNE, UMAP and so on.
#' @export
#' @param obj a SRTProject object or SpatialExperiment object.
#' @param item a string, the item used to mapped as colors.
#' @param plotEmbeddings a string, the name of plotembeddings, sucha as tSNE and UMAP.
#' @param axis_names a 2-dimensional string vector, the names of axis.
#' @param cols a string vector, the selected colors.
#' @param no_guides a logical value, whether add guides.
#' @param no_axis_name a logical value, whether display the axis names.
#' @param pt_size a positive real, the point size in plot.
#' @param point_alpha a positive real between 0 and 1, the point transparency in plot.
#' @param nrow.legend a positive number, the number of columns for legend.
#' @param base_size  a positive real, the baseline text size in plot.
#' @param border_col a string, the color of border in plot.
#' @param legend.position a string, the positive of legend.
#' @param no_axis a logical value, whether display the axis.
#' @param do_density a logical value, whether plot the density plot.
#' @param add_border_box a logical value, whether add border box in plot.
#'
#' @return returns a ggplot object.
#'
#' @seealso \code{\link{SRTProject-class}}.
#'




EmbedPlot <- function(obj, item="cluster", plotEmbeddings='tSNE',  ...){
  UseMethod(generic = 'EmbedPlot', object = obj)
}

#' @rdname EmbedPlot
#' @method EmbedPlot SRTProject
#' @export
EmbedPlot.SRTProject <- function(obj, item="cluster", plotEmbeddings=NULL,
                                 axis_names = NULL,
                                 cols = NULL, no_guides = FALSE,
                                 no_axis_name = FALSE, pt_size = 1, point_alpha = 1,
                                 nrow.legend = NULL, base_size=10,border_col='gray',
                                 legend.position = "bottom", no_axis = FALSE, do_density = FALSE){


  if(is.null(axis_names)){
    axis_names <- paste0(plotEmbeddings, 1:2)
  }
  if(is.null(plotEmbeddings)){
    plotEmbeddings <- tail(names(obj@plotEmbeddings), 1)
  }
  pEmbedName <-plotEmbeddings
  if(!pEmbedName %in% names(obj@plotEmbeddings)){
    stop("EmbedPlot: check argument `model` or `plotEmbeddings`, there is no ", pEmbedName, " in `plotEmbeddings` slot!")
  }

  plotEmbed <- obj@plotEmbeddings[[pEmbedName]]


  if(tolower(item)== 'cluster'){
    p1 <- .plot_scatter(plotEmbed, labels = factor(obj@clusters), label_name = "Cluster",
                        axis_names = axis_names,
                        cols = cols, no_guides = no_guides,
                        no_axis_name = no_axis_name, pt_size = pt_size, point_alpha = point_alpha,
                        nrow.legend = nrow.legend, base_size=base_size,border_col=border_col,
                        legend.position = legend.position, no_axis = no_axis, do_density = do_density)

  }else if(tolower(item)=='batch'){
    p1 <- .plot_scatter(plotEmbed, labels = factor(obj@cellMetaData$batch), label_name = "Batch",
                        axis_names = axis_names,
                        cols = cols, no_guides = no_guides,
                        no_axis_name = no_axis_name, pt_size = pt_size, point_alpha = point_alpha,
                        nrow.legend = nrow.legend, base_size=base_size,border_col=border_col,
                        legend.position = legend.position, no_axis = no_axis, do_density = do_density)
  }else{
     stop("EmbedPlot.SRTProject: Check argument: item! Unsupported item name!")
  }

  return(p1)
}

#' @rdname EmbedPlot
#' @method EmbedPlot SpatialExperiment
#' @export
EmbedPlot.SpatialExperiment <- function(spe, plotEmbeddings=NULL,
                                            colour_by = NULL,
                                            base_size = 15,
                                            legend.position='right', add_border_box=FALSE, ...){
  require(scater)
  if(is.null(plotEmbeddings)){
    plotEmbeddings <- tail(reducedDimNames(spe), 1)
  }
  if(is.null(colour_by)){
    colour_by <-  'clusters'
  }
  p <- plotReducedDim(spe, dimred = plotEmbeddings, colour_by = colour_by,theme_size=base_size, ...)

  if(add_border_box){
   p <- p +theme(panel.background= element_rect(fill = 'white', color='black'))
  }
  return(p)

}

# Downstream plots --------------------------------------------------------
#' Barplot for the enrichment terms
#'
#' @param top_dat a data frame with the information of enrichment terms
#' @param source a string, the source of database.
#' @param term_name a string, the name of y-axis.
#' @param nlog10P a string, the colname in top_dat that represents the negative log10 p-values.
#' @param bar_width a positive real, the width of bars.
#' @param base_size  a positive real, the  baseline text size in plot.
#' @param font_family a string, the font.
#' @param cols a string vector, the used colors in plot.
#' @return return a ggplot object.
#' @export
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


#' Visualize the inferred pseudotime on  embeddings.
#'
#' @param SRTProj a SRTProject object with cellMetaData Traject.commonPseudotime.
#' @param plotEmbeddings a string, specify the embeddings.
#' @param no_axis  a logical value, whether display the axis in plot.
#' @param base_size a positive real, the baseline text size in plot.
#' @param border_col a string, the color of border in plot.
#' @param legend.position a string, the position of legend.
#' @export
#'
TrajectEmbedPlot <- function(SRTProj, plotEmbeddings = "tSNE",
                             no_axis= FALSE,base_size=12,border_col='gray',
                             legend.position='right'){
  require(ggplot2)


  pEmbedName <- plotEmbeddings
  if(!pEmbedName %in% names(SRTProj@plotEmbeddings)){
    stop("EmbedPlot: check argument  `plotEmbeddings`, there is no ", pEmbedName, " in `plotEmbeddings` slot!")
  }
  plotEmbed <- SRTProj@plotEmbeddings[[pEmbedName]]
  pseudotime <-SRTProj@cellMetaData$Traject.commonPseudotime
  tmp_dat <- data.frame(embed1=plotEmbed[,1], embed2=plotEmbed[,2], pseudotime=pseudotime)
  plt <- ggplot(tmp_dat, aes(x=embed1, y=embed2, color=pseudotime)) + geom_point(size=1.5) +
    scale_color_gradientn(colours = c( "#FFBB78", "#96B5E9", "#F54902"))+
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.position = legend.position,
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col)) +
    xlab(paste0(pEmbedName, ' 1')) + ylab(paste0(pEmbedName, ' 2'))

  if (no_axis) {
    plt <- plt + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                       axis.title.y = element_blank(), axis.ticks.y = element_blank())
  }

  return(plt)
}

##
#'  Plot the embeddings for each SRT data
#'
#' @param obj an object with class \code{SpatialExperiment }.
#' @param reduction a string, speficify which reduction slot is used to plot. If \code{reduction=NULL}, the last one is used.
#' @param batch a string vector, denoting the names of batches to be plotted, or a integer vector denoting the ID numbers of batches to be plotted.
#' @param colour_by a string, specify the colname in cell metadata to be mapped to colors.
#' @param combine a logical value, whether put all plots in one figure.
#' @param base_size a positive real, the baseline text size.
#' @param layout.dim a 2-dimension integer vector, specify the layout of figure.
#' @param common.legend a logical value, whether use a common legend for all plots.
#' @param legend.position a string, the position of legend.
#' @param add_border_box a logical value, whether add the border box.
#' @param no_axis a logical value, whether display the axis.
#' @param ... other arguments that pass to other function.
#' @return return a ggplot object of a list of ggplot objects.
#'
#' @export
EachEmbedPlot <- function(obj, reduction=NULL, batch=NULL,  ...){
  UseMethod(generic = 'EachEmbedPlot', object = obj)
}

#' @rdname EachEmbedPlot
#' @method EachEmbedPlot SpatialExperiment
#' @export
EachEmbedPlot.SpatialExperiment <- function(obj, reduction=NULL,
                                            colour_by = NULL,batch=NULL, combine=TRUE,
                          base_size = 15, layout.dim = NULL, common.legend=TRUE,
                          legend.position='right', add_border_box=FALSE, no_axis=TRUE, ...){
  require(scater)
  if(is.null(reduction)){
    reduction <- tail(reducedDimNames(obj), 1)
  }
  if(is.null(colour_by)){
    colour_by <-  'clusters'
  }
  samplenames <- unique(colData(obj)[,'sample_id'])
  names(samplenames) <- samplenames
  if(!is.null(batch)){
    samplenames <- samplenames[batch]
  }
  if(is.null(layout.dim) && length(samplenames)>1){
    layout.dim <- c(2, round(length(samplenames)/2) )
  }
  if(is.null(layout.dim) && length(samplenames) == 1){
    layout.dim <- c(1,1)
  }
  pList <- lapply(samplenames, function(id){

    ret_obj <- obj[, colData(obj)[,'sample_id'] == id]
    p <- plotReducedDim(ret_obj, dimred = reduction, colour_by = colour_by,theme_size=base_size, ...)
    if(add_border_box){
      p <- p +  theme(panel.background= element_rect(fill = 'white', color= "black"))
    }
    if (no_axis) {
      p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                         axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                         axis.title.y = element_blank(), axis.ticks.y = element_blank())
    }


    return(p)
  })


  if (combine) {
    pList <- ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                               nrow = layout.dim[1], common.legend = common.legend,
                               legend = legend.position, ...)
  }
  return(pList)
}




## Others to use the Seurat Plot as interface


##
#' Draw a figure using a group of ggplot objects
#'
#' @param pList a list with component ggplot objects.
#' @param layout.dim a integer vector with length 2, the layout of subplots in rows and columns.
#' @param common.legend a logical value, whether use common legend for all subplots.
#' @param legend.position a string, the position of legend.
#' @param ... other arguments that pass to \code{\link{ggarrange}}.
#' @return return a new ggplot object.
#' @export
drawFigs <- function(pList, layout.dim = NULL, common.legend=FALSE,legend.position='right',  ...){
  if(!is.list(pList)) stop('drawFigs: pList must be a list!')

  if(is.null(layout.dim) && length(pList)>1){
    layout.dim <- c(2, round(length(pList)/2) )
  }
  if(is.null(layout.dim) && length(pList) == 1){
    layout.dim <- c(1,1)
  }
  ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                             nrow = layout.dim[1], common.legend = common.legend,
                             legend = legend.position, ...)

}

#' Save a ggplot object to a PNG figure
#'
#' @param plt a ggplot object, the object to be saved.
#' @param filename a string, the file names of the saved figure, default as `myfigure`.
#' @param y_reverse a logical value, whether reverse the graph by y-axis, default as \code{FALSE}.
#' @param x_reverse a logical value, whether reverse the graph by x-axis, default as \code{FALSE}.
#' @param width a positive real, the width of the saved figure, default as 7.
#' @param height a positive real, the height of the saved figure, default as 5.5.
#' @param dpi a positive integer, the resolution of the figure, default as 200.
#'
#' @export
write_fig <- function(plt, filename='myfigure', y_reverse=FALSE,
                      x_reverse=FALSE,width =7, height =5.5,  dpi=200){

  require(ggplot2)
  if(y_reverse)
    plt <- plt + scale_y_reverse()
  if(x_reverse)
    plt <- plt + scale_x_reverse()
  if(!is.null(filename)){
    if(!dir.exists("Figs")){
      dir.create("Figs")
    }
    ggsave(file=paste0("./Figs/",filename,".png"), plot = plt,
           width =width, height =height, units = "in", dpi = dpi)
  }
}
