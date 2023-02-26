#
# rm(list=ls())
# library(SRTpipeline)
# setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\")
# source(normalizePath("./SRTpipeline/R/setClass.R"))
# source(normalizePath("./SRTpipeline/R/main.R"))
# source(normalizePath("./SRTpipeline/R/embeddings.R"))
# source(normalizePath("./SRTpipeline/R/visualization.R"))
# source(normalizePath("./SRTpipeline/R/utility.R"))
# source(normalizePath("./SRTpipeline/R/model_newinterface.R"))
# source(normalizePath("./SRTpipeline/R/downstream.R"))

# setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\SRTpipeline")
# setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\SRTpipeline\\vignettes\\")
#
# setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\LearnH5file\\")
# load("mylast_data.rds")
# load("DR_SC_mylast_data.rds")
# library(Seurat)
# library(ggplot2)
# library(S4Vectors)
# library(dplyr)
# library(rhdf5)
# library(SRTpipeline)
# ## inputFiles = NULL, sampleNames = names(inputFiles), platform = "Visium"
# seuList <- PRECAST:::gendata_seulist()
# x <- seuList[[1]][['RNA']]@counts
#
#
# # Part I: Single sample analysis#### --------------------------------------------------
# #  -----------------------------------------------------
# library(S4Vectors)
# cntList <- SimpleList(x1=x)
# colnames(cntList[[1]])
# coordList <- lapply(1:1, function(x)  cbind(seuList[[1]]$row, seuList[[1]]$col))
# metadataList <- lapply(cntList, function(x) as.data.frame(t(x[1:2,])))
# filename <- "simuDRSC"
# sampleMetadata <- data.frame(sex=rep("Female", 1), age = c(24))
# row.names(sampleMetadata) <- names(cntList)
#
#
#
# SRTProj <- CreateSRTProject(cntList, coordList, projectName = "SimuSRTProject", metadataList,
#                             sampleMetadata, min.genes = 20, min.spots = 50, force = T)
#
# ## Check the data saved in the h5 file on the disk
# h5ls(SRTProj@projectMetadata$h5filePath)
#
#
# SRTProj <- normalizeSRT(SRTProj, force = F)
#
# ## Select top n variable genes
#
#
# SRTProj <- selectVariableFeatures(SRTProj, nfeatures = 50)
#
# ## Obtain adjacence matrix
#
# SRTProj <- AddAdj(SRTProj, platform = "ST")
#
# # SC.MEB fitting ----------------------------------------------------------
# SRTProj <- AddPCA(SRTProj)
# SRTProj <- SRTfitSCMEB(SRTProj, K= 4)
# SRTProj <- AddTSNE(SRTProj, n_comp = 2)
# SRTProj <- AddUMAP(SRTProj, n_comp = 3)
# # names(SRTProj@plotEmbeddings) <- c("tSNE", "UMAP3")
# EachRGBSpaHeatMap(SRTProj,plot_type = "UMAP", pt_size=4 ,title_name="UMAP RGB plot: ")
#
# EmbedPlot(SRTProj, item = "cluster")
# EmbedPlot(SRTProj, item = "batch")
#
# EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"), common.legend = F,
#           legend.position = 'right', y_text_size=20)
# CCHeatMap(SRTProj)
#
#
#
# # -------------------- DR-SC model fitting -------------------------------
#
# SRTProj <- SRTfitDRSC(SRTProj, K=4)
# SRTProj@reductions
# SRTProj <- AddTSNE(SRTProj, n_comp = 2)
#
# ## Plot for specified batches
# EachClusterSpaHeatMap(SRTProj, title_name="DR-SC: ")
#
# SRTProj <- AddTSNE(SRTProj, n_comp = 3)
# EachEachRGBSpaHeatMap(SRTProj, batch=1,plot_type = "tSNE", pt_size=4 ,title_name="tSNE RGB plot: ")
# EachEachRGBSpaHeatMap(SRTProj, batch=1:3,plot_type = "tSNE", pt_size=4 ,title_name="tSNE RGB plot: ")
#
# SRTProj <- AddUMAP(SRTProj, n_comp = 3)
# EachRGBSpaHeatMap(SRTProj, batch=1,plot_type = "UMAP", pt_size=4,
#               title_name="UMAP RGB plot: ")
#
# EachExprSpaHeatMap(SRTProj, features=c("gene2", "gene3"), quantVec=c(0.2,0.2), title_name=T)
#
# p12 <- EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"))
# EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"), common.legend = F,
#           legend.position = 'right', y_text_size=20)
#
#
#
# CCHeatMap(SRTProj, reduction = "PCA")
#
#
# ##
# EmbedPlot(SRTProj, item = "Cluster")
# EmbedPlot(SRTProj, item = "Cluster", legend.position='right')
# EmbedPlot(SRTProj, item = "Cluster", axis_names=c("t1", "t2"))
#
# save.image(file='DR_SC_mylast_data.rds')
#
# objectnames <- ls()
# flag_funcs <- sapply(objectnames, function(x) is.function(get(x)))
# rm(list=objectnames[flag_funcs])
#
#
#
# # SpatialAnno fitting -----------------------------------------------------
#
#
# SparseMatWrite(count_matrix, groupname='counts', hfile='tmp_count2')
# tmp_count <- .getSparseMatrixFromH5file(hfile='tmp_count2')
#
# tmp_count <- .getSparseMatrixFromH5file(hfile='tmp_count2')
#
#
#
# # Get started with SRTpipeline: --------------------------------------------
# # SRTpipeline is designed to xxx
# # create SRTProject object ------------------------------------------------
# library(Seurat)
# library(ggplot2)
# library(S4Vectors)
# library(dplyr)
# library(rhdf5)
# library(SRTpipeline)
# ## inputFiles = NULL, sampleNames = names(inputFiles), platform = "Visium"
# seuList <- PRECAST:::gendata_seulist()
# x <- seuList[[1]][['RNA']]@counts
#
# library(S4Vectors)
# cntList <- SimpleList(x1=x)
# colnames(cntList[[1]])
# coordList <- lapply(1:1, function(x)  cbind(seuList[[1]]$row, seuList[[1]]$col))
# metadataList <- lapply(cntList, function(x) as.data.frame(t(x[1:2,])))
# filename <- "simuDRSC"
# sampleMetadata <- data.frame(sex=rep("Female", 1), age = c(24))
# row.names(sampleMetadata) <- names(cntList)
#
#
#
# SRTProj <- CreateSRTProject(cntList, coordList, projectName = "SimuSRTProject", metadataList,
#                             sampleMetadata, min.genes = 20, min.spots = 50, force = T)
#
# ## Check the data saved in the h5 file on the disk
# h5ls(SRTProj@projectMetadata$h5filePath)
#
# # Normalization
# SRTProj <- normalizeSRT(SRTProj, force = F)
# ## Select top n variable genes
# SRTProj <- selectVariableFeatures(SRTProj, nfeatures = 50)
# ## Obtain adjacence matrix
# SRTProj <- AddAdj(SRTProj, platform = "ST")
#
# # Dimension reduction -----------------------------------------------------
# ## approximated PCA
# SRTProj <- AddPCA(SRTProj)
#
# ## accurate PCA
# SRTProj <- AddPCA(SRTProj, Method='PCA')
# ## weighted PCA
# SRTProj <- AddPCA(SRTProj, Method='WPCA')
#
# # clustering --------------------------------------------------------------
# SRTProj <- SRTfitSCMEB(SRTProj, K= 4)
#
# # Joint dimension reduction and clustering --------------------------------
# SRTProj <- SRTfitDRSC(SRTProj, K=4)
#
# # Annotation --------------------------------------------------------------
# ## SpatialAnno
#
#
#
# # plot Embeddings ----------------------------------------------------------
# # To run tSNE in SRTpipeline we use the AddTSNE() function:
# SRTProj <- AddTSNE(SRTProj, n_comp = 2, reduction='PCA')
# # The `reduction` tells AddTSNE uses specified reduction in the SRTProj@reductions slot.
#
# # To plot the 2-dimensional tSNE results, we use the plotEmbedding() function and pass the name of the UMAP embedding we just generated (“UMAP”). We can tell SRTpipeline how to color the spots by using a combination of cols which tells ArchR which matrix to use to find the specified metadata column provided to name.
# p_pca_tsne2 <- EmbedPlot(SRTProj, item='cluster', plotEmbeddings = 'tSNE', cols=chooseColors(n_colors = 4))
# p_pca_tsne2
# # Instead of coloring by “Sample” as above, we can color by “Clusters” which were identified in a previous chapter.
# SRTProj <- AddTSNE(SRTProj, n_comp = 3, reduction='DR.SC')
# p_drsc_tsne3<- EachRGBSpaHeatMap(SRTProj, plot_type = "tSNE", pt_size=4 ,title_name="tSNE RGB plot: ")
# p_drsc_tsne3
#
# # To run UMAP in SRTpipeline we use the AddUMAP() function:
# SRTProj <- AddUMAP(SRTProj, n_comp = 2, reduction='PCA')
# p_pca_umap2 <- EmbedPlot(SRTProj, item='cluster', plotEmbeddings = 'UMAP', cols=chooseColors(n_colors = 4))
# p_pca_umap2
# SRTProj <- AddUMAP(SRTProj, n_comp = 3, reduction='DR.SC')
# p_drsc_umap3<- EachRGBSpaHeatMap(SRTProj, plot_type = "UMAP", pt_size=4 ,title_name="UMAP RGB plot: ")
# p_drsc_umap3
# library(patchwork)
# p_pca_tsne2 + p_drsc_tsne3 + p_pca_umap2 + p_drsc_umap3 + plot_layout(nrow=2, ncol=2, byrow=F)
#
# ## to save the plot, we can use write_fig()
#
# # Visualization -----------------------------------------------------------
#
# EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"), common.legend = F,
#           legend.position = 'right', y_text_size=20)
#
# CCHeatMap(SRTProj)
#
#
# # DEG analysis ------------------------------------------------------------
#
#
#
# # Trajectory inference ----------------------------------------------------
#
#
#
#
# # Part II: Multiple sample analysis##### -------------------------------------------
#
# # create SRTProject object ------------------------------------------------
# library(Seurat)
# library(ggplot2)
# library(S4Vectors)
# library(dplyr)
# library(rhdf5)
# library(SRTpipeline)
# ## inputFiles = NULL, sampleNames = names(inputFiles), platform = "Visium"
# seuList <- PRECAST:::gendata_seulist()
# x <- seuList[[1]][['RNA']]@counts
# library(S4Vectors)
# cntList <- SimpleList(x1=x,x2=x, x3=x)
# colnames(cntList[[1]])
# coordList <- lapply(1:3, function(x)  cbind(seuList[[1]]$row, seuList[[1]]$col))
# metadataList <- lapply(cntList, function(x) as.data.frame(t(x[1:2,])))
# sampleMetadata <- data.frame(sex=rep("Female", 3), age = c(24, 50, 60))
# row.names(sampleMetadata) <- names(cntList)
#
#
# SRTProj <- CreateSRTProject(cntList, coordList, projectName = "SimuMultiSRTProject", metadataList,
#                             sampleMetadata, min.genes = 20, min.spots = 50, force = T)
#
# ## Check the data saved in the h5 file on the disk
# h5ls(SRTProj@projectMetadata$h5filePath)
#
#
# # Data preprocessing ------------------------------------------------------
# ## Data normalization
# SRTProj <- normalizeSRT(SRTProj, force = F)
#
# ## Select top n variable genes
# setSRTverbose(verbose=T)
# SRTProj <- selectVariableFeatures(SRTProj, nfeatures = 50)
# getSRTVerbose()
# ## Obtain adjacence matrix
# setSRToutputPrefix(SRTprefix = "====")
# SRTProj <- AddAdj(SRTProj, platform = "ST",  force = T)
#
# # Dimension reduction -----------------------------------------------------
# ## approximated PCA
# SRTProj <- AddPCA(SRTProj)
#
# ## accurate PCA
# SRTProj <- AddPCA(SRTProj, Method='PCA')
# ## weighted PCA
# SRTProj <- AddPCA(SRTProj, Method='WPCA')
#
#
# # clustering --------------------------------------------------------------
# ## Fit iSC.MEB
# SRTProj <- Integrate_iSCMEB(SRTProj, K=4, reduction = 'PCA')
#
#
# # Joint dimension reduction and clustering --------------------------------
#
# ## Fit PRECAST
# SRTProj <- Integrate_PRECAST(SRTProj, K=4)
#
#
# # plotEmbeddings ----------------------------------------------------------
# SRTProj <- AddTSNE(SRTProj, n_comp = 2, reduction = 'aligned.iSC.MEB')
#
# cbp <- chooseColors("Nature 10", n_colors = 4)
# p_all <- EachClusterSpaHeatMap(SRTProj, layout.dim=c(2,2), cols=cbp,
#                            no_axis=T, title_name='Batch: ', point_size=4,base_size=20 ,
#                            nrow.legend=1, no_guides=F)
# ## Plot for specified batches
# p12 <- EachClusterSpaHeatMap(SRTProj,batch=c("x1", "x2"), title_name="", cols=cbp)
#
# SRTProj <- AddUMAP(SRTProj, n_comp = 3)
# EachRGBSpaHeatMap(SRTProj, plot_type = "UMAP" ,combine = T, no_axis = T, pt_size = 5,title_name = "UMAP RGB plot: ")
# EachRGBSpaHeatMap(SRTProj, batch=1:2,plot_type = "UMAP", pt_size=4 ,title_name="UMAP RGB plot: ")
#
# SRTProj <- AddTSNE(SRTProj, n_comp = 3)
# EachRGBSpaHeatMap(SRTProj, batch=1:2,plot_type = "tSNE", title_name="tSNE RGB plot: ")
#
# ##
# EmbedPlot(SRTProj, item = "Cluster")
# EmbedPlot(SRTProj, item = "Cluster", legend.position='right')
# EmbedPlot(SRTProj, item = "Batch", axis_names=c("t1", "t2"))
#
# # Visualization -----------------------------------------------------------
# EachExprSpaHeatMap(SRTProj, batch=1:2, features=c("gene2", "gene3"), title_name=T)
#
# p12 <- EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"), common.legend = T)
# EachGCHeatMap(SRTProj, batch=NULL, features=c("gene2", "gene4"), common.legend = F,
#           legend.position = 'right', y_text_size=10)
#
# CCHeatMap(SRTProj, fill_legend_title="Pearson' correlation")
#
#
#
# speInt <- getIntegratedData(SRTProj)
#
#
# # DimPlot(subset(speInt,colData(speInt)$batch=='x1'), reduction='Coords')
# scater::plotReducedDim(speInt, dimred = 'aligned.iSC.MEB', colour_by = "clusters")
# scater::plotReducedDim(speInt, dimred = 'Coord', colour_by = "clusters")
#
# # Downstream analyses -----------------------------------------------------
# spe <- getGeneSpotData(SRTProj)
# spe <- runTSNE(spe)
# library(scater)
# plotPCA(spe, colour_by = "sample_id")
# plotPCA(spe, colour_by = "clusters")
# plotReducedDim(spe, dimred = 'aligned.iSC.MEB', colour_by = "clusters")
# plotReducedDim(spe, dimred = 'TSNE', colour_by = "clusters")
# deg1 <- FindDEGs(spe, use_cluster = 'clusters', cluster1='1', cluster2='2')
#
# seu <- Seurat::as.Seurat(spe)
# Idents(seu) <- seu$clusters
# DimPlot(seu)
#
# speInt <- getIntegratedData(SRTProj = SRTProj, only.var.features = F)
# head(colData(speInt))
# deg2 <- FindDEGs(speInt, use_cluster = 'clusters', cluster1='1', cluster2='2')
#
#
# ## 基于spe来写下游分析的函数
#
# # DEG analysis ------------------------------------------------------------
# deg_all <- FindAllDEGs(speInt, use_cluster = 'clusters', only.pos = T)
#
#
# # Trajectory inference ----------------------------------------------------
# speInt <- AddTrajectory(speInt, name = "")
# # speInt <- runTSNE(speInt)
# # tSNE <- calculateTSNE(x=t(reducedDim(speInt, "aligned.iSC.MEB")))
# # reducedDim(speInt, 'tSNE') <- tSNE
# plotReducedDim(speInt, dimred = 'PCA', colour_by = "clusters", text_size = 50)
# EmbedPlot(speInt, reduction = 'tSNE', colour_by = "clusters", base_size = 16)
#
# base_size <- 14; border_col <- "black";legend.position <- 'right'
# plotReducedDim(speInt, dimred = 'PCA', colour_by = "Traject.commonPseudotime") +
#   theme(axis.text.x=element_text(size=base_size, color=1),
#         axis.text.y=element_text(size=base_size, color=1),
#         axis.title.x = element_text(size=base_size+2, color='black'),
#         axis.title.y = element_text(size=base_size+2, color='black'),
#         strip.text =  element_text(size=base_size, color='black'),
#         strip.background = element_rect(
#           linetype = 'solid', color='gray3'
#         ),
#         legend.position = legend.position,
#         legend.text=element_text(size=base_size+1),
#         legend.title=element_text(size=base_size+2),
#         panel.background= element_rect(fill = 'white', color=border_col))
#
# EmbedPlot(speInt, reduction = 'PCA', colour_by = "PT")
# EachEmbedPlot(speInt, reduction = 'PCA', base_size = 20,
#               batch=c("x1", 'x2'),colour_by = "Traject.commonPseudotime")
# EachEmbedPlot(speInt, reduction = 'PCA', base_size = 20,
#               batch=c("x1", 'x2'),colour_by = "clusters")
#
# EachEmbedPlot(speInt, reduction = 'PCA', base_size = 20,
#               batch=c("x1", "x2"),colour_by = "PT")
# EmbedPlot(speInt, reduction = "tSNE", colour_by = 'clusters', base_size = 15)
#
# EachGCHeatMap(spe, features=features)
#
# # Downstream analyses functions -------------------------------------------
#
# ## If users are more familar to Seurat object.
# SRTProj2Seurat <- function(){
#
# }
#
# plot(as.matrix(posList[[3]]), col=as.numeric(clusterList[[3]]))
#
# # save(SRTProj, speInt, spe, file=paste0(SRTProj@projectMetadata$outputPath, "\\multi_SRTProj.rds"))
#
#
# # 说明 ----------------------------------------------------------------------
#
# ## 对每一个batch画图的方法：
# EachClusterSpaHeatMap # 空间聚类图
# EachExprSpaHeatMap #空间基因表达的分布图
# EachRGBSpaHeatMap # 空间的RGB图
# EachGCHeatMap   # gene*cell的热力图
# EachEmbedPlot  #embed plot
# ## 对combine的数据画的图有：
# EmbedPlot  #
