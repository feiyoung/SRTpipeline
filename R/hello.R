# #
# # # library(DR.SC)
# # # dir.file <- "F:\\Research paper\\IntegrateDRcluster\\iDR-SC\\Human_Breast_Cancer\\BlockA_Section"
# # # seuList <- list()
# # # for (r in 1:2) {
# # #   message("r = ", r)
# # #   seuList[[r]] <- read10XVisium(paste0(dir.file, r))
# # # }
# # # bc2 <- seuList
# # # ## add data
# # # usethis::use_data(bc2)
# #
# # # ###
# # # setwd("D:\\LearnFiles\\Research paper\\ProPCA\\DR-SC.Analysis\\data\\DLPFC_data\\")
# # # dat_tmp <- readRDS("151672.rds")
# # # suppressPackageStartupMessages(library(SingleCellExperiment))
# # # library(Seurat)
# # # meta_data <- as.data.frame(colData(dat_tmp))
# # # dlpfc_151672 <- CreateSeuratObject(counts = dat_tmp@assays@data$counts, meta.data = meta_data)
# # # setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\DataPRECAST")
# # # usethis::use_data(dlpfc_151672)
# #
# #
# #
# # # Simulated data ----------------------------------------------------------
# #
# #
# # ## load data
# # setwd("D:/LearnFiles/Research paper/IntegrateDRcluster/simu/simuR/RealDataBasedSimu")
# # source("help_func.R")
# # load("brain1507_1669_1673_yPosV.rds")
# #
# #
# # yList <- lapply(yList, as.numeric)
# # for(r in 1:3){
# #   y <- yList[[r]]
# #   y[is.na(y)] <- 7
# #   yList[[r]] <- y
# # }
# #
# # sapply(yList, table)
# #
# # ## generate simulated data:W, W_r, Lambda_r,
# # p <- 2000
# # q = 10
# # sigmaW=c(0.5,0.8,1);
# # sigmaZ = c(1,2, 0.5);
# # qvec=rep(2, 3); # dimension of lantent features of batch effects
# # require(MASS)
# # n_vec <- sapply(yList, length)
# # ## generate deterministic parameters, fixed after generation
# # set.seed(1)
# # sigma2 <- 2
# # LamMat <- rbind(sigma2*(1 + 0.5 * abs(rnorm(p, sd=3))),
# #                 sigma2*(1 + 0.2*(runif(p))),
# #                 sigma2*(1 + 1*(runif(p))))
# #
# #
# # W <- matrix(rnorm(p*q), p, q)
# # W <- qr.Q(qr(W))
# # Wlist <- list()
# # for(r in 1:3){
# #   set.seed(r+1) # sigma12 control the correlation strength, if sigma12 is small, correlation is strong
# #   Wtt1 <- matrix(rnorm(p* qvec[r]), p, qvec[r])
# #   W1 <- Wtt1 + sigmaW[r] * matrix(rnorm(p* qvec[r]), p, qvec[r])
# #   W1 <- qr.Q(qr(W1))
# #   Wlist[[r]] <- W1
# #   cat('cor(W,W1)=', mean(cancor(W, Wlist[[r]])$cor), '\n')
# # }
# #
# # K <- 7
# # mu <- matrix(c(rep(0,q),
# #                rep(0,q),
# #                rep(0,q),
# #                rep(0,q),
# #                rep(0,q),
# #                rep(0,q),
# #                rep(0,q)), ncol = K)
# # mu_value <- 8
# # mu[1,1] = mu_value
# # mu[2,2] = mu_value
# # mu[3,3] = mu_value
# # mu[4,4] = mu_value
# # mu[5,5] = -mu_value
# # mu[6,6] = -mu_value
# # mu[7:q,7] = -mu_value
# #
# # ## heter covariance components
# # diagmat = array(0, dim = c(q, q, K))
# # for(k in 1:K){
# #   diag(diagmat[,,k]) <- 1
# # }
# # diag(diagmat[,,1]) = c(10,rep(1,q-1))
# # diag(diagmat[,,2]) = c(1,10,rep(1,q-2))
# # diag(diagmat[,,3]) = c(rep(1,2),10, rep(1,q-3))
# # diag(diagmat[,,4]) = c(rep(1,3),10, rep(1,q-4))
# #
# #
# #
# # Mu <- t(mu)
# # Sigma <- diagmat
# # tauMat <- matrix(0, 3, q)
# # tauMat[2,1] <- 10; tauMat[2,2] <- -10
# # tauMat[3, ] <- rep(5, q);
# #
# # tau0Mat <- matrix(NA, 3, p)
# # for(r in 1:3){
# #   set.seed(r+5)
# #   tau0 <- rnorm(p, sd=2)
# #   tau0Mat[r, ] <- tau0
# # }
# #
# #
# # i <- 1
# # cat('i = ', i, '\n')
# # set.seed(i)
# #
# # ## generate low-dimensional embedding with biological effects
# # Zlist <- list()
# # for(r in 1:3){
# #   Z_tmp <- matrix(0, n_vec[r], q)
# #   for(k in 1:K){
# #     nk <- sum(yList[[r]]==k) # include conditional and sequencing batch
# #     if(nk > 0)
# #       Z_tmp[yList[[r]]==k, ] <- MASS::mvrnorm(nk, Mu[k,]+tauMat[r,],
# #                                               Sigma[,,k]+3*(r-1)*diag(q))
# #   }
# #   Zlist[[r]] <- Z_tmp
# # }
# # sapply(Zlist, dim)
# #
# #
# # ## batch effect
# # Zrlist <- list()
# # for(r in 1:3){
# #   set.seed(r+10)
# #   Zrlist[[r]] <- matrix(rnorm(n_vec[r]* qvec[r], sd=sigmaZ[r]), n_vec[r], qvec[r])
# # }
# # sapply(Zrlist, dim)
# # XtList <- list()
# # for(r in 1:3){
# #   message("r = ", r)
# #   X1 <- (Zlist[[r]] + VList[[r]] ) %*% t(W) + Zrlist[[r]] %*% t(Wlist[[r]])+
# #     MASS::mvrnorm(n_vec[r], rep(0,p), diag(LamMat[r,]))
# #
# #   tauMat0 <- matrix(tau0Mat[r, ], n_vec[r], p, byrow = T)
# #   Eta <- exp((X1 + tauMat0 ))
# #   summary(colSums(Eta))
# #   #X1 <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
# #   XtList[[r]] <- matrix(rpois(n_vec[r]*p, Eta), n_vec[r], p)
# # }
# # for(r in 1:3){
# #   row.names(XtList[[r]]) <- paste0("S1_spot",1:nrow(XtList[[r]]))
# #   colnames(XtList[[r]]) <- paste0("gene_",1:ncol(XtList[[r]]))
# # }
# #
# # XList <- lapply(XtList, function(x) log(1+x))
# #
# # ##### Start to integration
# # hK <- 7
# # hq <- 15
# # #### Method iDR-SC
# # XList1 <- lapply(XList, scale, scale=FALSE)
# # tic <- proc.time() #
# # set.seed(1)
# # resList <- iDR.SCv3::idrsc(XList1,posList=posList, q=hq, K=hK,
# #                            platform = 'Visium',maxIter = 30,
# #                            Sigma_equal =F, coreNum=length(hK))
# # toc <- proc.time()
# # timeMat[i,1] <- toc[3] - tic[3]
# # library(iDR.SCv3)
# # reslist <- selectModel(resList, pen_const=1)
# # str(reslist)
# #
# # ## cluster performance #
# # #### Clustering metric
# # cluster_metric <- function(hy, y, type='ARI'){
# #
# #   require(mclust)
# #   require(aricode)
# #   switch(type,
# #          ARI= adjustedRandIndex(hy, y),
# #          NMI = NMI(as.vector(hy), y))
# # }
# # cluster_metric(unlist(reslist$cluster), unlist(yList))
# # cluster_metric(unlist(reslist$cluster), unlist(yList), type='NMI')
# #
# # plot(posList[[1]], col=reslist$cluster[[1]])
# #
# # library(Seurat)
# # data_simu <- list()
# # for(r in 1:3){
# #   ## r <- 1
# #   meta_data <- data.frame(row=posList[[r]][,1], col=posList[[r]][,2])
# #   meta_data$sample <- paste0("S",r)
# #   row.names(meta_data) <- row.names(XtList[[r]])
# #   seu_tmp <- CreateSeuratObject(counts= t(XtList[[r]]), meta.data = meta_data)
# #   seu_tmp[["RNA"]]@data <- t(XList[[r]])
# #   data_simu[[r]] <- seu_tmp
# # }
# # data_simu
# # setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\DataPRECAST")
# # usethis::use_data(data_simu)
# #
# # library(DataPRECAST)
# # data("data_simu")
# # library(Seurat)
# # row.names(data_simu[[1]])[1:10]
# # data_simu[[1]][["RNA"]]@counts[1:4,1:5]
# #
# # data_simu1 <- data_simu
# # library(Seurat)
# # data_simu <- list()
# # for(r in 1:3){
# #   ## r <- 1
# #   meta_data <- data_simu1[[r]]@meta.data
# #   meta_data$true_cluster <- yList[[r]]
# #  count <- data_simu1[[r]][["RNA"]]@counts
# #  data <- data_simu1[[r]][["RNA"]]@data
# #  row.names(count) <- row.names(data) <- paste0("gene", 1:nrow(count))
# #   seu_tmp <- CreateSeuratObject(counts= count, meta.data = meta_data)
# #   seu_tmp[["RNA"]]@data <- data
# #   data_simu[[r]] <- seu_tmp
# # }
# # usethis::use_data(data_simu, overwrite = T)
# #
# #
# # true_cluster <- lapply(data_simu, function(x) x$true_cluster)
# # str(true_cluster)
# # mclust::adjustedRandIndex(unlist(PRECASTObj@resList$cluster), unlist(true_cluster))
# #
# ##
# # library(Seurat)
# # githubURL <- "https://github.com/feiyoung/DataPRECAST/blob/main/data/data_simu.rda?raw=true"
# # load(url(githubURL), verbose = T)
# # library(RCurl)
# # ## paste URL to make it easier to read code (cosmetic!)
# # dat_url <- githubURL
# # f <- getBinaryURL(dat_url)
# # L <- load(rawConnection(f))
# #
# # githubURL <- ("https://github.com/feiyoung/DataPRECAST/blob/main/data/data_simu.rda?raw=true")
# # download.file(githubURL,"my.rds",mode='wb')
# # load(file="./my.rds" )
# # load("data_simu.rda")
# # file_url <- "https://github.com/TarekDib03/ExploratoryDataAnalysisCoursera/blob/master/maacs.Rda?raw=true"
# # load(url(file_url))
# #
# #
# # load("./data/data_simu.rda")
#
#
#
# # library(SpatialExperiment)
# # dir <- system.file(
# #   file.path("extdata", "10xVisium", "section1"),
# #   package = "SpatialExperiment")
# #
# # # read in counts
# # fnm <- file.path(dir, "raw_feature_bc_matrix")
# # sce <- DropletUtils::read10xCounts(fnm)
# #
# # # read in image data
# # img <- readImgData(
# #   path = file.path(dir, "spatial"),
# #   sample_id="foo")
# #
# # # read in spatial coordinates
# # fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
# # xyz <- read.csv(fnm, header = FALSE,
# #                 col.names = c(
# #                   "barcode", "in_tissue", "array_row", "array_col",
# #                   "pxl_row_in_fullres", "pxl_col_in_fullres"))
# #
# # # construct observation & feature metadata
# # rd <- S4Vectors::DataFrame(
# #   symbol = rowData(sce)$Symbol)
# #
# # # construct 'SpatialExperiment'
# # (spe <- SpatialExperiment(
# #   assays = list(counts = assay(sce)),
# #   colData = colData(sce), rowData = rd, imgData = img,
# #   spatialData=DataFrame(xyz),
# #   spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"),
# #   sample_id="foo"))
# #
# # colnames(spe) <- paste0("spot", 1:ncol(spe))
# #
# #
# # seu <- spe_to_seurat(spe)
# # head(seu)
# # seu2 <- seu
# # seu2$sample_id <- paste0(seu2$sample_id, "2")
# # library(SeuratObject)
# # seu2 <- RenameCells(seu2,  paste0(colnames(seu2), "2") )
# # seu_all <- merge(seu, seu2)
# # head(seu_all@meta.data)
# # table(seu_all$sample_id)
# #
# # spe_to_seuratList(spe, batch='sample_id')
# #
# # spe_to_seuratList(list(spe))
# #
# # sce <- SingleCellExperiment(spe)
# # sce <- SingleCellExperiment(seu_all)
# # sce
# #
# # seuList <- seu2seuList(seu_all, batch='sample_id')
# # seu <- NormalizeData(seu)
# # seu <- FindVariableFeatures(seu)
# # DR.SC(seu = seu, K=4, platform = "Visium", approxPCA=T)
# #
# # library(DR.SC)
# # data("dlpfc151510")
# # dlpfc151510
# # head(dlpfc151510@meta.data)
# # ## Remove the unannotated spots
# # seu <- dlpfc151510[,!is.na(dlpfc151510$annotation)]
# # seulist_sim <- simu_seuList(seu, annotated_label = 'annotation', NumBatches = 2,
# #                             NumSpotsRatio = 0.2,
# #                             spatial_coords=c("row", "col"), ngenes=50)
# #
# #
# # seulist_sim <- lapply(seulist_sim, AddCoord2Reduc)
# # seulist_sim2 <- lapply(seulist_sim, function(x){
# #   Idents(x) <- factor(x$Group)
# #   return(x)
# # } )
# #
# # library(patchwork)
# # pList <- lapply(seulist_sim2, DimPlot, reduction = "position", pt.size = 2)
# # wrap_plots(pList, nrow=1)
# # FeaturePlot(seulist_sim2[[1]], features= row.names(seulist_sim2[[1]])[10], reduction = "position", pt.size = 2)
#
#
# # ##spe2seuratList
# # colData(spe)$batch_id <- rep(c("a", "b"), each=25)
# # seuList <- spe2seuratList(spe, batch = 'batch_id')
# # library(PRECAST)
# # Pobject <- CreatePRECASTObject(seuList=seuList, selectGenesMethod = "HVGs",
# #                                premin.spots = 0, premin.features = 0,
# #                                postmin.features = 0, postmin.spots = 0, verbose = F)
# #
# # Pobject <- AddAdjList(Pobject, platform = 'Visium')
# # Pobject <- AddParSetting(Pobject)
# # Pobject <- PRECAST(Pobject, K= 4)
# #
# #
# # ## Change getAdj_auto
# # getneighborhood_fast <- DR.SC:::getneighborhood_fast
# # pos <- cbind(1:50, 2:51)
# # getAdj_auto <- function(pos, lower.med=4, upper.med=6, radius.upper= NULL){
# #   if (!inherits(pos, "matrix"))
# #     stop("method is only for  matrix object!")
# #
# #   if(is.null(radius.upper)){
# #     n <- nrow(pos)
# #     idx <- sample(n, min(n, 100))
# #     radius.upper <- max(dist(pos[idx,]))
# #   }
# #
# #   radius.lower <- 1
# #   Adj_sp <- getneighborhood_fast(pos, radius=radius.upper)
# #   Med <- summary(Matrix::rowSums(Adj_sp))['Median']
# #   if(Med < lower.med) stop("The radius.upper is too smaller that cannot find median neighbors greater than 4.")
# #   start.radius <- 1
# #   Med <- 0
# #   message("Find the adjacency matrix by bisection method...")
# #   maxIter <- 30
# #   k <- 1
# #   while(!(Med >= lower.med && Med <=upper.med)){ # ensure that each spot has about 4~6 neighborhoods in median.
# #
# #     Adj_sp <- getneighborhood_fast(pos, radius=start.radius)
# #     Med <- summary(Matrix::rowSums(Adj_sp))['Median']
# #     if(k > 1){
# #       message("Current radius is ", round(start.radius, 2)," and its Median of neighborhoods is ", Med)
# #     }
# #     if(Med < lower.med){
# #       radius.lower <- start.radius
# #       start.radius <- (radius.lower + radius.upper)/2
# #     }else if(Med >upper.med){
# #       radius.upper <- start.radius
# #       start.radius <- (radius.lower + radius.upper)/2
# #     }
# #
# #     if(k > maxIter) {
# #       message("Reach the maximum iteration but can not find a proper radius!")
# #       break;
# #     }
# #     k <- k + 1
# #   }
# #
# #   return(Adj_sp)
# # }
# # getAdj_auto(pos)
