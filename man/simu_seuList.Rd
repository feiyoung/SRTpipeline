\name{simu_seuList}
\alias{simu_seuList}
\title{SpatialExperiment object to Seurat object}
\description{
  Transfer SpatialExperiment object to a Seurat object for preparation for \link{DR.SC} model fitting.
}
\usage{
  simu_seuList(seu, annotated_label, spatial_coords= c("row", "col"),
                             NumBatches = 3, NumSpotsRatio=0.9, ngenes = 2000, sim_seed=1,
                             batch_facLoc=0.1, batch_facScale = 0.1)
}
\arguments{
  \item{seu}{seu is an object named Seurat, which can easily created by R package Seurat.}
  \item{annotated_label}{a string, the colname of annotations in meta.data of seu about spatial domains/cell types.}
  \item{spatial_coords}{a string vector, the spatial coordinates information in the metadata of seu, default as named "row" and "col".}
  \item{NumBatches}{an optional integer, the number of batches to simulate; default as 3.}
  \item{NumSpotsRatio}{a real between 0.1 to 1. It controls the number of spots for each data batch.}
  \item{ngenes}{}
  \item{sim_seed}{}
  \item{batch_facLoc}{an optional logical value, whether output the information.}
  \item{batch_facScale}{}
}

\details{
  The number of spots for each data batch is approximately 80\%~100\% of the number of spots in seu, the baseline seurat object.
}

\value{
  Return a  Seurat object, where the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col".
}

\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  None
}

\examples{
\dontrun{
  library(DR.SC)
  data("dlpfc151510")
  dlpfc151510
  head(dlpfc151510@meta.data)
  ## Remove the unannotated spots
  seu <- dlpfc151510[,!is.na(dlpfc151510$annotation)]
  seulist_sim <- simu_seuList(seu, annotated_label = 'annotation', NumBatches = 2,
                          NumSpotsRatio = 0.5,spatial_coords=c("row", "col"), ngenes=50)

  seulist_sim <- lapply(seulist_sim, AddCoord2Reduc)
  seulist_sim2 <- lapply(seulist_sim, function(x){
    Idents(x) <- factor(x$Group)
    return(x)
  } )

  library(patchwork)
  pList <- lapply(seulist_sim2, DimPlot, reduction = "position", pt.size = 2)
  wrap_plots(pList, nrow=1)
  }
}

