\name{sce2seurat}
\alias{sce2seurat}
\title{SingleCellExperiment object to Seurat object}
\description{
  Transfer SingleCellExperiment object to a Seurat object for preparation for \link{DR.SC} model fitting; see our [DR.SC package website](https://feiyoung.github.io/DR.SC/index.html) for more usage of \link{DR.SC}.
}
\usage{
  sce2seurat(sce, verbose= TRUE)
}
\arguments{
  \item{sce}{a SingleCellExperiment object, at least including the raw gene count expression matrix.}
  \item{verbose}{an optional logical value, whether output the information.}
}

\value{
  Return a  Seurat object.
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

  library(SingleCellExperiment)
  dir <- system.file(
    file.path("extdata", "10xVisium", "section1"),
    package = "SpatialExperiment")

  # read in counts
  fnm <- file.path(dir, "raw_feature_bc_matrix")
  sce <- DropletUtils::read10xCounts(fnm)
  colnames(sce) <- paste0("cell", 1:ncol(spe))


  seu <- sce2seurat(sce)
  head(seu)

  ## Fit DR-SC model
  library(DR.SC)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- DR.SC(seu = seu, K=4, platform = "scRNAseq")
}

