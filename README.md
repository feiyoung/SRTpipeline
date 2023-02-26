<!-- README.md is generated from README.Rmd. Please edit that file -->

SRTpipeline
===========

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/feiyoung/SRTpipeline.svg?branch=master)](https://travis-ci.com/feiyoung/SRTpipeline)
<!-- badges: end -->

SRTpipeline is a powerful and easy-to-use R package for processing and
analyzing spaitally resolved transcriptomics (SRT) data by providing
step-by-step tutorials. SRTpipeline integrates a series of our developed
methods and commonly used analyses tools. SRTpipeline is able to handle
single data batch and multiple data batches by considering the
non-cellular effects such as batch effects.

For single data batch, it has the following functions:

-   (Spatial) dimension reduction
-   Spatial clustering
-   Joint dimension reduction and spatial clustering
-   SRT embeddings for visualization
-   DEG analysis
-   Trajectory inference
-   Conditional spaital variation analysis
-   SRT deconvolution analysis

For multiple data batches, it has the following functions:

-   (Spatial) dimension reduction by removing batch effects
-   Spatial clustering by removing batch effects
-   SRT embeddings for visualization
-   Recover comparable gene expression matrices among datasets
-   Combined DEG analysis
-   Combined trajectory inference

------------------------------------------------------------------------

Installation
------------

You can install the development version of SRTpipeline from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("feiyoung/SRTpipeline")

------------------------------------------------------------------------

<!-- # SRTpipeline <img src='man/figures/logo.png' align="right" height="139" /> -->
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/SRTpipeline)](https://cran.r-project.org/package=SRTpipeline)
[![](https://cranlogs.r-pkg.org/badges/SRTpipeline?color=orange)](https://cran.r-project.org/package=SRTpipeline)
[![Coverage
Status](https://img.shields.io/codecov/c/github/feiyoung/SRTpipeline/master.svg)](https://codecov.io/github/feiyoung/SRTpipeline?branch=master)
[![](https://badges.feiyoung.org/184_status.svg)](https://github.com/feiyoung/software-review/issues/184)
<!-- badges: end -->

<!-- [`SRTpipeline`](https://docs.feiyoung.org/SRTpipeline/) is a `R` package devoted to -->
<!-- analyzing spatially resovled transcriptomics data. -->
<!-- <!-- <a href="http://www.irea.cnr.it/en/"> <img src="man/figures/irea_logo.png" height="40" align="left" /></a>  -->

–&gt;

<!-- <span style='font-style:italic;font-weight:bold;'>`{SRTpipeline}` was developed by XX. -->
<!-- [It is dedicated to the memory of Lorenzo](https://docs.feiyoung.org/SRTpipeline/articles/lorenzo.html).</span> -->

Citation
--------

To cite `{SRTpipeline}` please use:

Yi Yang, Xingjie Shi, Wei Liu, Qiuzhong Zhou, Mai Chan Lau, Jeffrey Chun
Tatt Lim, Lei Sun, Cedric Chuan Young Ng, Joe Yeong, Jin Liu, SC-MEB:
spatial clustering with hidden Markov random field using empirical
Bayes, Briefings in Bioinformatics, Volume 23, Issue 1, January 2022,
bbab466,
<a href="https://doi.org/10.1093/bib/bbab466" class="uri">https://doi.org/10.1093/bib/bbab466</a>

Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie
Shi, Jin Liu, Joint dimension reduction and clustering analysis of
single-cell RNA-seq and spatial transcriptomics data, Nucleic Acids
Research, Volume 50, Issue 12, 8 July 2022, Page e72,
<a href="https://doi.org/10.1093/nar/gkac219" class="uri">https://doi.org/10.1093/nar/gkac219</a>

Wei Liu, Xu Liao, Ziye Luo, Yi Yang, Mai Chan Lau, Yuling Jiao, Xingjie
Shi, Weiwei Zhai, Hongkai Ji, Joe Yeong, Jin Liu. Probabilistic
embedding and clustering with alignment for spatial transcriptomics data
integration with PRECAST. Nat Commun 14, 296 (2023).
<a href="https://doi.org/10.1038/s41467-023-35947-w" class="uri">https://doi.org/10.1038/s41467-023-35947-w</a>

Xiao Zhang, Wei Liu, Fangda Song, Jin Liu, iSC.MEB: an R package for
multi-sample spatial clustering analysis of spatial transcriptomics
data, Bioinformatics Advances, 2023;, vbad019,
<a href="https://doi.org/10.1093/bioadv/vbad019" class="uri">https://doi.org/10.1093/bioadv/vbad019</a>

Xingjie Shi, Yi Yang, Xiaohui Ma, Zhenxing Guo, Jin Liu, Probabilistic
cell/domain-type assignment of spatial transcriptomics data with
SpatialAnno, bioRxiv 2023.02.08.527590; doi:
<a href="https://doi.org/10.1101/2023.02.08.527590" class="uri">https://doi.org/10.1101/2023.02.08.527590</a>

------------------------------------------------------------------------

Website
-------

For more information, documentation and examples of use, see also the
`{SRTpipeline}` website at
[SRTpipeline](https://github.com/SRTpipeline/SRTpipeline).

------------------------------------------------------------------------
