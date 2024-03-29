% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_newinterface.R
\name{model_set_PRECAST}
\alias{model_set_PRECAST}
\title{Set parameters for PRECAST model}
\usage{
model_set_PRECAST(
  Sigma_equal = FALSE,
  Sigma_diag = TRUE,
  mix_prop_heter = TRUE,
  error_heter = TRUE,
  Sp2 = TRUE,
  wpca_int = FALSE,
  int.model = "EEE",
  coreNum = 5,
  coreNum_int = coreNum,
  beta_grid = seq(0.2, 4, by = 0.2),
  maxIter_ICM = 6,
  maxIter = 20,
  epsLogLik = 1e-05,
  verbose = TRUE,
  seed = 1
)
}
\arguments{
\item{Sigma_equal}{a logical value, whether set the mixture covariance matrices equal, default as FALSE.}

\item{Sigma_diag}{a logical value, whether set the mixture covariance matrices diagonal, default as TRUE.}

\item{mix_prop_heter}{a logical value, whether set the smoothing parameter of each data batch to be different, default as TRUE.}

\item{error_heter}{a logical value, whether use the heterogenous error for PRECAST model, default as TRUE.}

\item{Sp2}{a logical value, whether add the intrisical CAR component in the PRECAST model, default as TRUE.}

\item{wpca_int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as FALSE which means the ordinary PCA is used.}

\item{int.model}{an optional string, specify which Gaussian mixture model is used in evaluting the initial values for DR-SC, default as "EEE"; and see Mclust for more models' names.}

\item{coreNum}{an integer value to decide how many cores are used to run  in parallel, default as 5.}

\item{coreNum_int}{an integer value to decide how many cores are used in parallel computing initial values.}

\item{beta_grid}{a numeric vector, specify the smoothness parameter of Random Markov Field. The default is seq(0.2,4,0.2).}

\item{maxIter_ICM}{the maximum iteration of ICM step, default as 6.}

\item{maxIter}{the maximum iteration of ICM-EM algorithm. The default is 25}

\item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}

\item{verbose}{a logical value, whether output the information in iteration.}

\item{seed}{a postive integer, the random seed to be set.}
}
\description{
Prepare parameters setup for PRECAST model fitting.
}
\examples{
model_set_PRECAST(Sigma_equal=TRUE, coreNum=5)

}
