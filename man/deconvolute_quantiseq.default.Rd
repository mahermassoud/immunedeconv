% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantiseq.R
\name{deconvolute_quantiseq.default}
\alias{deconvolute_quantiseq.default}
\title{Use quanTIseq to deconvolute a gene expression matrix.}
\usage{
deconvolute_quantiseq.default(
  mix.mat,
  arrays = FALSE,
  signame = "TIL10",
  tumor = FALSE,
  mRNAscale = TRUE,
  method = "lsei",
  btotalcells = FALSE,
  rmgenes = "unassigned"
)
}
\arguments{
\item{mix.mat}{table with the gene TPM (or microarray expression values) for all samples to be deconvoluted
(Gene symbols on the first column and sample IDs on the first row). Expression data must be on non-log scale}

\item{arrays}{specifies whether expression data are from microarrays (instead of RNA-seq).
If TRUE, the "--rmgenes" parameter is set to "none".}

\item{signame}{name of the signature matrix. Currently only \code{TIL10} is available.}

\item{tumor}{specifies whether expression data are from tumor samples. If TRUE, signature genes
with high expression in tumor samples are removed.
Default: FALSE.}

\item{mRNAscale}{specifies whether cell fractions must be scaled to account for cell-type-specific
mRNA content.
Default: TRUE.}

\item{method}{deconvolution method to be used: "hampel", "huber", or "bisquare" for robust regression
with Huber, Hampel, or Tukey bisquare estimators, respectively, or "lsei" for constrained
least squares regression. The fraction of uncharacterized cells ("other") is computed only
by the "lsei" method.
Default: "lsei".}

\item{btotalcells}{compute cell densities instead of fractions
Default: FALSE}

\item{rmgenes}{Default: "default" for RNAseq, "none" for microArray data}
}
\description{
Source code from https://github.com/FFinotello/quanTIseq
}
\details{
F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder,
H. Hackl, A. Krogsdam, W. Posch, D. Wilflingseder, S. Sopper, M. Jsselsteijn,
D. Johnsons, Y. Xu, Y. Wang, M. E. Sanders, M. V. Estrada, P. Ericsson-Gonzalez,
J. Balko, N. F. de Miranda, Z. Trajanoski. "quanTIseq: quantifying immune contexture of human tumors".
bioRxiv 223180. https://doi.org/10.1101/223180.
}
