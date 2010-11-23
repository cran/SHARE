\name{share-class}
\docType{class}
\alias{share-class}
\alias{share}
\alias{dplot,share-method}
\alias{dplot}

\title{SHARE output}
\description{A class to store the result of SHARE algorithm}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("share", ...)}.
%%	 ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{uhap}:}{Object of class \code{"vector"};  the
      concatenated haplotype sequences}
    \item{\code{weight}:}{Object of class \code{"vector"};  the weight
    of each haplotype}
    \item{\code{nFold}:}{Object of class \code{"numeric"};  the fold
      size specified while calling the \code{cshare} function to
    perform the cross-validation approach}
    \item{\code{maxSNP}:}{Object of class \code{"numeric"};  the
      maximum SNP size specified while calling the \code{cshare} functionused}
    \item{\code{deviance}:}{Object of class \code{"vector"};  the
      prediction deviance from cross-validation approach, or the BIC
      from BIC approach}
    \item{\code{bestsize}:}{Object of class \code{"numeric"};  the
      size of most informative SNPs}
    \item{\code{finalHap}:}{Object of class \code{"data.frame"};  the
      result of estimated haplotype sequences with the selected SNPs}
    \item{\code{finalHapFreq}:}{Object of class \code{"vector"};  the
      frequencies of the estimated haplotype sequences with the selected SNPs }
    \item{\code{finalHapTest}:}{Object of class \code{"data.frame"}; 
      the hypothesis test results of the estimated haplotype sequences
      with the selected SNPs }
    \item{\code{globalP}:}{Object of class \code{"numeric"}; the
      global p-values}
    \item{\code{modelmethod}:}{Object of class \code{"character"};  the
      mothod of model selection specified while calling the
      \code{cshare} function}
    \item{hap1 }{the code of first haplotype of each subject.}
    \item{hap2 }{the code of first haplotype of each subject.}
    \item{poolHapPair }{Subject IDs.}
    \item{nPosHapPair }{vector for the count of haplotype pairs that
      map to each subject's marker genotypes.}
    \item{post }{vector of posterior probabilities of pairs of haplotypes
      for a person, given their marker phenotypes.}
    \item{\code{haploSeq}:}{Object of class \code{"data.frame"}; the
      original haplotype sequences with all SNPs}
    \item{\code{haploFreq}:}{Object of class \code{"vector"}; the
      frequencies of the original haplotype sequences with all SNPs}
    \item{\code{pheno}:}{Object of class \code{"data.frame"} to store
      the phenotype information}
    \item{\code{nngcov}:}{Object of class \code{"numeric"}; the number
		of non-genetic covariates}
    \item{\code{ngcov}:}{Object of class \code{"vector"}; the vector 
		of non-genetic covariates}
    \item{\code{inherit}:}{Character string to show which mode of
      inheritance for haplotype effects is specified.
    }
  }
}
\section{Extends}{
Class \code{"\linkS4class{haplo}"}, directly.
Class \code{"\linkS4class{haploSet}"}, by class "haplo", distance 2.
}
\section{Methods}{
  \describe{
    \item{dplot}{\code{signature(shareObj = "share")}: to create the
      deviance plot to show the estimation of SNP size }
  }
}
%%\references{ ~put references to the literature/web site here ~ }
\author{James Y. Dai & Ting-Yuan Liu}
%%\note{ ~~further notes~~ }
%%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%%
%%\seealso{
%%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%%}
\seealso{\code{\link{cshare}}, \linkS4class{haplo}, \linkS4class{haploSet}}
\examples{
showClass("share")

## See vignette for more details
\dontrun{
dplot(unphasedKerem[["Cross-Val"]])
dplot(unphasedKerem[["BIC"]])
}
}
\keyword{classes}
