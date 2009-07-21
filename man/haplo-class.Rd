\name{haplo-class}
\docType{class}
\alias{haplo-class}

\title{Haplotype Estimation Set}
\description{A class for storing haplotype estimation information.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("haplo", ...)}.
%%	 ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{hap1}:}{the code of first haplotype of each subject.  This is from
      the "hap1code" component in the result of haplo.em funciton.}
    \item{\code{hap2}:}{the code of first haplotype of each subject.  This is from
      the "hap2code" component in the result of haplo.em funciton.}
    \item{\code{poolHapPair}}{Subject IDs. This is from the "subj.id"
      component in the result of haplo.em funciton.}
    \item{\code{nPosHapPair}:}{vector for the count of haplotype pairs that
      map to each subject's marker genotypes.  This is from the "nreps"
      component in the result of haplo.em funciton.}
    \item{\code{post}:}{vector of posterior probabilities of pairs of haplotypes
      for a person, given their marker phenotypes.  This is from the "post"
      component in the result of haplo.em funciton.}
    \item{\code{haploSeq}:}{Object of class \code{"data.frame"} to store
      the unique haplotype sequences.  This is from the "haplotype"
      component in the result of haplo.em funciton.}
    \item{\code{haploFreq}:}{Object of class \code{"vector"} to store
      the MLE's of haplotype probabilities.  This is from
      the "hap.prob" component in the result of haplo.em funciton.}
    \item{\code{pheno}:}{Object of class \code{"data.frame"} to store
      the phenotype information}
  }
}
\section{Extends}{
Class \code{"\linkS4class{haploSet}"}, directly.
}
%%\section{Methods}{
%%  \describe{
%%    \item{show}{\code{signature(object = "haplo")}: ... }
%%	 }
%%}
%%\references{ ~put references to the literature/web site here ~ }
\author{Ting-Yuan Liu}
%%\note{ ~~further notes~~ }
%%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%%
%%\seealso{
%%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%%}
\seealso{\code{\link{haplo}}, \linkS4class{haploSet}}
\examples{
showClass("haplo")
}
\keyword{classes}
