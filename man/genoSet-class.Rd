\name{genoSet}
\docType{class}
\alias{genoSet}
\alias{genoSet-class}
\alias{genoSeq,genoSet-method}
\alias{genoSeq}
\alias{phenoData,genoSet-method}
\alias{phenoData}
\alias{haplo,genoSet-method}

\title{Genotype Sequence Set}
\description{A class for storing genotype sequences and phenotype information.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("genoSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{genoSeq}:}{Object of class \code{"data.frame"} to store
      allelic data with sequence names as row names and SNP names
      as column names.}
    \item{\code{phenoData}:}{Object of class \code{"data.frame"} to
      store phenotype information for each genotype sequence.  The names
    of the phenotype are used as the column names, and the row names
    should match the sequence names in the \code{genoSeq} slot. }
  }
}
\section{Methods}{
  \describe{
    \item{genoSeq}{\code{signature(.Object = "genoSet")}: extract the
      \code{"data.frame"} object in the \code{genoSeq} slot.}
    \item{haplo}{\code{signature(.Object = "genoSet")}: to estimate the
      haplotype sequences by EM algorithm in the \code{haplo.stats} package. }
%%    \item{nameSeq}{\code{signature(.Object = "genoSet")}: output a
%%      vector of the names of the genotype sequences. }
%%    \item{nameSNP}{\code{signature(.Object = "genoSet")}: output a
%%      vector of the SNP names. }
%%    \item{nSeq}{\code{signature(.Object = "genoSet")}: how many genotype
%%    sequences. }
%%    \item{nSNP}{\code{signature(.Object = "genoSet")}: how many SNPs }
    \item{phenoData}{\code{signature(.Object = "genoSet")}: return the
      \code{"data.frame"} object in the \code{phenoData} slot.}
  }
}
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
\examples{
showClass("genoSet")
## See vignette for more details
\dontrun{
unphasedGeno <- new("genoSet",
                    genoSeq = data.frame(keremRandAllele),
                    phenoData = data.frame(CF=keremRandStatus)
                    )
unphasedGeno
}
}
\keyword{classes}
