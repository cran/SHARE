\name{haploSet-class}
\docType{class}
\alias{haploSet}
\alias{haploSet-class}
\alias{haploFreq,haploSet-method}
\alias{haploFreq}
\alias{haploSeq,haploSet-method}
\alias{haploSeq}
%%\alias{nameSeq,haploSet-method}
%%\alias{nameSNP,haploSet-method}
%%\alias{nSeq,haploSet-method}
%%\alias{nSNP,haploSet-method}

\title{Haplotype Set}
\description{A class to store the haplotype sequences and their frequencies.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("haploSet",
  ...)}.  
}
\section{Slots}{
  \describe{
    \item{\code{haploSeq}:}{Object of class \code{"data.frame"} to store
    the haplotype sequences.  Names of haplotype sequences are used as
    the row names of the \code{"data.frame"}, and names of SNPs are used
    as the column names.}
    \item{\code{haploFreq}:}{Object of class \code{"vector"} to store
      the frequencies of the haplotype sequences in the \code{haploSeq}
      slot.  Names of each elements must match the names of haplotype sequences.}
  }
}
\section{Methods}{
  \describe{
    \item{haploFreq}{\code{signature(.Object = "haploSet")}: the
      \code{"vector"} object in the \code{haploFreq} slot.}
    \item{haploSeq}{\code{signature(.Object = "haploSet")}: the
      \code{"data.frame"} object in the \code{haploSeq} slot.}
%%    \item{nameSeq}{\code{signature(.Object = "haploSet")}: name of the
%%      haplotype sequences. }
%%    \item{nameSNP}{\code{signature(.Object = "haploSet")}: name of the SNPs. }
%%    \item{nSeq}{\code{signature(.Object = "haploSet")}: number of
%%      haplotype sequences }
%%    \item{nSNP}{\code{signature(.Object = "haploSet")}: number of SNPs }
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
\seealso{\code{\link{haplo-class}}}
\examples{
showClass("haploSet")
}
\keyword{classes}
