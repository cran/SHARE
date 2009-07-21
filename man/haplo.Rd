\name{haplo}
\alias{haplo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Haplotype Estimation}
\description{
  Estimate haplotype sequences and frequencies from genotype sequences by EM algorithm
}
\usage{
haplo(.Object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{A object of calss "genoSet".}
  %%\item{wgt}{A vector of weights.}
  \item{\dots}{Other arguments used in haplo.em from haplo.stats package.}
}
\details{
  By using the haplo.em function in haplo.stats package, \code{haplo}
  convert unphased genotype sequences into phased haplotype sequences by
  EM algorithm.  Haplotype frequencies are also estimated by haplp.em function.
}
\value{
  \code{haplo} returns an object of class "haplo" which
  results in the following slots:
  \item{haploSeq }{unique haplotypes.  This is from the "haplotype"
    component in the result of haplo.em funciton.}
  \item{haploFreq }{ MLE's of haplotype probabilities.  This is from
    the "hap.prob" component in the result of haplo.em funciton.}
  \item{hap1 }{the code of first haplotype of each subject.  This is from
    the "hap1code" component in the result of haplo.em funciton.}
  \item{hap2 }{the code of first haplotype of each subject.  This is from
    the "hap2code" component in the result of haplo.em funciton.}
  \item{poolHapPair }{Subject IDs. This is from the "subj.id"
    component in the result of haplo.em funciton.}
  \item{nPosHapPair }{vector for the count of haplotype pairs that
    map to each subject's marker genotypes.  This is from the "nreps"
    component in the result of haplo.em funciton.}
  \item{post }{vector of posterior probabilities of pairs of haplotypes
    for a person, given their marker phenotypes.  This is from the "post"
    component in the result of haplo.em funciton.}
}
\references{haplo.em}
\author{James Y. Dai}
%%\note{ ~~further notes~~ 
%%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%%}
%%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\seealso{\code{\link{genoSet}}, \code{\link{haplo-class}}}
\examples{
## See vignette for more details
\dontrun{
unphasedHaplo <- haplo(unphasedGeno)
unphasedHaplo
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
