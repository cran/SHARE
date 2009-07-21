\name{SHARE-package}
\alias{SHARE-package}
\alias{SHARE}
\docType{package}
\title{
  An adaptive algorithm to select the most informative set of SNPs for genetic association  
}
\description{
  This is the R package to perform the adaptive algorithm, developed by
  James Dai, et al., to select the most informative SNP set for genetic association.
}
\details{
%%\tabular{ll}{
%%Package: \tab SHARE\cr
%%Type: \tab Package\cr
%%Version: \tab 0.0.20\cr
%%Date: \tab 2009-07-08\cr
%%License: \tab GPL (>= 2)\cr
%%LazyLoad: \tab no\cr
%%}
Association studies have been widely used to identify genetic liability
variants for complex diseases.  While scanning the chromosomal region
one SNP at a time may not fully explore linkage disequilibrium (LD),
haplotype analyses tend to require a fairly large number of parameters,
thus potentially losing power.  Clustering algorithms, such as the
cladistic approach, have been proposed to reduce the dimensionality, yet
they have important limitations.  We propose the SHARE algorithm
that seeks the most informative set of SNPs for genetic
association in a targeted region by growing/shrinking haplotypes with
one more/less SNP in a stepwise fashion, and comparing prediction errors
of different models via cross-validation.
}
\author{
  James Y. Dai & Ting-Yuan Liu
}
\references{
Dai, J. Y., LeBlanc, M., Smith, N. L., Psaty, B. M. and Kooperberg, C. (2009). SHARE: 
an adaptive algorithm to select the most informative set of SNPs for genetic associat ion. 
Biostatistics, In Press.
}
\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
%%\examples{
%%~~ simple examples of the most important functions ~~
%%}
