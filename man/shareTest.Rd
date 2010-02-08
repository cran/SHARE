\name{shareTest}
\alias{shareTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation Test for the Results from SHARE Algorithm}
\description{
  Permutation tests to compute the experimentwise p-values that account
  for model searching.
}
\usage{
shareTest(outObj, haploObj, status, tol = 1e-08, verbose = FALSE,
          nperm = 1000, seed=38329832)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outObj}{ the \code{share} object outputed from
    \code{cshare} function}
  \item{haploObj}{The \code{haplo} object \code{cshare} applied to}
  %%\item{phenoData}{The case-control status of the subjects}
  \item{status}{A character string indicating the column name of the phenotype in
    \code{haploObj@pheno} to be used as the clinical status in the analysis.}
  \item{tol}{The convergence parameter for the haplotype logistic regression.}
  \item{verbose}{TRUE/FALSE to decide whether to create log file for debug}
  \item{nperm}{maximal number of permutation tests }
  \item{seed}{Seed for randomization in permutation tests }
}
\details{
If the best model size is zero, there appears to be no genetic association in the region of interest.  
There is no need to perform a permutation test.  
For final models with at least 1 SNPs, we permute case-control labels 1000 times regardless of the genotypic data, carry out model searching for each permuted dataset, and compute the nominal p-value using a Wald test.  
Finally the experiment-wise p-value is computed by comparing the observed p-value to its null distribution.
}
\value{
  The experiment-wise p-value from the permutation test will be returned.
}
\references{
J. Y. Dai, M. LeBlanc, N. L. Smith, B. M. Psaty, and C. Kooperberg. SHARE: an
adaptive algorithm to select the most informative set of SNPs for genetic association.
Biostatistics, 2009. In press.

J. Besag and P. Clifford. Sequential monte carlo p-values. Biometrika, 78(2):301,
June 1, 1991.
}
\author{James Y. Dai}
\seealso{\code{\link{cshare}}}
\examples{
\dontrun{
## See vignette for more details
permuPValue <- shareTest(outObj=kerem[["Cross-Val"]],
  			 haploObj=keremHaplo,
			 status = "CF",
  			 nperm=1000 
  			 )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
