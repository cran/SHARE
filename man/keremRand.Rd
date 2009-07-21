\name{keremRand}
\alias{keremRandSeq}
\alias{keremRandStatus}
\alias{keremRandAllele}

\docType{data}
\title{Psudo-subjects from Kerem's Cystic Fibrosis Data}
\description{
  This datasets contains the psudo-subjects created from cystic fibrosis
  data in Kerem et al. (1989).  
}
\usage{data(keremRand)}
\format{
  Here is the list of the 23 alleles:
  \describe{
    \item{\code{locus_01}}{Probe: metD;		Enzyme: Ban I}
    \item{\code{locus_02}}{Probe: metD;		Enzyme: Taq I}
    \item{\code{locus_03}}{Probe: metH;		Enzyme: Taq I}
    \item{\code{locus_04}}{Probe: E6;		Enzyme: Taq I}
    \item{\code{locus_05}}{Probe: E7;		Enzyme: Taq I}
    \item{\code{locus_06}}{Probe: pH131;	Enzyme: Hinf I}
    \item{\code{locus_07}}{Probe: W3Dl.4;	Enzyme: Hind III}
    \item{\code{locus_08}}{Probe: H2.3A (XV2C);	Enzyme: Taq I}
    \item{\code{locus_09}}{Probe: EG1.4;	Enzyme: Hinc II}
    \item{\code{locus_10}}{Probe: EG1.4;	Enzyme: Bgl II}
    \item{\code{locus_11}}{Probe: JG2El (KM19);	Enzyme: Pst I}
    \item{\code{locus_12}}{Probe: E2.6 (E.9);	Enzyme: Msp I}
    \item{\code{locus_13}}{Probe: H2.8A;	Enzyme: Nco I}
    \item{\code{locus_14}}{Probe: E4.1 (Mp6d.9);Enzyme: Msp I}
    \item{\code{locus_15}}{Probe: J44;		Enzyme: Xba I}
    \item{\code{locus_16}}{Probe: 10-1X.6;	Enzyme: Acc I}
    \item{\code{locus_17}}{Probe: 10-lX.6;	Enzyme: Hae III}
    \item{\code{locus_18}}{Probe: T6/20;	Enzyme: Msp I}
    \item{\code{locus_19}}{Probe: H1.3;		Enzyme: Nco I}
    \item{\code{locus_20}}{Probe: CEL.0;	Enzyme: Nde I}
    \item{\code{locus_21}}{Probe: J32;		Enzyme: Sac I}
    \item{\code{locus_22}}{Probe: J3.11;	Enzyme: Msp I}
    \item{\code{locus_23}}{Probe: J29;		Enzyme: Pu II}
  }
}
\details{
  SHARE algorithm requires subject-level information, i.e., it needs to
  know the haplotype/genotype sequences of every subjects in both case
  and control groups.  However, the original data in Kerem et
  al. (1989) only provide the sequence-level information, meaning that
  we only know what group (case/control) each haplotype sequence belongs
  to.  We need to simulate subject-level information to demostrate SHARE
  algorithm.  Two haplotypes with the same clinical status (having cystic fibrosis
  or not) are then ramdonly paired to form a psudo-subject with the that
  status.  
  
  Three objects will be attached after loading the dataset keremRand:

  The data.frame object \code{keremRandSeq} contains 186 sequences with
  23 SNPs.  The row names show the subject id and the sequence id within
  this subject.  The SNPs are coded as 1 referring to the large allele
  of the RFLP, and 2 referring to the smaller allele.

  The vector object \code{keremRandStatus} provides the CF/control
  status of each subject.  1 indicates subjects in case group (i.e.,
  CF), and 0 indicates control group. There are 47 subjects in CF group
  and 46 in control group.
  
  The data.frame object \code{keremRandAllele} contains allelic data for
  23 SNPs, coded as 0, 1, 2 as the number of minor alleles.

  How these three objects were created is shown in the example section.
}

\source{
  This dataset was originally released in Kerem et al. (1989), and was
  converted to R objects in Browning (2006). Browning's dataset could be found
  in the HapVLMC package
  (\url{http://www.stat.auckland.ac.nz/~browning/HapVLMC/index.htm}).
}
\references{
  S. R. Browning. Multilocus association mapping using variable-length markov chains.
  American Journal of Human Genetics, 78(6):903-913, Jun 2006.

  B. Kerem, J. M. Rommens, J. A. Buchanan, D. Markiewicz, T. K. Cox, A. Chakravarti,
  M. Buchwald, and L. C. Tsui. Identification of the cystic fibrosis gene: genetic analysis.
  Science (New York, N.Y.), 245(4922):1073-1080, Sep 8 1989.
}
\examples{
\dontrun{
## Here are how the psudo-subjects are simulated
#### loading HapVLMC package and the dataset
library(HapVLMC)
data(Kerem)
set.seed(20090313)
randOrder <- runif(nrow(kerem.snps.data))
keremRandSeq <- rbind(## randomly order the TRUE part
                   kerem.snps.data[kerem.status, ][order(randOrder[kerem.status]), ],
                   ## randomly order the FALSE part
                   kerem.snps.data[!kerem.status, ][order(randOrder[!kerem.status]), ]
                   )

nLoci <- ncol(keremRandSeq)
lociNum <- unlist(sapply(1:nLoci, 
           function(x){
		paste(paste(
			rep("0", ceiling(log10(nLoci)) - nchar(as.character(x))), collapse=""),
			x, sep="", collapse="")
                        })
                 )
colnames(keremRandSeq) <- paste("locus_", lociNum, sep="")

nSubj <- nrow(keremRandSeq)/2
subjNum <- unlist(sapply(1:nSubj, 
           function(x){
		paste(paste(
			rep("0", ceiling(log10(nSubj)) - nchar(as.character(x))), collapse=""),
			x, sep="", collapse="")
                        })
                 )
subjLabel <- paste("subj_", subjNum, sep="")
seqLabel <- paste("seq", 1:2, sep="_")
rownames(keremRandSeq) <- paste(rep(subjLabel, each=2), seqLabel, sep="_")

keremRandStatus <- c(rep(1, sum(kerem.status)/2), rep(0, sum(!kerem.status)/2))

keremRandAllele <- NULL
for(i in seq(1, nrow(keremRandSeq), by=2)){
    keremRandAllele <- rbind(keremRandAllele,
                       apply(keremRandSeq[c(i, i+1), ], 2,
                             function(x){
                                 ## counting how many small alleles
                                 sum(x==2)
                             }
                             )
                       )
}
rownames(keremRandAllele) <- unique(gsub("^(subj_.*)_seq_(.*)$", "\\1", rownames(keremRandSeq)))
}

## load keremRand
data(keremRand)

## check which objects are attached 
ls()

## dimention of psedu-subject data
dim(keremRandSeq)

## number of CF (TRUE) and control (FALSE) subjects
table(keremRandStatus)
}
\keyword{datasets}
