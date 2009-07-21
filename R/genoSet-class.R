## Code:        genoSet-class.R
## Purpose:     definition of genoSet class and associated methods
##      
## Revision
##      2009-03-10
##              - move from S4-class.R & haplo.R

setClass("genoSet",
         representation(genoSeq="data.frame",       ## k genotype sequences (0-1-2 format) * n SNPs
                        phenoData="data.frame"  ## k genotype sequences * p variables (phenotype)
                        )
         )

## show
setMethod("show", signature=signature(object="genoSet"),
          function(object){
              cat(paste("An object of class \"genoSet\"",
                        "\n# of SNP:", nSNP(object),
                        "\n# of genotype sequence: ", nSeq(object),
                        "\nPhenotype:\n", paste(colnames(object@phenoData), sep=", ", collapse=""),
                        "\n"
                        )
                  )
          }
          )

## accessor
setGeneric("genoSeq",
           function(.Object){
               standardGeneric("genoSeq")
           }
           )
setMethod("genoSeq", signature=signature(.Object="genoSet"),
          function(.Object){
              .Object@genoSeq
          }
          )

setGeneric("phenoData",
           function(.Object){
               standardGeneric("phenoData")
           }
           )
setMethod("phenoData", signature=signature(.Object="genoSet"),
          function(.Object){
              .Object@phenoData
          }
          )

## general methods
setMethod("nSNP", signature=signature(.Object="genoSet"),
          function(.Object){
              ncol(.Object@genoSeq)
          }
          )

setMethod("nameSNP", signature=signature(.Object="genoSet"),
          function(.Object){
              colnames(.Object@genoSeq)
          }
          )

setMethod("nSeq", signature=signature(.Object="genoSet"),
          function(.Object){
              nrow(.Object@genoSeq)
          }
          )

setMethod("nameSeq", signature=signature(.Object="genoSet"),
          function(.Object){
              rownames(.Object@genoSeq)
          }
          )

## specific methods
setGeneric("haplo",
           function(.Object, ...){
               standardGeneric("haplo")
           }
           )

setMethod("haplo", signature=signature(.Object="genoSet"),
          function(.Object, wgt=NULL, ...){
              #### check if wgt is a vector of length nrow(dat)
              if(is.null(wgt)){
                  wgt <- rep(1, nrow(.Object@genoSeq))
              }else{
                  wgtTry <- try(is.vector(wgt))
                  if(inherits(wgtTry, "try-error")){
                      stop("The argument wgt must be vector or NULL.")
                  }
                  if(length(wgt)!=nrow(.Object@genoSeq)){
                      stop("The length of argument wgt must equal to nrow(dat).")
                  }
              }

              ## randomize the order
              set.seed(37273)
              newOrder <- sample(1:nrow(.Object@genoSeq))
              .Object@genoSeq <- .Object@genoSeq[newOrder, ]
              if(ncol(.Object@phenoData)==1){
                  newPheno <- data.frame(.Object@phenoData[newOrder, ])
                  colnames(newPheno) <- colnames(.Object@phenoData)
                  .Object@phenoData <- newPheno
              }else{
                  .Object@phenoData <- data.frame(.Object@phenoData)[newOrder, ]
              }
              
              require(haplo.stats)
              ## construct genotype data.frame
              geno <- geno1to2(genoSeq(.Object), locus.label=nameSNP(.Object))
              
              ## haplotype reconstruction & frequency estimation
              em <- haplo.em(geno=geno, locus.label=nameSNP(.Object), 
                             miss.val=NA, weight=wgt,
                             control=haplo.em.control(loci.insert.order=NULL, insert.batch.size=3, min.posterior=0.01),
                             ...
                             )

              ## return the required values to a haplo-class object
              return(new("haplo",
                         haploSeq=em$haplotype,
                         haploFreq=em$hap.prob,
                         hap1=em$hap1code,
                         hap2=em$hap2code,
                         poolHapPair=em$subj.id,
                         nPosHapPair=em$nreps,
                         post=em$post,
                         pheno=.Object@phenoData
                         )
                     )              

          }
          )

