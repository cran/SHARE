## Revision
##      2009-03-30
##              - add new slot "inherit" to store the value of "Minherit"
## Revision
##      2010-07-10
##              - add new slots "nngcov" and "ngcov_vec" for non-genetic covariate information
##                nngcov = number of non-genetic covariates, 
##						ngcov_vec = vector of all covariate values

setClass("share",
         representation(uhap="vector",
                        weight="vector",
                        nFold="numeric",
                        maxSNP="numeric",
                        deviance="vector",
                        bestsize="numeric",
                        finalHap="data.frame",
                        finalHapFreq="vector",
                        finalHapTest="data.frame",
                        globalP="numeric",
                        modelmethod="character",
                        inherit="character",
                        nngcov="numeric",
                        ngcov="vector"
                        ),
         contains="haplo"
         )


##FIXME:  works, but mask graphics::plot
##        need to use namespace to prevent the conflict with graphics::plot
##
##setGeneric("plot",
##           function(.Object, ...){
##               standardGeneric("plot")
##               }
##           )
##
##setMethod("plot", signature("share"),
##          function(.Object, ...){
##              graphics::plot(.Object@deviance,
##                   type="b", xaxt="n",
##                   xlab="Model Size", ylab="Deviance")
##              axis(side=1, at=1:(length(.Object@deviance)),
##                   labels=as.character(c(0, 1:.Object@maxSNP, (.Object@maxSNP-1):1 )))
##              }
##          )

## show
setMethod("show", signature=signature(object="share"),
          function(object){

              if (object@modelmethod=="Cross-Val") {
                  cat(paste("An object of class \"share\"",
                            "\nGiven nFold=", object@nFold, " cross-validation", " and maxSNP=", object@maxSNP, 
                            ", the best number of SNP to select is ", object@bestsize, ".",
                            "\n",
                            sep=""
                            )
                      )

                  if (object@bestsize>0) {
                      testDF <- object@finalHapTest
                      for(col in colnames(testDF)){
                          testDF[, col] <- round(testDF[, col], 6)
                      }
                      ##cat("An object of class \"finalsubset\"")
                      cat(paste("\n# of selected SNP:", ncol(object@finalHap),
                                "\n# of haplotype sequence: ", nrow(object@finalHap),
                                "\n"
                                )
                          )
                      cat("\nHaplotype sequences in the final subset:\n")
                      print(cbind(object@finalHap, Frequency=round(object@finalHapFreq, 6)))
                      cat("\nTesting:\n")
                      print(testDF)
                      cat(paste("\nGlobal p-value:", round(object@globalP, 6), "\n"))
                  }
              } ## end of if (object@modelmethod=="Cross-Val") 
              else if (object@modelmethod=="BIC") {
                  cat(paste("An object of class \"share\"",
                            "\nGiven BIC criterion", " and maxSNP=", object@maxSNP, 
                            ", the best number of SNP to select is ", object@bestsize, ".",
                            "\n",
                            sep=""
                            )
                      )

                  if (object@bestsize>0) {
                      testDF <- object@finalHapTest
                      for(col in colnames(testDF)){
                          testDF[, col] <- round(testDF[, col], 6)
                      }
                      ##cat("An object of class \"finalsubset\"")
                      cat(paste("\n# of selected SNP:", ncol(object@finalHap),
                                "\n# of haplotype sequence: ", nrow(object@finalHap),
                                "\n"
                                )
                          )
                      cat("\nHaplotype sequences in the final subset:\n")
                      print(cbind(object@finalHap, Frequency=round(object@finalHapFreq, 6)))
                      cat("\nTesting:\n")
                      print(testDF)
                      cat(paste("\nGlobal p-value:", round(object@globalP, 6), "\n"))
                  }
              } ## end of if (object@modelmethod=="BIC")
          }
          )       

setGeneric("dplot",
           function(shareObj, ...){
               standardGeneric("dplot")
           })

setMethod("dplot", signature=c("share"),
          function(shareObj, ...){
              if (shareObj@modelmethod=="Cross-Val") {   
                  plot(shareObj@deviance,
                       type="b", xaxt="n",
                       xlab="Model Size", ylab="Prediction Deviance", ...)
                  axis(side=1, at=1:(length(shareObj@deviance)),
                       labels=as.character(c(0, 1:shareObj@maxSNP, (shareObj@maxSNP-1):1 )))
              } else if (shareObj@modelmethod=="BIC") {
                  plot(shareObj@deviance,
                       type="b", xaxt="n",
                       xlab="Model Size", ylab="BIC", ...)
                  axis(side=1, at=1:(length(shareObj@deviance)),
                       labels=as.character(c(0, 1:shareObj@maxSNP, (shareObj@maxSNP-1):1 )))
              }
          }
          )

## FIXME:  re-declare summary will lost the old definition
##setGeneric("summary",
##           function(.Object, ...)
##           standardGeneric("summary")
##           )
##
##
##setMethod("summary", "finalsubset",
##          function(.Object){
##              glmout <- data.frame(round(.Object@finalHapCoef, 6),
##                                   c(round(.Object@expPValue, 6), rep("", length(.Object@finalHapCoef)-1))
##                                   )
##              colnames(glmout) <- c("Coefficient", "p-Value")
##              rownames(glmout) <- rownames(.Object@finalHap)
##              return(glmout)
##              }
##          )
