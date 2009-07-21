setClass("haploSet",
         representation(haploSeq="data.frame",      ## m haplotype sequences * n SNPs
                        haploFreq="vector"      ## m haplotype sequences * 1
                        )
         )

## show
setMethod("show", signature=signature(object="haploSet"),
          function(object){
              cat(paste("An object of class \"haploSet\"",
                        "\n# of SNP:", nSNP(object),
                        "\n# of haplotype sequence: ", nSeq(object),
                        "\n"
                        )
                  )
          }
          )

## accessor
setGeneric("haploSeq",
           function(.Object){
               standardGeneric("haploSeq")
           }
           )
setMethod("haploSeq", signature=signature(.Object="haploSet"),
          function(.Object){
              .Object@haploSeq
          }
          )

setGeneric("haploFreq",
           function(.Object){
               standardGeneric("haploFreq")
           }
           )
setMethod("haploFreq", signature=signature(.Object="haploSet"),
          function(.Object){
              .Object@haploFreq
          }
          )

## general methods
setMethod("nSNP", signature=signature(.Object="haploSet"),
          function(.Object){
              ncol(.Object@haploSeq)
          }
          )

setMethod("nameSNP", signature=signature(.Object="haploSet"),
          function(.Object){
              colnames(.Object@haploSeq)
          }
          )

setMethod("nSeq", signature=signature(.Object="haploSet"),
          function(.Object){
              nrow(.Object@haploSeq)
          }
          )

setMethod("nameSeq", signature=signature(.Object="haploSet"),
          function(.Object){
              rownames(.Object@haploSeq)
          }
          )

