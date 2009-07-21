setClass("haplo",
         contains="haploSet",         
         representation(##hap == haploSeq @ haploSet class
                        ##hap="data.frame",       ## unique haplotype X loci (haplo.em$haplotype)
                        ## hapProb == haploFreq @ haploSet class
                        ##hapProb="vector",       ## haplotype frequencies for each unique haplotype (haplo.em$hap.prob)
                        hap1="vector",          ## index of first haplotype of each subject (haplo.em$hap1code)
                        hap2="vector",          ## index of second haplotype of each subject (haplo.em$hap2code)

                        poolHapPair="vector",   ## pooled haplotype ID (haplo.em$subj.id)
                        nPosHapPair="table",    ## number of possible haplotype pairs (haplo.em$nreps)

                        post="vector",          ## posterior probability for the pooled haplotype pairs
                        
                        pheno="data.frame"      ## pheno type data after re-ordering
                        )
         )
## show
setMethod("show", signature=signature(object="haplo"),
          function(object){
              cat(paste("An object of class \"haplo\" (extended from class \"haploSet\")",
                        "\n# of SNP:", nSNP(object),
                        "\n# of haplotype sequence: ", nrow(object@haploSeq),
                        "\n"
                        )
                  )
          }
          )
