## Code:        S4-class.R
## Purpose:     Define general methods
## Revision
##      2008-12-03
##              - define new class: genoSet & haploSet
##      2008-12-12
##              - remove pheno data from haplo-class
##      2008-12-23
##              - rownames() instead of rowname()
##      2009-03-10
##              - combine "finalsubset" class into "share" class
##              - move class-related definitions to separate code


##---Generic Function of General Methods ---##

##--- for SNP (loci) info ---##
## nSNP
setGeneric("nSNP",
           function(.Object){
               standardGeneric("nSNP")
           }
           )

## nameSNP
setGeneric("nameSNP",
           function(.Object){
               standardGeneric("nameSNP")
           }
           )


##--- for Sequence (genotype of haplotype) info ---##
## nSeq
setGeneric("nSeq",
           function(.Object){
               standardGeneric("nSeq")
           }
           )

## nameSeq
setGeneric("nameSeq",
           function(.Object){
               standardGeneric("nameSeq")
           }
           )

