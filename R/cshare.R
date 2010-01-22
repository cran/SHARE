## Revision
##      2008-12-12
##              - add pheno data into argument
##              - change the sign of deviance
##              - add bestsize to the output
##      2008-12-23
##              - add 'PACKAGE="SHARE"' into .C()
##      2009-03-03
##              - James D: add two arguments: ModSelMethod & Minherit
##      2009-03-09
##              - change the acceptable value of ModSelMethod to "Cross-Val" and "BIC"
##              - change the acceptable value of Minherit to "additive", "dominant", and "recessive"
##      2009-03-30
##              - report standard error in the testing section
##              - add (Intercept) to hap1
##              - in shareTest, remove arguments "nfold", "maxsnps", "ModSelMethod", because they are included in outObj already
##              - compare global p-values only if bestsize != 0
##      2009-05-15
##              - make the default value of verbose=FALSE
##      2009-07-02
##              - use "status" argument to determine which clinical status in phenoData for analysis 

cshare <- function(haploObj,    ## haploObj
                   status,      ## label of status in the haploObj@pheno data.frame
                   nfold=10,    ## number of fold (numeric), NOT USED IF MODEL SELECTION IS THROUGH BIC
                   maxsnps,     ## maximum SNP selected in the model (numeric) 
                   tol=1e-8,
                   verbose=FALSE,   ## TRUE:  turn on log file
	           ModSelMethod=c("Cross-Val", "BIC"),    ## default model selection method is cross-val (=1), however can use BIC (=2) for faster computation
                   Minherit=c("additive", "dominant", "recessive")          ## the default mode of inheritance is "additive",  the other options is "dominant" and  "recessive" 
                   ){

    ## checking argument
    ModSelMethod <- match.arg(ModSelMethod)
    Minherit <- match.arg(Minherit)
    
    ## Validate arguments 
    ## check if dat is a haplo object
    if(class(haploObj)!="haplo"){
        stop("The argument dat must be object of the class 'haplo'.")
    }

    ## check if nfold is a numberic number
    nfoldTry <- try(is.numeric(nfold))
    if(inherits(nfoldTry, "try-error")){
        stop("The argument nfold must be numeric.")
    }
    if(length(nfold)!=1){
        stop("The length of argument nfold must be 1.")
    }

    ## check if maxsnp is a numberic number
    maxsnpsTry <- try(is.numeric(maxsnps))
    if(inherits(maxsnpsTry, "try-error")){
        stop("The argument maxsnps must be numeric.")
    }
    if(length(maxsnps)!=1){
        stop("The length of argument maxsnps must be 1.")
    }
    if(maxsnps>(ncol(haploObj@haploSeq)-1)){
        stop("The value of argument maxsnps must be less than the total number of SNPs.")
    }

    uhap <- rep(0, length(haploObj@haploFreq)* ncol(haploObj@haploSeq)) 
    count<-1
    for (i in 1:(length(haploObj@haploFreq))){
        for (j in 1:(ncol(haploObj@haploSeq))) {      
            uhap[count]<- as.numeric(haploObj@haploSeq[i,j])
            count<-count+1
        }
    }

    deviance <- rep(0,2*maxsnps)
    phase <- 0

    if (ModSelMethod=="Cross-Val") {
        temp <- .C("xshare",
                   id=as.integer(haploObj@poolHapPair),
                   nSubj=as.integer(length(haploObj@nPosHapPair)),
                   nObs=as.integer(length(haploObj@poolHapPair)),
                   nPair=as.integer(haploObj@nPosHapPair),
                   ccStatus=as.integer(rep(unlist(haploObj@pheno[, status]), haploObj@nPosHapPair)),
                   nLoci=as.integer(ncol(haploObj@haploSeq)),
                   nFold=as.integer(nfold),
                   nHap=as.integer(length(haploObj@haploFreq)),
                   hap1=as.integer(haploObj@hap1),
                   hap2=as.integer(haploObj@hap2),               
                   uhap=as.integer(uhap),
                   hapFreq=as.numeric(haploObj@haploFreq),
                   postHap=as.numeric(haploObj@post),
                   maxSNP=as.integer(maxsnps),
                   deviance=as.numeric(deviance),
                   converge=as.numeric(tol),
                   verbose=as.integer(sum(verbose)),
                   phase=as.integer(phase),
                   Minherit=match(Minherit, c("additive", "dominant", "recessive")),
                   PACKAGE="SHARE"
                   )

        outObj <- haploObj
        class(outObj) <- "share"
        outObj@uhap=uhap
        outObj@weight=temp$postHap
        outObj@nFold=nfold
        outObj@maxSNP=maxsnps
        outObj@deviance= - temp$deviance
        outObj@bestsize=c(0, 1:maxsnps, (maxsnps-1):1 )[which.min(outObj@deviance)]

        if (outObj@bestsize==0) {
            fSubset <- outObj
            fSubset@modelmethod <- "Cross-Val"
        } else {
            fSubset <- finalsubset(outObj, phenoData=unlist(haploObj@pheno[, status]), match(Minherit, c("additive", "dominant", "recessive")))
        }
        
        return(fSubset)
	
    } ## end of (ModSelMethod=="Cross-Val")
    else if (ModSelMethod=="BIC") {
	alpha <- 2*log(length(haploObj@nPosHapPair))
	outputsnp <- rep(0,ncol(haploObj@haploSeq))
	outhapvec <- rep(0,1000)
	outhapfreq <- rep(0,100)
	varstore <- rep(0,1000)
	coef <- rep(0,100)
	outnhaps <- 0
        bestsize <- 0
	deviance <- rep(0,2*maxsnps)
        phase <- 0
	temp <- .C("stepwise_search_alpha", 
                   id=as.integer(haploObj@poolHapPair),
                   nSubj=as.integer(length(haploObj@nPosHapPair)),
                   nObs=as.integer(length(haploObj@poolHapPair)),
                   nPair=as.integer(haploObj@nPosHapPair),
                   postHap=as.numeric(haploObj@post),
                   ccStatus=as.integer(rep(unlist(haploObj@pheno[, status]), haploObj@nPosHapPair)),
                   nLoci=as.integer(ncol(haploObj@haploSeq)),
                   nHap=as.integer(length(haploObj@haploFreq)),
                   hap1=as.integer(haploObj@hap1),
                   hap2=as.integer(haploObj@hap2),
                   uhap=as.integer(uhap),
                   hapFreq=as.numeric(haploObj@haploFreq),
                   bestsize=as.integer(bestsize),
                   Minherit=match(Minherit, c("additive", "dominant", "recessive")),
                   deviance=as.numeric(deviance),
                   maxSNP=as.integer(maxsnps),
                   converge=as.numeric(tol),
                   alpha=as.numeric(alpha),
                   verbose=as.integer(sum(verbose)),
                   phase=as.integer(phase), 
                   outputsnp=as.integer(outputsnp), 
                   outhapvec=as.integer(outhapvec), 
                   outhapfreq=as.numeric(outhapfreq), 
                   varstore=as.numeric(varstore), 
                   coef=as.numeric(coef), 
                   outnhaps=as.integer(outnhaps),
                   PACKAGE="SHARE"
                   )

        outObj <- haploObj
        class(outObj) <- "share"
        outObj@uhap=uhap
        outObj@weight=temp$postHap
        outObj@nFold=nfold
        outObj@maxSNP=maxsnps
        outObj@deviance= temp$deviance
        outObj@bestsize=temp$bestsize
	
        if (outObj@bestsize==0) {
            outObj2 <- outObj
            outObj2@modelmethod <- "BIC"
            
        } else {
            
	    ## var-covar matrix
            amat <- matrix(temp$varstore[1:(temp$outnhaps ^2)], nrow=temp$outnhaps,
                           ncol=temp$outnhaps, byrow=TRUE)
            bmat <- matrix(temp$varstore[((temp$outnhaps ^2)+1):(2*(temp$outnhaps ^2))],
                           nrow=temp$outnhaps, ncol=temp$outnhaps, byrow=TRUE)

            require(MASS) ## for ginv funciton
            var<-ginv(amat)%*% bmat %*% ginv(amat)
            
            pval <- drop(1-pchisq(t(temp$coef[2:temp$outnhaps]) %*% ginv(var[2:temp$outnhaps,2:temp$outnhaps]) %*% temp$coef[2:temp$outnhaps],temp$outnhaps-1))

            ## name of haplotypes & SNPs
            snpNames <- colnames(outObj@haploSeq)[sort(temp$outputsnp)]
            hapName <- paste("hap", 1:temp$outnhaps, sep="")
            ## add (Intercept) to hap1
            ##hapName[1] <- paste(hapName[1], "(Intercept)")
            
            ## data frame for the output haplotypes
            hapDF <- data.frame(matrix(temp$outhapvec[1:(temp$bestsize*temp$outnhaps)],
                                       nrow=temp$outnhaps, ncol=temp$bestsize, byrow=TRUE)
                                )
            colnames(hapDF) <- snpNames
            rownames(hapDF) <- hapName

            ## haplotype frequency
            hapFreq <- temp$outhapfreq[1:temp$outnhaps]
            names(hapFreq) <- hapName

            ## coefficient
            coef <- temp$coef[1:temp$outnhaps]
            names(coef) <- hapName

            ## s.e.
            stdErr <- sqrt(diag(var))

            ## z-score
            zScore <- coef / stdErr

            ## p-values
            pValue <- 1-pnorm(abs(zScore))

            ## GLM-like output
            haploTest <- data.frame(coef, stdErr, zScore, pValue)
            colnames(haploTest) <- c("Coefficient", "Standard Error", "z-Score", "p-Value")
            hapbase <- which.max(temp$outhapfreq[1:temp$outnhaps])
            hapName1 <- paste("hap",c(hapbase,(1:temp$outnhaps)[-hapbase]))
            hapName1[1] <- paste(hapName1[1], "(Intercept)")
            rownames(haploTest) <- hapName1

            outObj2 <- outObj
            outObj2@modelmethod <- "BIC"
            outObj2@finalHap <- hapDF
            outObj2@finalHapFreq <- hapFreq
            outObj2@finalHapTest <- haploTest
            outObj2@globalP <- pval
            outObj2@inherit <- c("additive", "dominant", "recessive")[Minherit]
        }  
        
        return(outObj2)
    } ## end of (ModSelMethod=="BIC")
}

shareTest <- function(outObj,            ## the returned object from chare function
                      haploObj,          ## haploObj
                      status,
                      tol=1e-8,
                      verbose=FALSE,     ## TRUE: turn on log file
                      nperm=1000         ## the number of permutations 
                      ) {
    phase <- 0

    count <- 0
    for (i in 1:nperm) {
        set.seed(38329832+i)
        if (i%%10==0) cat(i, "... ")
        haploObj@pheno[, "y"] <- sample(as.vector(haploObj@pheno[, status]))
        out <- cshare(haploObj, "y", outObj@nFold, outObj@maxSNP, tol, verbose,
                      outObj@modelmethod, outObj@inherit)

        if(out@bestsize!=0){
            if (outObj@globalP > out@globalP) count <- count+1
        }
        pval <- count/(i+1)
        if (count>=50) {
            cat("Permutation test stops early due to insignificance.","\n")
            break
        }
    }
    return(pval)
}
