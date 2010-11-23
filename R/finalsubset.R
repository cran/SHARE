## Revision
##      2008-12-03
##              - calculate the z-score = coeff / s.e., where s.e. is at the var-covar matrix
##      2008-12-12
##              - conditionally output the result by the value of basesize
##      2008-12-23
##              - not using stop(); otherwise it will faile the creation of vignette
##              - add more info in the warning message
##              - add 'PACKAGE="SHARE"' into .C()
##      2009-03-30
##              - change "STE" to "Standard Error"
##              - add (Intercept) to hap1
##      2010-07-15
##              - add handling of non-genetic covariates 

finalsubset <- function(shareObj,
                        phenoData,
                        Minherit,
                        tol=1e-8,
                        verbose=0,
                        phase=0  
                        ){

    ## if not significant, report NULL for most of slots
    if (shareObj@bestsize==0) {
        ##pval<-0.051
        cat("Warning:\nIt is unable to determine the final subset if the best size is 0.\nNo result was outputted.\n")
    }else{ ## otherwise, find the final subset
        if (shareObj@bestsize > shareObj@maxSNP){
            outputsnp<-rep(0, 2 * shareObj@maxSNP - shareObj@bestsize)
        }else{
            outputsnp<-rep(0, shareObj@bestsize)
        }

        varstore <- rep(0,1000)
        coef <- rep(0,100)
        outnhaps=0
        outhapvec <- rep(0,1000)
        outhapfreq<- rep(0,100)

        temp1 <- .C("finalsubset",
                    id=as.integer(shareObj@poolHapPair),
                    nSubj=as.integer(length(shareObj@nPosHapPair)),
                    nObs=as.integer(length(shareObj@poolHapPair)),
                    nPair=as.integer(shareObj@nPosHapPair),
                    ccStatus=as.integer(rep(phenoData, shareObj@nPosHapPair)),
                    nLoci=as.integer(ncol(shareObj@haploSeq)),
                    nHap=as.integer(length(shareObj@haploFreq)),
                    hap1=as.integer(shareObj@hap1),
                    hap2=as.integer(shareObj@hap2),

                    isUniHap=as.integer(shareObj@uhap),
                    hapFreq=as.numeric(shareObj@haploFreq),
                    postHap=as.numeric(shareObj@post),
                    nngcov=as.integer(shareObj@nngcov),
                    ngcov_vec=as.numeric(shareObj@ngcov),

                    maxSNP=as.integer(shareObj@maxSNP),
                    as.integer(shareObj@bestsize),
                    outputsnp=as.integer(outputsnp),
                    outhapvec=as.integer(outhapvec),
                    outhapfreq=as.numeric(outhapfreq),

                    as.numeric(tol),
                    as.integer(phase),
                    varstore=as.double(varstore),
                    coef=as.double(coef),
                    outnhaps=as.integer(outnhaps),
                    Minherit=as.integer(Minherit),
                    as.integer(verbose),
                    PACKAGE="SHARE"
                    )

        ## var-covar matrix
        amat <- matrix(temp1$varstore[1:((temp1$outnhaps+temp1$nngcov) ^2)], 
							  nrow=temp1$outnhaps+temp1$nngcov,
                       ncol=temp1$outnhaps+temp1$nngcov, 
							  byrow=TRUE)
        bmat <- matrix(temp1$varstore[( ( (temp1$outnhaps+temp1$nngcov)^2 )+1):(2*( (temp1$outnhaps+temp1$nngcov) ^2))],
                       nrow=temp1$outnhaps+temp1$nngcov, ncol=temp1$outnhaps+temp1$nngcov, byrow=TRUE)


        require(MASS) ## for ginv function
        var<-ginv(amat)%*% bmat %*% ginv(amat)
        
        pval <- drop(
			1-pchisq(t(temp1$coef[2:(temp1$outnhaps+temp1$nngcov)]) %*% ginv(var[2:(temp1$outnhaps+temp1$nngcov),2:(temp1$outnhaps+temp1$nngcov)]) %*% temp1$coef[2:(temp1$outnhaps+temp1$nngcov)],(temp1$outnhaps+temp1$nngcov)-1))

        ## name of haplotypes, SNPs and non-genetic covariates
        snpNames <- colnames(shareObj@haploSeq)[sort(temp1$outputsnp)]
        hapName <- paste("hap", 1:temp1$outnhaps, sep="")
        ngcovName <- paste("ngCov", 1:temp1$nngcov, sep="")
        ## add (Intercept) to hap1
        ## hapName[1] <- paste(hapName[1], "(Intercept)")
        
        ## data frame for the output haplotypes
        hapDF <- data.frame(matrix(temp1$outhapvec[1:(length(temp1$outputsnp)*temp1$outnhaps)],
                                   nrow=temp1$outnhaps, ncol=length(temp1$outputsnp), byrow=TRUE)
                            )
        colnames(hapDF) <- snpNames
        rownames(hapDF) <- hapName

        ## haplotype frequency
        hapFreq <- temp1$outhapfreq[1:temp1$outnhaps]
        names(hapFreq) <- hapName

        ## coefficients (haplotypes  + non-genetic covariates)
        coef <- temp1$coef[1:(temp1$outnhaps+temp1$nngcov)]
        if(temp1$nngcov!=0) {names(coef) <- c(hapName, ngcovName)} else {names(coef) <- hapName}

        ## currently not returning non-genetic covariates 
        ## returning matrix of repeated values for subjects (ie, not original input covariates) 
      	#  if(temp1$nngcov!=0) {
			#  	ngcovMat <- data.frame(matrix(shareObj@ngcov, nrow=length(shareObj@poolHapPair), ncol=shareObj@nngcov, byrow=FALSE))
      	#  names(ngcovMat) <- ngcovName 
			#	} else { ngcovMat <- data.frame() }

        ## s.e.
        stdErr <- sqrt(diag(var))

        ## z-score
        zScore <- coef / stdErr

        ## p-values
        pValue <- 1-pnorm(abs(zScore))

        ## GLM-like output

        haploTest <- data.frame(coef, stdErr, zScore, pValue)
        colnames(haploTest) <- c("Coefficient", "Standard Error", "z-Score", "p-Value")
        hapbase <- which.max(temp1$outhapfreq[1:temp1$outnhaps])
        hapName1 <- paste("hap",c(hapbase,(1:temp1$outnhaps)[-hapbase]))
        hapName1[1] <- paste(hapName1[1], "(Intercept)")
        if(temp1$nngcov!=0) {rownames(haploTest) <- c(hapName1, ngcovName)} else {rownames(haploTest) <- hapName1} 

        outObj <- shareObj
        class(outObj) <- "share"
        outObj@finalHap <- hapDF
        outObj@finalHapFreq <- hapFreq
        outObj@finalHapTest <- haploTest
        outObj@globalP <- pval
        outObj@modelmethod <- "Cross-Val"
        outObj@inherit <- c("additive", "dominant", "recessive")[Minherit]

        return(outObj)
    }
    
}
