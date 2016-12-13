meta.abf<-function(betas,ses,prior.sigma,prior.cor="indep",prior.rho=NA,cryptic.cor=NA,log=FALSE,log10=FALSE,tolerance=1e-1000){
    #betas is a vector, matrix, or dataframe of observed effect sizes of a single SNP in a set of studies. If betas is a vector, then there is a single SNP and each element of the vector is assumed to correspond to a study. If betas is a matrix, then the rows should be SNPs and the columns should be studies.
    #ses is a vector, matrix, or dataframe of standard errors corresponding to those in betas. It should have the same dimension of betas.
    #prior.sigma is the prior on true effect sizes for each SNP in each study. It can be a flat value, set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as betas).
    #prior.cor is a square matrix whose row and column numbers are the same as the number of studies. Its elements are the pairwise correlations between true effect sizes of the studies. It can take values "indep" (independent effects), "fixed" (fixed effects), "correlated" (correlated effects, which requires the prior.rho parameter to be set), as well as individual matrices. If betas and ses are matrices, the same prior.cor will be applied to every row (representing every SNP).
    #prior.rho is either a single value or the upper triangle of a correlation matrix for the prior.cor matrix when it is set to "correlated". If this value is set, but prior.cor is not set to "correlated", this parameter will be ignored.
    #cryptic.cor a square matrix whose row and coumn numbers are the same as the number of studies. If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.
    #log sets whether the answer should be given in log space.
    #log10 sets whether the answer should be given in log10 space.
    #tolerance for the ABF calculation, this can be lowered (or raised, if necessary) if the answers are not what was expected. Should probably never be altered, but is there in case it is needed.


    #Initial error checking
    if(!class(betas) %in% c("data.frame","matrix","numeric")){
        stop("betas must be a data frame, a matrix, or a numeric vector.")
    }
    if(!class(ses) %in% c("data.frame","matrix","numeric")){
        stop("ses must be a data frame, a matrix, or a numeric vector.")
    }
    if(!class(prior.sigma) %in% c("data.frame","matrix","numeric")){
        stop("prior.sigma must be a data frame, a matrix, a numeric vector, or single numeric value.")
    }
    if(log && log10){
        stop("can give the approximate Bayes factor in log space or log_10 space, but not both at the same time.")
    }
    if(tolerance>sqrt(.Machine$double.eps)){
        warning(paste0("Your tolerance might be too high. The standard value for the internal functions is ",sqrt(.Machine$double.eps),"."))
    }

    #Coerces data frame data to matrix
    if(class(betas) %in% c("data.frame","matrix")){
        if(all(apply(betas,2,class)=="numeric")){
            betas<-as.matrix(betas)
        } else {
            stop("elements of \"betas\" are not numeric.")
        }
    }
    if(class(ses) %in% c("data.frame","matrix")){
        if(all(apply(ses,2,class)=="numeric")){
            ses<-as.matrix(ses)
        } else {
            stop("elements of \"ses\" are not numeric.")
        }
    }
    if(class(prior.sigma) %in% c("data.frame","matrix")){
        if(all(apply(prior.sigma,2,class)=="numeric")){
            prior.sigma<-as.matrix(prior.sigma)
        } else {
            stop("elements of \"prior.sigma\" are not numeric.")
        }
    }
    if(class(betas)!=class(ses)){
        stop(paste0("betas belongs to class ",class(betas),", but ses belongs to class ",class(ses),". These two sets of values must belong to the same class."))
    }

    #Get the number of studies
    if(class(betas)=="matrix"){
        if(dim(betas)[1]!=dim(ses)[1] || dim(betas)[2]!=dim(ses)[2]){
            stop("betas and ses do not have the same dimensions.")
        }
        nstudies<-ncol(betas)
    } else {
        if(length(betas)!=length(ses)){
            stop("betas and ses do not have the same length.")
        }
        nstudies<-length(betas)
    }

    ##Check prior.sigma isn't erroneous.
    if(class(prior.sigma)=="matrix"){
        if(class(betas)!="matrix"){
            stop("prior.sigma must be a vector if betas and ses are.")
        }
        if(dim(prior.sigma)[1]!=dim(betas)[1] || dim(prior.sigma)[2]!=dim(betas)[2]){
            stop("the dimensions of prior.sigma do not match the dimensions of betas.")
        }
    } else if(!length(prior.sigma) %in% c(1,nstudies)){
        stop("prior.sigma should either be a single value or a vector whose length is equal to the number of studies.")
    }
    if(any(prior.sigma<=0)){
        stop("all values of prior.sigma must be > 0.")
    }

    ##Get the prior correlation matrix.
    if(class(prior.cor)=="character"){
        if(prior.cor=="correlated"){
            if(is.na(prior.rho)){
                stop("prior.rho must be set or prior.cor must be changed to something other than \"correlated.\"")
            }
            if(!is.numeric(prior.rho)){
                stop("prior.rho must be numeric.")
            }
            if(!length(prior.rho) %in% c(1,choose(nstudies,2))){
                stop("prior.rho should be either a flat value or the full upper triangle of the desired correlation matrix.")
            }
            if(length(prior.cor)==1){
                prior.cor.mat<-matrix(prior.rho,nrow=nstudies,ncol=nstudies)
                diag(prior.cor.mat)<-1
            } else {
                prior.cor.mat<-diag(nstudies)
                prior.cor.mat[upper.tri(prior.cor.mat,diag=FALSE)]<-prior.cor
                prior.cor.mat[lower.tri(prior.cor.mat)]<-t(prior.cor.mat)[lower.tri(prior.cor.mat)]
            }
        } else if(prior.cor=="indep"){
            prior.cor.mat<-diag(nstudies)
        } else if(prior.cor=="fixed"){
            prior.cor.mat<-matrix(1,nrow=nstudies,ncol=nstudies)
        }
    } else if(class(prior.cor)!="matrix"){
        stop("prior.cor should be a matrix, or one of the following values: \"indep\", \"fixed\", or \"correlated\"")
    } else {
        if(dim(prior.cor)[1]!=dim(prior.cor)[2]){
            stop("prior.cor is not a square matrix.")
        }
        if(dim(prior.cor[1])!=nstudies){
            stop("the dimensions of prior.cor do not match the number of studies in the meta-analysis.")
        }
        if(!isSymmetric.matrix(prior.cor)){
            stop("prior.cor is not a symmetric matrix.")
        }
        if(any(eigen(prior.cor)$values<0)){
            stop("prior.cor is not positive semidefinite.")
        }
        if(!all(diag(prior.cor) %in% c(0,1))){
            stop("the diagonal of prior.cor should be 1 for all studies with true effects and 0 elsewhere.")
        }
        if(all(prior.cor[which(diag(prior.cor==0)),]!=0) || all(prior.cor[,which(diag(prior.cor==0))]!=0)){
            if(length(which(diag(prior.cor==0)))==1){
                stop(paste0("Row ",which(diag(prior.cor==0))," and column ",which(diag(prior.cor==0))," of prior.cor should all be 0, or else prior.cor[",which(diag(prior.cor==0)),",",which(diag(prior.cor==0)),"] should be 1."))
            } else {
                stop(paste0("Rows ",paste(which(diag(prior.cor==0)),collapse=",")," and columns ",paste(which(diag(prior.cor==0)),collapse=",")," of prior.cor should all be 0, or else prior.cor[c(",paste(which(diag(prior.cor==0)),collapse=","),"), c(",paste(which(diag(prior.cor==0)),collapse=","),")] should all be 1."))
            }
        }
        prior.cor.mat<-prior.cor
    }

    ##Get the cryptic correlation matrix
    if(all(is.na(cryptic.cor)) && length(cryptic.cor)==1){
        cryptic.cor.mat<-diag(nstudies)
    } else if(class(cryptic.cor)!="matrix"){
        stop("cryptic.cor should be a square matrix.")
    } else if(dim(cryptic.cor)[1]!=dim(cryptic.cor)[2]){
        stop("cryptic.cor should be a square matrix.")
    } else if(dim(cryptic.cor)[1]!=nstudies){
        stop("the dimensions of cryptic.cor do not match the number of studies in the meta-analysis.")
    } else if(!isSymmetric.matrix(cryptic.cor)){
        stop("cryptic.cor is not a symmetric matrix.")
    } else if(any(eigen(cryptic.cor)$values<0)){
        stop("cryptic.cor is not positive definite.")
    } else if(any(diag(cryptic.cor)!=1)){
        stop("the diagonal of cryptic.cor should be 1 uniformly.")
    } else if(length(cryptic.cor)>1 && any(is.na(cryptic.cor))){
        stop("If cryptic.cor is defined, then no NAs are permitted.")
    } else {
        cryptic.cor.mat<-cryptic.cor
    }


    ##Internal function to calculate the multivariate ABF.
    get.abf<-function(info,pcm,ccm,log=FALSE,log10=FALSE){
        m<-length(info)/3
        b<-info[1:m]
        se<-info[(m+1):(2*m)]
        prior.se<-info[((2*m)+1):(3*m)]

        ind<-intersect(which(!is.na(b)),which(!is.na(se)))
        n<-length(ind)

        b<-b[ind]
        se<-se[ind]
        prior.se<-prior.se[ind]
        prior.V<-pcm[ind,ind]*matrix(prior.se,nrow=n,ncol=n)*t(matrix(prior.se,nrow=n,ncol=n))
        ccm<-ccm[ind,ind]

        if(n<=1){
            return(NA)
        }

        V<-diag(se^2)
        for(i in 1:(nrow(V)-1)){
            for(j in (i+1):ncol(V)){
                V[i,j]<-sqrt(V[i,i])*sqrt(V[j,j])*ccm[i,j]
                V[j,i]<-V[i,j]
            }
        }
        A<-prior.V+V
        invA<-ginv(A,tol=tolerance)
        quad.form<-t(b) %*% (ginv(V,tol=tolerance) - invA) %*% b
        lbf<-(-0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-as.numeric(determinant(V,logarithm=TRUE)$modulus)))
        lbf<-(lbf + 0.5 * quad.form)
        if(log){
            return(lbf)
        } else if(log10){
            return(lbf/log(10))
        } else {
            return(exp(lbf))
        }
    }

    ##Doing the calculation
    if(class(betas)=="matrix"){
        if(class(prior.sigma)=="matrix"){
            if(length(!which(is.na(prior.sigma)) %in% union(which(is.na(betas)),which(is.na(ses))))>0){
                stop("prior.sigma is not set in at least one instance where both the observed effect and its standard error are recorded.")
            }
            BF<-apply(cbind(betas,ses,prior.sigma),1,get.abf,prior.cor.mat,cryptic.cor.mat,log=log,log10=log10)
        } else if(length(prior.sigma)==1){
            prior.sigma<-matrix(prior.sigma,nrow=nrow(betas),ncol=nstudies)
            BF<-apply(cbind(betas,ses,prior.sigma),1,get.abf,prior.cor.mat,cryptic.cor.mat,log=log,log10=log10)
        } else {
            prior.sigma<-matrix(prior.sigma,nrow=nrow(betas),ncol=nstudies,byrow=TRUE)
            BF<-apply(cbind(betas,ses,prior.sigma),1,get.abf,prior.cor.mat,cryptic.cor.mat,log=log,log10=log10)
        }
    } else {
        if(length(prior.sigma)==1){
            BF<-get.abf(c(betas,ses,rep(prior.sigma,nstudies)),prior.cor.mat,cryptic.cor.mat,log=log,log10=log10)
        } else {
            BF<-get.abf(c(betas,ses,prior.sigma),prior.cor.mat,cryptic.cor.mat,log=log,log10=log10)
        }
    }
    return(BF)
}
