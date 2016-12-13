wakefield.abf<-function(beta,se,prior.sigma,log=FALSE,log10=FALSE){
    #beta is a single value or numeric vector. It should represent the observed effect size of a SNP on a trait from a genome-wide association study.
    #se is a single value or numeric vector of the same length as beta. It should represent the standard error of the effect size of a SNP on a trait from a genome-wide association study.
    #prior.sigma is a single value or numeric vector. If it is a single value, it does not need to be the same length as se and beta. This is the prior on true effect sizes.
    #log sets whether the Bayes factor is returned in log (natural log) space. Cannot be set true if log10 is also true.
    #log10 sets whether the Bayes factor is return in log10 space. Cannot be set true if log is also true.


    #Error checking
    if(class(beta)!="numeric"){
        stop("beta values are not numeric.")
    }
    if(class(se)!="numeric"){
        stop("se values are not numeric.")
    }
    if(class(prior.sigma)!="numeric"){
        stop("prior.sigma values are not numeric.")
    }
    if(all(is.na(beta))){
        stop("all values of beta are NA.")
    }
    if(all(is.na(se))){
        stop("all values of se are NA.")
    }
    if(any(se<0,na.rm=TRUE)){
        stop("cannot accept negative values of se.")
    }
    if(any(se==0,na.rm=TRUE)){
        stop("one of the values of se is 0. This yields an approximate Bayes factor of 0.")
    }
    if(length(beta)!=length(se)){
        stop("beta and se do not have the same length.")
    }
    if(!length(prior.sigma) %in% c(1,length(beta))){
        stop("prior.sigma must be either a single value or a numeric vector of equal length to beta and se.")
    }
    if(log && log10){
        stop("both log and log10 have been set as true. Only one of these may be be true at a time.")
    }
    #Setting up variables for Wakefield's calculation
    v<-se^2
    w<-prior.sigma^2
    zsq<-(beta^2)/v
    if(log){
        #Performing the calculation in natural log space
        log_wabf<-((zsq/2)*(w/(v+w)))-log(sqrt((v+w)/v))
        return(log_wabf)
    } else if(log10){
        #Performing the calculation in natural log space
        log_wabf<-((zsq/2)*(w/(v+w)))-log(sqrt((v+w)/v))
        #Converting to log10 space
        return(log_wabf/log(10))
    } else {
        #Calculating Wakefield's approximate Bayes factor
        wabf<-exp((zsq/2)*(w/(v+w)))/sqrt((v+w)/v)
        return(wabf)
    }
}
