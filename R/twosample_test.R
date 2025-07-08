#' Multivariate Two-Sample Tests
#' 
#' This function runs a number of two sample tests using Rcpp and parallel computing.
#' 
#' @param  x  a matrix of numbers if data is continuous or a vector of counts  if data is discrete, or a list of x and y
#' @param  y a matrix of numbers if data is continuous or a vector of counts  if data is discrete.
#' @param  vals_x =NA, a vector of values for discrete random variable
#' @param  vals_y =NA, a vector of values for discrete random variable
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  B =5000, number of simulation runs for permutation test
#' @param  nbins =c(5,5), for chi square tests (2D only).
#' @param  minexpcount =5, lowest required count for chi-square test
#' @param  Ranges =matrix(c(-Inf, Inf, -Inf, Inf),2,2) a 2x2 matrix with lower and upper bounds
#' @param  DoTransform =TRUE, should data be transformed to interval (0,1)?
#' @param  samplingmethod ="Binomial" for Binomial sampling or "independence" for independence sampling
#' @param  rnull function to generate new data sets for simulation
#' @param  SuppressMessages =FALSE, should messages be printed
#' @param  LargeSampleOnly =FALSE should only methods with large sample theories be run?
#' @param  maxProcessor =1 maximum number of cores to use. If missing (the default) no parallel processing is used.
#' @param  doMethods ="all" Which methods should be included?
#' @return A list of two numeric vectors, the test statistics and the p values. 
#' @examples
#' #Two continuous data sets from a multivariate normal:
#' x = mvtnorm::rmvnorm(100, c(0,0))
#' y = mvtnorm::rmvnorm(120, c(0,0))
#' twosample_test(x, y, B=500)
#' #Using a new test, this one is an (included) chi square test. 
#' #Also enter data as a list:
#' TSextra=list(which="statistics", nbins=rbind(c(3,3), c(4,4)))
#' dta=list(x=x, y=y)
#' twosample_test(dta, TS=chiTS.cont, TSextra=TSextra, B=500)
#' #Two discrete data sets from some distribution:
#' x = table(sample(1:4, size=1000, replace = TRUE))
#' y = table(sample(1:4, size=1000, replace = TRUE, prob=c(1,2,1,1)))
#' vals_x=rep(1:2,2)
#' vals_y=rep(1:2, each=2)
#' twosample_test(x, y, vals_x, vals_y, B=500)
#' #Run a discrete chi square test and enter the data as a matrix:
#' TSextra=list(which="statistics")
#' dta=cbind(x=x, y=y, vals_x=vals_x, vals_y=vals_y)
#' twosample_test(dta, TS=chiTS.disc, TSextra=TSextra, B=500)
#' @export 
twosample_test=function(x, y, vals_x=NA, vals_y=NA, TS, TSextra, B=5000, 
          nbins=c(5,5),minexpcount =5, Ranges=matrix(c(-Inf, Inf, -Inf, Inf),2,2),
          DoTransform=TRUE, samplingmethod="Binomial", rnull,
          SuppressMessages=FALSE, LargeSampleOnly=FALSE, 
          maxProcessor=1, doMethods="all") {
    if(!is.numeric(samplingmethod))
      samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(missing(y)) {
       if(is.list(x)) {#Continuous Data
         if(!SuppressMessages) message("Data is assumed to be continuous")
         Continuous=TRUE
         dta=x
         y=x$y
         x=x$x 
         Dim=ncol(x)
       }
       else {#Discrete data
         if(!SuppressMessages) message("Data is assumed to be discrete")
         Continuous=FALSE
         dta=x
         x=matrix(1:4,2,2)#just some dummy numbers
         y=matrix(1:4,2,2)
         dta=list(x=dta[,"x"], y=dta[,"y"], vals_x=dta[,"vals_x"], vals_y=dta[,"vals_y"])
         DoTransform=FALSE
       }
    }
    else {
       Continuous=ifelse(any(is.na(c(vals_x, vals_y))), TRUE, FALSE)
       if(Continuous) {
          if(!SuppressMessages) message("Data is assumed to be continuous")
          dta=list(x=x, y=y)
          Dim=ncol(x)
       } 
       else {
          if(!SuppressMessages) message("Data is assumed to be discrete")
          dta=list(x=x, y=y, vals_x=vals_x, vals_y=vals_y)
          x=matrix(1:4,2,2)#just some dummy numbers
          y=matrix(1:4,2,2)
          DoTransform=FALSE
       }    
    }
    if(test_methods(doMethods, Continuous)) return(NULL)
    if(Continuous) {
      if(nrow(y)<nrow(x)) { #switch x and y so that nx>ny
        tmp=y
        y=x
        x=tmp
        dta=list(x=x, y=y)
      }
      rawdta=dta
      if(missing(DoTransform)) {
        z=c(x, y)
        DoTransform=FALSE
        if(min(z)<0 | max(z)>1) DoTransform=TRUE
      }  
      if(DoTransform) {
          dta=transform01(dta)
          x=dta$x
          y=dta$y
          Ranges=matrix(c(0, 1, 0, 1),2,2)
      }
    }
    if(missing(TSextra)) TSextra=list(aaa=0)
    if(Continuous) 
        TSextra = c(TSextra, 
                knn=function(x) FNN::get.knn(x, 5)$nn.index,
                dist=function(dta) find_dist(dta),
                distances=list(find_dist(dta)),
                DoTransform=DoTransform)
    else TSextra = c(TSextra, 
                     dist=function(dta) NULL,
                     organize=function(dta) dta=dta[order(dta[,1], dta[,2]), ],
                     samplingmethod=samplingmethod)
                     
    if(!missing(rnull)) {
       if(Continuous) 
           TSextra=c(TSextra, rnull=rnull, rawdta=list(rawdta))
       else TSextra=c(TSextra, rnull=rnull)
    }   
    outchi=list(statistics=NULL, p.value=NULL)
    outpvals=list(statistics=NULL, p.value=NULL)
# what methods are to be run?  
    CustomTS=TRUE
    if(missing(TS)) { # do included methods
      CustomTS=FALSE
      if(Continuous) {
          if(Dim==2) {
            if(length(nbins)==1) nbins=c(nbins, nbins)
            if(missing(rnull))
               outchi = chisq2D_test_cont(x,y, Ranges, nbins, minexpcount)
          }  
          if(missing(rnull))
             outpvals=TS_cont_pval(x, y) #Methods that find p values
          typeTS = 1
          TS=TS_cont
          dta=list(x=x, y=y)
      }
      else {
          typeTS=4
          if(missing(rnull))
             outchi = chisq2D_test_disc(dta, minexpcount)
          TS=TS_disc
      }    
    }  
    else { # do user-supplied tests
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
         if(maxProcessor>1) {
            if(!SuppressMessages) message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
            maxProcessor=1
         }  
      }
      if(Continuous) typeTS=length(formals(TS))
      else typeTS=ifelse(length(formals(TS))==5, 6, 5)  
    }
    TS_data=calcTS(dta, TS, typeTS, TSextra)
    if(any(is.null(names(TS_data)))) {
      if(!SuppressMessages) message("output of TS routine has to be a named vector!")
      return(NULL)
    }
#  set number of processors for parallel programming
#  and do a timing check    
    if(missing(maxProcessor))
       maxProcessor=max(parallel::detectCores(logical = FALSE)-1, 1)  
    if(maxProcessor>1) {
      tm=timecheck(dta, TS, typeTS, TSextra)
      if(2*tm[1]*B<20 | B<2*maxProcessor) {
        maxProcessor=1
        if(!SuppressMessages) message("maxProcessor set to 1 for faster computation")
      }
      else if(!SuppressMessages) message(paste("Using ",maxProcessor," cores.."))  
    } 
# if either only one core is present, B=0 or maxProcessor=1, run testC. 
    if(!LargeSampleOnly) {
      if(maxProcessor==1) outTS = testC(dta, TS, typeTS, TSextra, B=B)
      else {
# run testC in parallel. Use one less core than is present, at most maxProcessor.    
        cl=parallel::makeCluster(maxProcessor)
        z=parallel::clusterCall(cl, testC, dta=dta, TS=TS, 
                              typeTS=typeTS, TSextra=TSextra, 
                              B=round(B/maxProcessor))
        parallel::stopCluster(cl)
# average p values over cores    
        p=z[[1]]$p.values
        for(i in 2:maxProcessor) p=p+z[[i]]$p.values
        p = round(p/maxProcessor, 4)  
        outTS = list(statistics=z[[1]]$statistics, p.values=p)
      }
    }
    else outTS=list(statistics=NULL, p.values=NULL)
    if(CustomTS) return(signif.digits(outTS))
    s = c(outTS$statistics, outpvals$statistics, outchi$statistic)
    p = c(outTS$p.values, outpvals$p.values, outchi$p.value)
    if(doMethods[1]!="all") {
      s=s[doMethods]
      p=p[doMethods]
    }  
    signif.digits(list(statistics=s, p.values=p))
}
