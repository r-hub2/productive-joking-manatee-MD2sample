#' Change Marginals of 2D Data
#' 
#' This function modifies the marginals of some data set
#' 
#' @param dta a list with two 2D data sets
#' @param which method for modifying the data sets.
#' @param theta how much to modify
#' @keywords internal
#' @return a list of functions
#' @export
change.marginals=function(dta, which, theta) {
  if(which=="Exponential") {
    dta$x[,2] =  qexp(dta$x[,2], 1)
    dta$x[,2][dta$x[,2]>3]=3/dta$x[,2][dta$x[,2]>3]
    dta$y[,2] =  qexp(dta$y[,2], theta)
    dta$y[,2][dta$y[,2]>3]=3/dta$y[,2][dta$y[,2]>3]
  }
  if(which=="Linear") {
    if(theta>0)
       dta$y[,2]=theta-1+sqrt((1-theta)^2+4*theta*dta$y[,2])/2/theta
  }
  if(which=="Normal") {
    dta$x[,2] =  qnorm(dta$x[,2])
    dta$y[,2] =  qnorm(dta$y[,2], 0, theta)
  }
  dta
}