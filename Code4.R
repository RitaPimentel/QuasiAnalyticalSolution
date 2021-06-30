###
### To create the values to plot the value functions
###
#
##################################
### Creating auxiliary objects ###
##################################
#
### All possible values of nbar
vnbar <- 2:maxnbar
#
### Values in the X-axis of the plots
xx <- seq(minvalue,maxvalue,funcinc)
lenxx <- length(xx)
rr <- 1:lenxx
#
#
### To plot the value functions, we create a matrix to save the values
### Each column represents one of the values of XX, which the image will be saved on Fxx
### Each row represents one of the functions:
### First maxnbar-1 functions are V(nbar,0) for nbar=2,...,nbar=maxnbar
### The next maxnbar-1 functions are V(nbar,1) for nbar=2, nbar=3,...,nbar=maxnbar
### The next maxnbar-2 functions are V(nbar,2) for nbar=3, nbar=4,...,nbar=maxnbar
### ...
### Until the last function, which is V(mavnbar,nbar-1) for nbar=maxnbar.
# 
### Total number of functions 
nfunctions <- sum(1:maxnbar)-1
#
### Matrix to save the values of all functions
Fxx <- matrix(0,nfunctions,lenxx)
#
### Auxiliary function to locate where the functions for each n start and end
### For example, for V(nbar,i) the functions are available from the row posfunc[j] until posfunc[j+1]-1 with j=i+1
posfunc <- c(0,cumsum(rev(1:maxnbar)))
#
#
######################
### Value function ###
######################
#
### Function in the stopping region for each level n (see Equation (A.1) in the paper)
g <- function(n,p){
  gresult <- p/(r-mu) - I*phi^n
  return(gresult)
}
#
### Part of the function in the continuation region, for each level n and k (see Equation (A.2) in the paper)
### Here we only have the second part which does not involve the sum with As and Bs
f <- function(nbar,n,k,p){
  fresult <- (lamb^(nbar-n-k))*(p/((r-mu)*(s2^(nbar-n-k))) - I*phi^(nbar-k)/(s1^(nbar-n-k)))
  return(fresult)
}
#
#
### Fill in the functions with the values in the stopping region
### Now we are calculation the value in the stopping for all values of xx
### Afterwards we will substitute the values in the continuation region by the correct ones
for (xi in 1:lenxx){
  for (vi in 1:maxnbar){
    ### for vi equal to 1 this will start in 0, which is not an acceptable index, but it is ignored and it only starts in 1, as we intend
    Fxx[posfunc[vi]:(posfunc[vi+1]-1),xi]<-g(vi-1,xx[xi])
  }
}
#
### Calculate the values of the functions in the continuation region
for (nbar in vnbar){
  auxps <- c(0,rev(M[(1:nbar),nbar]))
  for (i in 2:length(auxps)){
    k <- i-2
    for (xi in rr[(auxps[i-1] <= xx) & (xx < auxps[i])]){
      for (tt in 1:(nbar-k)){
        n <- tt-1
        aux <- 0
        for (j in 0:(nbar-1-n-k)){
          cc1 <- list(c(nbar,n,k,j))
          aux <- aux + As[cc1]*((log(xx[xi]))^j)*(xx[xi]^d1) + Bs[cc1]*((log(xx[xi]))^j)*(xx[xi]^d2)
        }
        Fxx[(posfunc[tt]+nbar-tt),xi] <- aux + f(nbar,n,k,xx[xi])
      }
    }
  }
}
#
#
### For plotting 
### we need to consider the first element as 1 and not 0, because it is the index where the values of the first function are saved 
posfunc[1] <- 1
#
#
# END