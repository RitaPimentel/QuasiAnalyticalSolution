###
### Functions to calculate As, Bs and pstars
###
#
### Given the recursive nature of the solutions, this code should be only used inside of Code3
### Indeed, we are assuming that for each nbar, ps is created, and the needed values in As, Bs and M are already available
#
#
##############################
### Auxiliary calculations ###
##############################
#
### the highest precision is chosen
### if one has a very high maxnbar this level of precision can be needed
options(digits=22)
#
### auxiliary parameters
s1 <- r + lamb
s2 <- r + lamb - mu
s3 <- r + lamb*(1-phi)
#
### roots of the characteristic polynomial (see Section 4 in the paper for details)
d1 <- ((sigma^2/2-mu)+sqrt((sigma^2/2-mu)^2 + 2*sigma^2*s1))/sigma^2
d2 <- ((sigma^2/2-mu)-sqrt((sigma^2/2-mu)^2 + 2*sigma^2*s1))/sigma^2
#
#
#########
### A ###
#########
#
### Auxiliary function to calculate A(nbar, n, k, 0) (see Equation (A.3) in the paper)
### Here we calculate the part that is added for each k (i.e. we do not have the A(nbar, n, k+1, 0))
AuxA0 <- function(nbar,n,k){
  ### part without indicative function
  aux1 <- lamb^(nbar-1-n-k)*((1-d2)*ps[(nbar-k)]/(s2^(nbar-n-k)) + d2*I*phi^(nbar-k-1)*s3/(s1^(nbar-n-k)))
  aux2 <- 0
  aux3 <- 0
  if (k < nbar-n-1){ #k<nbar-n-1
    c0 <- list(c(nbar,n,k,nbar-n-k-1))
    aux21 <- As[c0]*((d2-d1)*(log(ps[(nbar-k)]))-(nbar-n-k-1))*((log(ps[(nbar-k)]))^(nbar-n-k-2))*(ps[(nbar-k)])^d1
    aux22 <- (nbar-n-k-1)*Bs[c0]*((log(ps[(nbar-k)]))^(nbar-n-k-2))*(ps[(nbar-k)]^d2)
    aux2 <- aux2 + aux21 - aux22
    if (k < nbar-n-2){ #k<nbar-2-n
      for (j in 1:(nbar-n-k-2)){
        c1 <- list(c(nbar,n,k,j))
        c2 <- list(c(nbar,n,(k+1),j))
        aux31 <- (As[c1]-As[c2])*((d2-d1)*log(ps[(nbar-k)])-j)*((log(ps[(nbar-k)]))^(j-1))*(ps[(nbar-k)])^d1
        aux32 <-(Bs[c1]-Bs[c2])*j*((log(ps[(nbar-k)]))^(j-1))*(ps[(nbar-k)])^d2
        aux3 <- aux3 + aux31 - aux32
      }
    }
  }
  aux0 <- (ps[(nbar-k)]^(-d1)/(d1-d2))*(aux1 + aux3 + aux2)
  return(aux0)
}
#
### Function to calculate A(nbar, n, 0, 0) (see Equation (A.5) in the paper)
AuxA00 <- function(nbar,n){
  ### part without indicative function
  aux1 <- lamb^(nbar-1-n)*ps[nbar]/(s2^(nbar-n))
  aux2 <- 0
  aux3 <- 0
  if (n < nbar-1){ #n<nbar-1
    c21 <- list(c(nbar,n,1,0))
    c22 <- list(c(nbar,n,0,nbar-1-n))
    aux2 <- d2*Bs[c21]*ps[nbar]^d2-As[c22]*(nbar-1-n+d1*log(ps[nbar]))*((log(ps[nbar]))^(nbar-2-n))*ps[nbar]^d1
    if (n < nbar-2){ #n<nbar-2
      for (l in 1:(nbar-2-n)){
        c31 <- list(c(nbar,n,1,l))
        c32 <- list(c(nbar,n,0,l))
        aux31 <- (As[c31]-As[c32])*(l + d1*log(ps[nbar]))*((log(ps[nbar]))^(l-1))*((ps[nbar])^d1)
        aux32 <- Bs[c31]*(l + d2*log(ps[nbar]))*((log(ps[nbar]))^(l-1))*((ps[nbar])^d2)
        aux3 <- aux3 + aux31 + aux32
      }
    }
  }
  #
  c0 <- list(c(nbar,n,1,0))
  aux0 <- As[c0] + ((ps[nbar]^(-d1))/d1)*(aux3 + aux2 + aux1)
  return(aux0)
}
#
### Function to calculate A(nbar, n, k, j)
AA <- function(nbar, n, k, j){
  if (j!=0){ #j!=0 (see Equation (A.7) in the paper)
    auxA1 <- 0
    for (l in (j-1):(nbar-2-n-k)){
      c1 <- list(c(nbar,n+1,k,l))
      auxA1 <- auxA1 + ((-1)^(l+1-j))*factorial(l)*As[c1]/(factorial(j)*(d1-d2)^(l+2-j))
    }
    auxA1 <- (-2*lamb/sigma^2)*auxA1
    return(auxA1)
  } else{ #j=0 
    if (k==0){ #k=0 and j=0 (see Equation (A.5) in the paper)
      auxA2 <- AuxA00(nbar,n)
      return(auxA2)
    } else{ #k!=0 and j=0 (see Equation (A.3) in the paper)
      ck1 <- list(c(nbar,n,k+1,0))
      auxA3 <- As[ck1] + AuxA0(nbar,n,k)
      return(auxA3)
    }
  }
}
#
#
#########
### B ###
#########
#
### Auxiliary function to calculate B(nbar, n, k, 0) (see Equation (A.4) in the paper)
### Here we calculate the part that is added for each k (i.e. we do not have the B(nbar, n, k+1, 0))
AuxB0 <- function(nbar,n,k){
  ### part without indicative function
  aux1 <- lamb^(nbar-1-n-k)*((d1-1)*ps[(nbar-k)]/(s2^(nbar-n-k)) - d1*I*(phi^(nbar-k-1))*s3/(s1^(nbar-n-k)))
  aux2 <- 0
  aux3 <- 0
  if (k < nbar-n-1){ #k<nbar-n-1
    c0 <- list(c(nbar,n,k,nbar-n-k-1))
    aux21 <- Bs[c0]*((d1-d2)*(log(ps[(nbar-k)]))-(nbar-n-k-1))*((log(ps[(nbar-k)]))^(nbar-n-k-2))*(ps[(nbar-k)])^d2
    aux22 <- (nbar-n-k-1)*As[c0]*((log(ps[(nbar-k)]))^(nbar-n-k-2))*(ps[(nbar-k)])^d1
    aux2 <- aux2 - aux21 + aux22
    if (k < nbar-n-2){ #k<nbar-n-2
      for (j in 1:(nbar-n-k-2)){
        c1 <- list(c(nbar,n,k,j))
        c2 <- list(c(nbar,n,(k+1),j))
        aux31 <- (Bs[c1]-Bs[c2])*((d2-d1)*log(ps[(nbar-k)])+j)*((log(ps[(nbar-k)]))^(j-1))*(ps[(nbar-k)])^d2
        aux32 <- (As[c1]-As[c2])*j*(log(ps[(nbar-k)]))^(j-1)*(ps[(nbar-k)])^d1
        aux3 <- aux3 + aux31 + aux32
      }
    }
  }
  aux0 <- (ps[(nbar-k)]^(-d2)/(d1-d2))*(aux1 + aux3 + aux2)
  return(aux0)
}
#
### Function to calculate B(nbar, n, k, j)
BB <- function(nbar, n, k, j){
  if (j!=0){ #j!=0 (see Equation (A.8) in the paper)
    auxB1 <- 0
    for (l in (j-1):(nbar-2-n-k)){
      c1 <- list(c(nbar,n+1,k,l))
      auxB1 <- auxB1 + ((-1)^(l+1-j))*factorial(l)*Bs[c1]/(factorial(j)*(d2-d1)^(l+2-j))
    }
    auxB1 <- (-2*lamb/sigma^2)*auxB1
    return(auxB1)
  } else{ #j=0 (see Equation (A.4) in the paper)
    ck1 <- list(c(nbar,n,k+1,0))
    auxB2 <- Bs[ck1] + AuxB0(nbar,n,k)
    return(auxB2)
  }
}
#
#
#############
### pstar ###
#############
#
### Auxiliary function to find pstar(nbar,n) 
### Here pstar(nbar,n) is represented by the variable p
### We calculate B(nbar, n, 1, 0) as a function of pstar(nbar,n) (see Equation (A.4) in the paper)
### For n=nbar-1, pstar(nbar,nbar-1) has an explicit formula
### This function is only called when n<nbar-1
BBp <-  function(nbar, n, p){
  ### This function is B summed for k from 1 to nbar-1-n 
  ### Note that the last term in the sum (i.e. when k=nbar-1-n) has the unknown pstar(nbar,n), represented here by p
  auxBp1 <- 0
  if (n < nbar-2){
    for (l in 1:(nbar-2-n)){
      auxBp1 <- auxBp1 + AuxB0(nbar,n,l)
    }
  }
  auxBp1 <- auxBp1 + ((d1-1)*p/s2 - d1*I*(phi^n)*s3/s1)*p^(-d2)/(d1-d2)
  return(auxBp1)
}
#
### Auxiliary function to find pstar(nbar,n) (see Equation (A.9) in the paper)
### This represents the right hand side of the equation as a function of p
### For n=nbar-1, pstar(nbar,nbar-1) has an explicit formula
### This function is only called when n<nbar-1
Topstar <- function(nbar,n, p){
  ### part without indicative function
  auxt <- lamb^(nbar-n-1)*((d1-1)*ps[nbar]/(s2^(nbar-n)) - d1*I*(phi^(nbar-1))*s3/(s1^(nbar-n)))
  ### the function is only called when n<nbar-1 
  c1 <- list(c(nbar,n,0,nbar-1-n))
  ### note that BBp depends on p
  auxt1 <- (d1-d2)*BBp(nbar,n,p)*ps[nbar]^d2 + As[c1]*(nbar-1-n)*((log(ps[nbar]))^(nbar-2-n))*(ps[nbar])^d1
  auxt <- auxt + auxt1
  if (n < nbar-2){ #n<nbar-2
    for (j in 1:(nbar-2-n)){
      c2 <- list(c(nbar,n,0,j))
      c3 <- list(c(nbar,n,1,j))
      auxt2 <- (As[c2]-As[c3])*j*((log(ps[nbar]))^(j-1))*(ps[nbar])^d1
      auxt3 <- Bs[c3]*((d1-d2)*((log(ps[nbar]))^j)-j*((log(ps[nbar]))^(j-1)))*(ps[nbar])^d2
      auxt <- auxt2 + auxt3 + auxt
    } 
  }
  return(auxt)  
}
#
### Function to find pstar(nbar,n) (see Equations (A.9) and (A.10) in the paper)
pstar <- function(nbar, n){
  ### The initial code (including two if cycles) is commented because in fact this function is only used for n<=nbar-2
  ### Indeed, for n=nbar and n=nbar-1, the pstars are calculated within Code3
  # if (n==nbar){ #when we are in the last level (n=nbar), we assume that pstar is zero
  #   return(0)
  # } else{
  #   if (n==nbar-1){ #For the level n=nbar-1, there is an explicit expression for pstar (see Equation (A.10) in the paper)
  #     auxp1 <- I*phi^(nbar-1)*s3*(d2-1)/d2
  #     return(auxp1)
  #   } else{ #n<nbar-1
  #   ### given that p(nbar,n) increases with nbar, then p(nbar-1,n) is a minimum boundary for p(nbar,n)
      xmin <- M[(n+1),(nbar-1)]
      ### the difference between consecutive pstars decreases with nbar, then p(nbar,n)-p(nbar-1,n) is lower than p(nbar-1,n)-p(nbar-2,n)  
      ### then a maximum boundary for p(nbar,n) is p(nbar-1,n) + (p(nbar-1,n) - p(nbar-2,n))
      ### however, p(0,0)=0 is not defined in the matrix M, and it is needed to calculate p(2,0), then we need to set it here equal to 0
      if ((nbar==2) & (n==0)){
        prevpstar <- 0
      } else{
        prevpstar <- M[(n+1),(nbar-2)]
      }
      ### we add a small increment to xmax (10^-15) because for large values of nbar the pstars can be already equal (with the maximum precision), thus xmin=xmax i.e. the uniroot routine will produce an error message
      xmax <- 2*xmin - prevpstar + 10^-15
      fzero <- uniroot(Topstar, c(xmin,xmax), tol=1e-30, maxiter=1000, extendInt="yes", nbar=nbar, n=n)
      auxp2 <- fzero$root
      return(auxp2)
  #   }
  # }
}
#
#
### END