###
### Procedure to create a matrix with all pstars
###
#
##################################
### Creating auxiliary objects ###
##################################
#
### All possible values of nbar
vnbar <- 2:maxnbar
#
### Matrix to save the pstar values for each nbar
### row 1: pstar0, row2: pstar1,...; column1: nbar=1, column2: nbar=2,...
### Thus, row i, column j has pstar with nbar=j and n=i-1
M <- matrix(0,maxnbar,maxnbar)
#
### Initializing the matrix for nbar=1, i.e. pstar(1,0) (see Equation (A.10) in the paper)
### In fact, this is not an "admissible" pstar because nbar=1 (and nbar >= 2)
### However, it is useful to calculate pstar(2,0)
M[1,1] <- I*s3*(d2-1)/d2
#
### create sparse tensors to save As and Bs values
dims <- c(maxnbar, maxnbar,maxnbar,maxnbar)
subs <- list(c(0,0,0,0))
vals <- 0
As <- sptensor(subs, vals, dims)
Bs <- sptensor(subs, vals, dims)
#
#
##############################
### Calculating the pstars ###
##############################
#
### For each value of nbar, calculate the pstars for n=0,1,...,nbar-1
for (nbar in vnbar){
  #
  ### vector of pstars for this nbar
  ps <- rep(0,nbar)
  #
  ### pstar(nbar,nbar-1) has a closed form, then we use it directly here (see Equation (A.10) in the paper)
  ps[nbar] <- I*s3*phi^(nbar-1)*(d2-1)/d2
  #
  ### A(nbar,nbar-1,0,0) has a closed form, then we use it directly here (see Equation (A.6) in the paper)
  cc0 <- list(c(nbar,nbar-1,0,0))
  As[cc0] <- ps[nbar]^(1-d1)/(d1*s2)
  #
  for (n in (nbar-2):0){
    ### fill all As and Bs with j!=0 (for all ks)
    for (k in 0:(nbar-2-n)){
      for (j in 1:(nbar-1-n-k)){
        cc1 <- list(c(nbar,n,k,j))
        As[cc1] <- AA(nbar,n,k,j)
        if (k==0){
          Bs[cc1] <- 0
        } else{
          Bs[cc1] <- BB(nbar,n,k,j)
        }
      }
    }
    ### calculate the pstar for this level, n
    pstar_now <- pstar(nbar,n)
    ps[(n+1)] <- pstar_now
    ### fill all As and Bs with j=0 (for all ks)
    for (k in (nbar-1-n):0){
      cc2 <- list(c(nbar,n,k,0))
      As[cc2] <- AA(nbar,n,k,0)
      if (k==0){
        Bs[cc2] <- 0
      } else{
        Bs[cc2] <- BB(nbar,n,k,0)
      }
    }
    # in case one wants to follow the status, the nbar, n and time can be printed
    # print(c(nbar,n))
    # print(Sys.time())
    rm(pstar_now)
  }
  ### fill in the matrix M with all pstars for this level nbar
  M[(1:nbar),nbar] <- ps
}
#
#
#########################################
### Saving the matrix with all pstars ###
#########################################
#
### Matrix with all the pstars
### row 1: pstar0, row2: pstar1,...; column1: nbar=2, column2: nbar=3,...
### Thus, row i, column j has pstar with nbar=j+1 and n=i-1
Mpstar <- M[,-1]
row.names(Mpstar) <- 0:(maxnbar-1) #n
colnames(Mpstar) <- vnbar          #nbar
#
#
### END