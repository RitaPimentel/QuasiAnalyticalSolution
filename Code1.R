#
### Implementation of the approximated solution to the general case, presented in Appendix A.
### The plots and tables in the paper are produced using this code.
#
### clear the environment
rm(list = ls())
#
### set the directory where the files are located
# setwd("...")
#
###################
### 1. Packages ###
###################
# 
### in case the packages are not installed yet
# install.packages("tensorr")
# install.packages("matrixStats")
# install.packages("xlsx")
# Note that for the xlsx package, one needs to have downloaded and installed the 64 bit version of JAVA
#
library(tensorr)
library(xlsx)
library(matrixStats)
#
#
##########################################
### 2. Initial parameters of the model ###
##########################################
#
### Note: to be chosen by the user (see Section 2 in the paper for details)
#
r <- 0.05       # interest rate (should be >0)
mu <- 0.03      # drift (should be such that r>mu)
sigma <- 0.1    # volatility (should be >0)
lamb <- 0.1     # Poisson process' rate (should be >0)
I <- 1          # initial investment cost (should be >0)
phi <- 0.9      # parameter leading the decreasing in the investment cost (should be between 0 and 1)
#
#
### If one only intends to plot the graphics regarding the comparative statics go to step 7
#
#
#########################################################
### 3. General parameters for the plots and/or tables ###
#########################################################
#
# Note: to be chosen by the user, depending on what is intended
#
### Maximum level considered for the truncated problem (see Section 3 in the paper)
### For Figures 3 and 4 and Tables C1 and C2 we have maxnbar <- 10
### For Figure 5 we have maxnbar <- 30
maxnbar <- 10   # should be larger than 2
#
### Level n for which we want to plot the pstars
### For Figures 3, 4 and 5 we have lev <- 0
lev <- 0       # should be between 0 and maxnbar-1
#
#
#####################################################
### 4. Upload functions and create needed objects ###
### to calculate pstars and plot value functions  ###
#####################################################
#
#
### Functions to calculate As, Bs and pstars
source("Code2.R")
#
### Create a matrix with the pstar values
source("Code3.R")
#
#
### Save the objects created in Code3 that will be used for Steps 5 and 6
### Note that these objects are valid only for the parameters chosen before (if the parameters are changed, then one needs to run Step 4 again)
#
save(M,file="M.RData")
save(As,file="As.RData")
save(Bs,file="Bs.RData")
#
#
################################################
### 5. For plots and tables regarding pstars ###
################################################
#
### define the precision
options(digits=22)
#
### To plot Figure 3 (a) in the paper
### Plot values p(nbar,n) for different values of nbar and a fixed level n
#
### auxiliary vector to indicate the indexes of the pstars that should be plotted
adnbar <- 1:(maxnbar-1)
if (lev>1){
  adnbar <- lev:(maxnbar-1)
}
### vector with the values of pstars
yy <- Mpstar[(lev+1),adnbar]
#
### it will save the plot in the same directory with the name Fig3a_pstar_maxnbar
fn1 <- paste("Fig3a_pstar_",maxnbar, ".jpg", sep = "")
jpeg(fn1, width = 714, height = 485)
plot((adnbar+1),yy, xlab=" ", ylab=" ")
dev.off()
#
#
### To plot Figure 3 (b) in the paper
### plot of the different between two consecutive values of pstar, i.e. pstar(nbar+1,n)-pstar(nbar,n)
dyy <- yy[2:length(yy)]-yy[1:(length(yy)-1)]
#
### it will save the plot in the same directory with the name Fig3b_diffpstar_maxnbar
fn2 <- paste("Fig3b_diffpstar_",maxnbar, ".jpg", sep = "")
jpeg(fn2, width = 714, height = 485)
plot((adnbar[-1]+1),dyy, xlab=" ", ylab=" ")
dev.off()
#
#
### To produce Table C1 in the paper
### represent p(nbar,n) for different possible combinations of nbar and n
### it will save an Excel file in the same directory with the name TableC1_pstar_maxnbar
dd1 <- paste("TableC1_pstar_",maxnbar, ".xlsx", sep = "")
write.xlsx(Mpstar, dd1, col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)
#
#
### To produce Table C2 in the paper
### p(nbar,n) - p(nbar-1,n) for n=lev and all admissible values of nbar
### it will save an Excel file in the same directory with the name TableC2_pstar_maxnbar
Dpstar <- rowDiffs(Mpstar[-nrow(Mpstar),])
for (i in 1:(ncol(Dpstar)-1)){
  Dpstar[i+2,i] <- 0 
}
row.names(Dpstar) <- head(row.names(Mpstar),-1)   #n
colnames(Dpstar) <- colnames(Mpstar)[-1]          #nbar
#
dd2 <- paste("TableC2_diffpstar_",maxnbar, ".xlsx", sep = "")
write.xlsx(Dpstar, dd2, col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE)
#
#
#################################################
### 6. For plots regarding the value function ###
#################################################
#
### Note: XX values to be chosen by the user, depending on what is intended
#
### values of XX that we want to calculate the images (min and max) to plot afterwards (and also the step between points)
minvalue <- 0.0001  # should be positive (and most likely close enough to zero)
maxvalue <- 0.070   # should be larger than the pstar(nbar,n) (in order to plot the function in the continuation and stopping regions)
funcinc <- 0.0005   # increment
#
#
### Create a matrix with the value functions
source("Code4.R")
#
#
### Save the matrix created in Code4 that will be used in this step
### Note that this object is only valid for the parameters chosen before (if the parameters are changed, then one needs to run Step 4 again)
#
save(Fxx,file="Fxx.RData")
#
#
### indexes in object Fxx of all the functions V(nbar,n) for level n=lev
ro <- (posfunc[lev+1]):(posfunc[lev+2]-1)
#
#
### To plot Figure 4 (a) in the paper
### plot the value functions V(nbar,n) for n=lev and all admissable values of nbar
### it will save the plot in the same directory with the name Fig4a_AllFunc_maxnbar_lev
vf1 <- paste("Fig4a_AllFunc_",maxnbar, "_", lev, ".jpg", sep = "")
jpeg(vf1, width = 714, height = 485)
plot(xx, Fxx[ro[1],], type='l', xlab=" ", ylab=" ")
for (wi in 2:length(ro)){
  lines(xx,Fxx[ro[wi],],lty=wi)
}
dev.off()
#
#
### To plot Figure 4 (b) in the paper
### plot the differences between two consecutive value functions
### V(nbar,n) - V(nbar-1,n) for n=lev and all admissible values of nbar
### it will save the plot in the same directory with the name Fig4b_AlldiffFunc_maxnbar_lev
vf2 <- paste("Fig4b_AlldiffFunc_",maxnbar, "_", lev, ".jpg", sep = "")
jpeg(vf2, width = 714, height = 485)
plot(xx, Fxx[ro[2],]-Fxx[ro[1],], type='l', xlab=" ", ylab=" ")
for (wi in 2:(length(ro)-1)){
  lines(xx,Fxx[ro[wi+1],]-Fxx[ro[wi],],lty=wi)
}
dev.off()
#
#
### To plot Figure 5 in the paper
### plot the value function V(maxnbar,n) with n=lev
### it will save the plot in the same directory with the name Fig5_ValueFunc_maxnbar_lev
vf3 <- paste("Fig5_ValueFunc_",maxnbar, "_", lev, ".jpg", sep = "")
jpeg(vf3, width = 714, height = 485)
plot(xx, Fxx[posfunc[lev+2]-1,], type='l', xlab=" ", ylab=" ")
dev.off()
#
#
###################################################################################
### 7. For plots regarding the comparative statics (see Section 5 in the paper) ###
###################################################################################
#
### Functions to be used in the calculations
source("Code5.R")
#
### Note: to be chosen by the user, in order to see the differences depending on phi
### For Figure 6 (a) is phi1 <- 0.90 and for Figure 6 (b) is phi1 <- 0.99
phi1 <- 0.90
phi2 <- 0.99
#
### Discretization of lambda that we consider in the plots
ll <- seq(0.0001,1,0.0001)
#
### To plot Figure 6 (a) in the paper
### it will save the plot in the same directory with the name Fig6a_phi1
fc1 <- paste("Fig6a_",phi1, ".jpg", sep = "")
jpeg(fc1, width = 714, height = 485)
plot(ll,p21(ll,phi1), xlab=" ", ylab=" ", type="l")
dev.off()
#
### To plot Figure 6 (b) in the paper
### it will save the plot in the same directory with the name Fig6b_phi2
fc2 <- paste("Fig6b_",phi2, ".jpg", sep = "")
jpeg(fc2, width = 714, height = 485)
plot(ll,p21(ll,phi2), xlab=" ", ylab=" ", type="l")
dev.off()
#
#
### END