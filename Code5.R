###
### Comparative statics 
###
### We present the plots for Section 5 in the paper
### Comparative statics with respect to lambda
### The plots represent pstar(2,1), i.e. nbar=2 and n=1, for values of lambda between 0 and 1, and for two possible values of phi (phi1 and phi2) 
#
#
### roots of the characteristic polynomial as a function of lambda
Fd1 <- function(lamb,phi){
  dd1 <- ((sigma^2/2-mu)+sqrt((sigma^2/2-mu)^2 + 2*sigma^2*(r + lamb)))/sigma^2
  return(dd1)
}
#
Fd2 <- function(lamb,phi){
  dd2 <- ((sigma^2/2-mu)-sqrt((sigma^2/2-mu)^2 + 2*sigma^2*(r + lamb)))/sigma^2
  return(dd2)
}
#
### pstar(2,1) with respect to lambda and phi (see Equation (A.10) or Equation (11) in the paper)
p21 <- function(lamb,phi){
  p <- I*phi*(r + lamb*(1-phi))*(Fd2(lamb)-1)/Fd2(lamb)
  return(p)
}
#
#
### END