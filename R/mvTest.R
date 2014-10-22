####################################### 
# shapeQTL mapping experiment with R
#
# Nicolas Navarro - 2013-2014
########################################  
# Utilities function for multivariate linear testing
Pillai.test <- function (SSCPef,SSCPer,dfef,dfer,p=qr(SSCPef+SSCPer)$rank) 
{ 
  q <- dfef
  v <- dfer
  s <- min(q,p)
  m <- 0.5*(abs(p-q)-1)
  n <- 0.5*(v-p-1)
  df1 <- s*(2*m+s+1)
  df2 <- s*(2*n+s+1)
  #!! HE^-1 not symmetric !!
  A <- eigen(SSCPef%*%ginv(SSCPer),only.values=TRUE)$values
  A <- Re(A)
  V <- sum(A/(1+A))
  Fapprox <- V/(s-V) * (df2/df1)
  return(-pf(Fapprox,df1,df2,lower.tail=F,log.p=T)/log(10))
}
LikelihoodRatio.test <- function(SSCPfull,SSCPred,dfef,dfer,p=qr(SSCPfull)$rank)
{
  #Following LOD definition of Haley & Knott (1992) in HK regression
  scale <- -(dfer - 0.5*(p - dfef + 1))
  eig <- eigen(SSCPfull,only.values=TRUE)$values
  det.full <- prod(eig[eig>.Machine$double.eps])
  eig <- eigen(SSCPred,only.values=TRUE)$values
  det.red <- prod(eig[eig>.Machine$double.eps])
  return(scale*log10(det.full/det.red))
}
Hotelling.test <- function(SSef,SSer,dfef,dfer,p=qr(SSef+SSer)$rank){
  #D'apres Claude 2008, p.252 
  k <- dfef
  w <- dfer
  s <- min(k,p)
  m <- (w-p-1)/2
  t1 <- (abs(p-k)-1)/2
  Ht <- sum(diag(SSef%*%ginv(SSer)))
  Fapprox <- Ht*(2*(s*m+1))/(s^2*(2*t1+s+1))
  ddfnum <- s*(2*t1+s+1)
  ddfden <- 2*(s*m+1)
  -pf(Fapprox,ddfnum,ddfden,lower.tail=F,log.p=T)/log(10)
}

