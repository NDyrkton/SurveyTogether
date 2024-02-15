

set.seed(1235)

posrate = 0.5
phi = 0.75
n = 1000
N = 10000000

Y <- numeric(10000)
Y2 <- numeric(10000)


for(i in 1:10000){
  P <- rbinom(n = 1,size = N, prob = posrate)
  Y[i] <- rnoncenhypergeom(n = 1, n1 = P,n2 = N-P, m1 = n, psi = phi)
  Y2[i] <- rbinom(1,n,prob = (1-(1-posrate)^phi))
}
mean(Y)
mean(Y2)

dbinom(Y,n,prob = (1-(1-posrate)^phi))
dnoncenhypergeom(Y, n1 = P,n2 = N-P, m1 = n, psi = phi)

hist(rgamma(10000,shape = 0.9,scale = 2),breaks = 100)


a = -5
b = 3


alpha = 12/(b-a)^3
beta = (b+a)/2

inv.cdf <- function(a,b,x){
  num <- ((a^2)* (  (a^4) - 3*(a^3)*b + 3*(a^2)*(b^2) - a*b^3 + 3*x))^(1/3)
  return((num/a) +b)
}


inv.cdf <- function(a,b,x){
  num <- ((a^2)* (  (a^4) - 3*(a^3)*b + 3*(a^2)*(b^2) - a*b^3 + 3*x))^(1/3)
  return((num/a) +b)
}

U <- runif(10000,min = 0,max = 1)

inv.2 <- function(x,k=1,mu=0,c=1){
  first <- mu + c*(2*x-1)
  return.val <- sign(first)*abs(first)^(1/(2*k+1))
  return(return.val)
}
vals <- inv.2(x = U, k = 0.5, c = 10, mu = 0)

hist(vals,breaks = 100)
