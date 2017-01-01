#part (b)
set.seed(123)
gibbs <- function(n.sims, y, burnin, thin)
{
  n<-length(y)
  x.draws <- matrix(NA, nrow=(n.sims-burnin)/thin,ncol=10)
  mu.draws <- c()   # initialize vector that will store draws from the full conditional 
  sig2.draws<- c() # initialize vector that will store draws from theh full conditional 
  sig2.cur = 2
  mu.cur = 5
  xi.update <- function(yi,mu,sig2) { # updates xi using the full conditional distribution
    xi = rtruncnorm(1, a=yi-.5, b=yi+.5, mean = mu, sd = sqrt(sig2))
    return(xi)
  }
  mu.update <- function(y,x,mu,sig2) { # updates using the MH-RW
    mu.mh.rw(y,x,mu,sig2)
  }
  sig2.update <- function(y,x,mu,sig2){
    sig2.mh.rw(y,x,mu,sig2)
  }
  for (i in 1:n.sims) {  # simulates and calls update functions to simulate parameters
    x1.cur <- xi.update(y[1],mu.cur,sig2.cur)
    x2.cur <- xi.update(y[2],mu.cur,sig2.cur)
    x3.cur <- xi.update(y[3],mu.cur,sig2.cur)
    x4.cur <- xi.update(y[4],mu.cur,sig2.cur)
    x5.cur <- xi.update(y[5],mu.cur,sig2.cur)
    x6.cur <- xi.update(y[6],mu.cur,sig2.cur)
    x7.cur <- xi.update(y[7],mu.cur,sig2.cur)
    x8.cur <- xi.update(y[8],mu.cur,sig2.cur)
    x9.cur <- xi.update(y[9],mu.cur,sig2.cur)
    x10.cur <- xi.update(y[10],mu.cur,sig2.cur)
    x.cur = c(x1.cur, x2.cur, x3.cur, x4.cur, x5.cur, x6.cur, x7.cur, x8.cur, x9.cur, x10.cur)
    mu.cur <- mu.update(y,x.cur,mu.cur,sig2.cur)
    sig2.cur <- sig2.update(y,x.cur,mu.cur,sig2.cur)
    if (i > burnin & (i - burnin)%%thin == 0) {  # applys burn-in and thining to the simulated data
      x.draws[(i - burnin)/thin,] <- x.cur
      mu.draws[(i - burnin)/thin] <- mu.cur
      sig2.draws[(i - burnin)/thin] <- sig2.cur
    }
  }
  sims <- cbind(mu.draws, sig2.draws, x.draws)

  return(sims)
}
mu.mh.rw<-function(y,x,mu,sig2){
  mu.accpt.cnt <- 0
  n<-length(y)
  mu.full = function(m){
    ldens = 0
    for(i in 1:n){
      ldens = ldens + log(  pnorm(y[i]+.5, m, sqrt(sig2) ) - pnorm(y[i]-.5, m, sqrt(sig2) )   )
    }
    ldens = ldens - (5/sig2)*(m-mean(x))^2
    return(ldens)
  }
  p.cur = mu.full(mu)
  mu.pro <- exp(log(mu) + rnorm(1, 0, 1))  ##generate a proposed value
  p.pro = mu.full(mu.pro)
  accpt.prob <- exp(p.pro - p.cur)
  if(runif(1) < accpt.prob)
  {
    mu <- mu.pro
    mu.accpt.cnt <- mu.accpt.cnt + 1
  }
  return(mu)
}
sig2.mh.rw = function(y,x,mu,sig2){
  sig2.accpt.cnt <- 0
  n<-length(y)
  sig2.full = function(s2){
    ldens2 = 0
    for(i in 1:n){
      ldens2 = ldens2 + log(pnorm(y[i]+.5,mu,sqrt(s2)) - pnorm(y[i]-.5,mu,sqrt(s2)))
    }
    for(i in 1:n){
      f = 0
      f = f + (x[i]-mu)^2
    }
    ldens2 = ldens2 + (-6)*log(s2) - (1/(2*s2))*f
    return(ldens2)
  }
  p.cur = sig2.full(sig2)
  sig2.pro <- exp(log(sig2) + rnorm(1, 0, 10) ) ##generate a proposed value
  p.pro = sig2.full(sig2.pro)
  accpt.prob <- exp(p.pro - p.cur)
  if(runif(1) < accpt.prob)
  {
    sig2 <- sig2.pro
    sig2.accpt.cnt <- sig2.accpt.cnt + 1
  }
  return(sig2)
}
n.sims <- 30000
y<-c(7,6,7,5,5,3,6,5,4,3)
sample= gibbs(n.sims, y, 1000, 5)

#data samples
mu.mcmc  = sample[,1]
sig2.mcmc = sample[,2]
x.mcmc = cbind(sample[,3],sample[,4],sample[,5],sample[,6],sample[,7],sample[,8],sample[,9],sample[,10])
x.mcmc = cbind(x.mcmc,sample[,11],sample[,12])

hist(sig2.mcmc)
hist(mu.mcmc)

