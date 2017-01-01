# part a
set.seed(123)
nsim <- 2000
y<-c(7,6,7,5,5,3,6,5,4,3)

mugrid <- seq(0,8,length=200)
logsig2grid <- seq(-1,4,length=200)
contours <- c(.0001, .001, .01, seq(.05, .95, .05))

post.b = function(mu,sig2,y){
  ldens = 0
  for(i in 1:length(y)){
    ldens = ldens + log(pnorm(y[i]+.5,mu,sqrt(sig2)) - pnorm(y[i]-.5,mu,sqrt(sig2)))
  }
  ldens = log(1/sig2) + ldens
  return(ldens)
}

summ <- function(x){
  c(mean(x),sqrt(var(x)), quantile(x, c(.025,.25,.5,.75,.975)))
}

logdens <- outer(mugrid, exp(logsig2grid), post.b, y)
dens <- exp(logdens - max(logdens)) 
contour(mugrid, logsig2grid, dens, levels=contours, 
        xlab='mu', ylab='log sigma')
mtext("Posterior density, accounting for rounding")
dens.mu <- apply(dens,1,sum)
muindex <- sample (1:length(mugrid),nsim, replace=T,prob=dens.mu)
mu.grid <- mugrid[muindex]
sig2.grid <- rep(NA, nsim)
for (i in (1:nsim)){
  sig2.grid[i] <- exp(sample(logsig2grid, 1, prob=dens[muindex[i],]))
}
print(rbind(summ(mu.grid),summ(sig2.grid)))

hist(mu.grid)
hist(sig2.grid)




