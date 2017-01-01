#part (c)

set.seed(123)
source("/%PATH%/rounded_data_ptA.R")
source("/%PATH%/rounded_data_ptB.R")


quantile(mu.grid,c(.025,.975))
#     2.5%    97.5% 
#  4.180905 5.989950 
quantile(sig2.grid,c(.025,.975))
#    2.5%     97.5% 
#  0.7434287 4.6567805 
quantile(mu.mcmc,c(.025,.975))
#   2.5%    97.5% 
#  4.614078 5.571675  
quantile(sig2.mcmc,c(.025,.975))
#   2.5%     97.5% 
#  0.5194269 2.0457027 

quantile(x.mcmc[,1],c(.025,.975))
#     2.5%    97.5% 
#  6.511127 7.430978  
quantile(x.mcmc[,2],c(.025,.975))
#     2.5%    97.5% 
#  5.517410 6.458669 
quantile(x.mcmc[,3],c(.025,.975))
#     2.5%    97.5% 
#  6.512171 7.431616 
quantile(x.mcmc[,4],c(.025,.975))
#   2.5%    97.5% 
#  4.530237 5.473577 
quantile(x.mcmc[,5],c(.025,.975))
#   2.5%    97.5% 
#  4.524440 5.474056 
quantile(x.mcmc[,6],c(.025,.975))
#   2.5%    97.5% 
#  2.579688 3.490216 
quantile(x.mcmc[,7],c(.025,.975))
#   2.5%    97.5% 
#  5.520413 6.456466 
quantile(x.mcmc[,8],c(.025,.975))
#   2.5%    97.5% 
#  4.528392 5.479156 
quantile(x.mcmc[,9],c(.025,.975))
#    2.5%    97.5% 
# 3.552256 4.485055 
quantile(x.mcmc[,10],c(.025,.975))
#   2.5%    97.5% 
#  2.580537 3.486689 

mean(x.mcmc[,1])
# 6.846235
mean(x.mcmc[,2])
# 5.928351
mean(x.mcmc[,3])
# 6.847636
mean(x.mcmc[,4])
# 5.011024
mean(x.mcmc[,5])
# 5.006686
mean(x.mcmc[,6])
# 3.160727
mean(x.mcmc[,7])
# 5.924822
mean(x.mcmc[,8])
# 5.010117
mean(x.mcmc[,9])
# 4.090461
mean(x.mcmc[,10])
# 3.157286

mean(mu.mcmc)
# 5.087816
mean(sig2.mcmc)
# 1.047941
mean(mu.grid)
# 5.090211
mean(sig2.grid)
# 1.955827




# DIC

diclikl = function(mu,sig2,y){
  dens = 1
  for(i in 1:length(y)){
    dens = dens * (pnorm(y[i]+.5,mu,sqrt(sig2)) - pnorm(y[i]-.5,mu,sqrt(sig2)))
  }
  return(dens)
}


all = log(diclikl(mean(mu.mcmc), mean(sig2.mcmc),y))
bll = log(diclikl(mean(mu.grid), mean(sig2.grid),y))

a.est = (1/length(mu.mcmc))*sum( log( diclikl(mu.mcmc,sig2.mcmc,y)  )   )
b.est = (1/length(mu.grid))*sum( log( diclikl(mu.grid,sig2.grid,y)  )   )
a.pdic = 2*(all - a.est)
b.pdic = 2*(bll - b.est)
a.dic = -2*a.est + 2*a.pdic
b.dic = -2*b.est + 2*b.pdic

#> a.dic
#[1] 41.35351
#> b.dic
#[1] 40.75112



