rm(list = ls())
library(plot.matrix)
library(ggplot2)
library(reshape2)

real.mat = function(rho, d) {
  # generates pd symmetric matrix with auto-correlation structure
  grid = expand.grid(1:d, 1:d)
  x = apply(grid, 1, function(vec) {
    rho ^ abs(vec[1] - vec[2])
  })
  matrix(x, nrow = d, byrow = T)
}

d = 100 #dimension of the matrix
invM1 = 100*real.mat(0.9, d) 

invM2 = invM1
for (k in 1:d) {
  for (j in 1:d) {
    if (abs(k - j) == 1) {
      invM2[k, j] = -invM1[k, j] 
    }
  }
}

plot(invM2-invM1) # the true difference

M1 = solve(invM1) #argument1 of D_trace
M2 = solve(invM2) #argument2 of D_trace

nlambda = 20 # length of lambda sequence/ the regularizer sequence
lambda_min_ratio = 0.05 
lambda_max <- max(abs(M2 - M1)) # this might be some other choice, e.g., a fixed value
lambda_min <- lambda_min_ratio * lambda_max
lambdas <-  exp(seq(log(lambda_min), log(lambda_max), length = nlambda))

# DO WE PLAY AROUNF WITH THE LAMBDAS TOO?


Dsolve = list(NULL)
for(i in 1:nlambda){
  Dsolve[[i]] =   D_trace(A = M1, B = M2, lambda = lambdas[i], rho = 1, maxIter = 5e4, tolerance = 1e-4)
  print(i)
}

out = sapply(3:nlambda, function(val) norm(Dsolve[[val]] - (invM2-invM1),type = 'F')) #checking how distant the estimates are from the true value
plot(lambdas[-(1:2)],out, xaxt = 'n')
axis(1, at = lambdas[-(1:2)],labels = round(lambdas[-(1:2)],3))

lnf = reshape2::melt(invM2-invM1)
head(lnf)
p0 = ggplot(lnf, aes(x=Var1, y=rev(Var2), fill=value)) + #the true difference
  geom_tile(colour="white") +
    scale_fill_gradientn(colours = terrain.colors(n = 10), breaks = seq(-400, 1000, 100), limits = c(-400, 1000))+
  scale_x_continuous(name = element_blank())+
  scale_y_continuous(name = element_blank())+
  theme(panel.background = element_blank())

p <- list()
for(i in 1:nlambda){ # the estimated difference
  p[[i]] <- ggplot(melt(Dsolve[[i]]), aes(x=Var1, y=rev(Var2), fill=value)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(colours = terrain.colors(n = 10), breaks = seq(-400, 1000, 100), limits = c(-400, 1000))+
    scale_x_continuous(name = element_blank())+
    scale_y_continuous(name = element_blank())+
    theme(panel.background = element_blank())
}

egg::ggarrange(p0, p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], nrow = 3, ncol = 4, byrow = T)
