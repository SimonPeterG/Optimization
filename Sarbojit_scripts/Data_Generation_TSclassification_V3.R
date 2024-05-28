# Simulating data from two d X 1 dimensional time series of length lenT
rm(list = ls())
init.time = Sys.time() #initial timing
set.seed(135)
PCKG = c(
  'Matrix',
  'mvtnorm',
  'doParallel',
  'foreach',
  'freqdom',
  'complexplus'
)

install.packages(setdiff(PCKG, rownames(installed.packages())))
lapply(PCKG, library, character.only = T)
##########################################################################
real.mat = function(rho, d){ # generates pd symmetric matrix with auto-correlation structure
  grid = expand.grid(1:d,1:d)
  x = apply(grid, 1, function(vec){ #generates matrix with entries rho^|i-j| (AR(1))
    rho^abs(vec[1]-vec[2])
  })
  matrix(x, nrow = d, byrow = T)
}

im.mat = function(rho, d){ # generates skew-symmetric matrix with auto-correlation structure
  grid = expand.grid(1:d,1:d)
  x = apply(grid, 1, function(vec){
    rho^abs(vec[2]-vec[1]) * (-1)^(vec[2]>vec[1])
  }) 
  matrix(x, nrow = d, byrow = F) - diag(d) 
  #generates matrix with entries rho^|i-j|*sgn(j-I) 
  #positive below main diagonal, negative above main diagonal
}
##########################################################################
real.rho1 = 0.5
im.rho1 = 0.5
sig1 = 1/5
sparse.real.rho2 = 0.5
sparse.im.rho2 = 0.5

ITER = 5 #100 #number of iterations.
N1.train = N2.train = 20 # training sizes of the two class
N1.test = N2.test = 20 # test sizes of the two class
N1 = N1.train + N1.test
N2 = N2.train + N2.test

true.train.lab = rep(1:2, c(N1.train, N2.train))
true.test.lab = rep(1:2, c(N1.test, N2.test))

d = 100 # dimension of the time series changed
varNo = d*(d+1) / 2 + d*(d-1)/2 # total number of unknown parameters
lenT = 50 #length of the time series changed

freqs = (1:lenT) * 2 * pi / lenT # let lenT be even
freq.len = length(freqs)

freqs.sipe = seq(-lenT/(2*lenT), lenT/(2*lenT), by = 1/lenT) #frequencies where the estimates will be obtained at
len.sipe = length(freqs.sipe)
##########################################################################
# parameters of the simulation exercise are constructed below
F1 = F2 = IF1 = IF2 = D = replicate(lenT, array(dim = c(d, d)), simplify = FALSE)
#filling the first half of the lists
for (i in 1:(freq.len * 0.5)) {
  IF1.R = real.mat(real.rho1, d) #creating real part of the matrix
  IF1.I = im.mat(im.rho1, d) #creating imaginary part of the matrix
  # tmp.IF1 = tmp.IF2 = matrix(0,d,d)
  #filling the entries of the ith IF1 and IF2 matrix
  for(k in 1:d){
    for(j in 1:d){
      #filling the entries of both matrices based on the real and imaginary matrices computed above
      IF1[[i]][k,j] = complex(length.out = 1, real = IF1.R[k,j], imaginary = IF1.I[k,j])
      IF2[[i]][k,j] = IF1[[i]][k,j]
      if(abs(k-j)==1 & (i <= round(0.1 * lenT))){
        #changing the sign of the elements immediately to the right (also left) and below (also above)
        #the main diagonal and if it is an observation below the 10% of the total observations
        IF2[[i]][k,j] = -IF1[[i]][k,j]
        #complex(length.out = 1, real = sparse.real.rho2, imaginary = (sparse.im.rho2)* (-1)^(j>k)) #^abs(k-j)
        # print(k-j); print(i)
      }
    }
  }
  #note that the I on the beginning stands for inverse
  IF1[[i]] = sig1 * (IF1[[i]]+diag(d)) 
  IF2[[i]] = sig1 * (IF2[[i]]+diag(d))
  #-------------------------------------------------------------------
  D[[i]] = IF2[[i]] - IF1[[i]]; #difference matrix on the ith slot of the list
  
  if (freqs[i] == pi) { #if the frequency is pi then the imaginary part is zero
    IF1[[i]] = Re(IF1[[i]])
    IF2[[i]] = Re(IF2[[i]])
  }
  #it is inverting the inverse of Fi, ie, getting Fi
  F1[[i]] = solve(IF1[[i]]) #spectral density matrix
  F2[[i]] = solve(IF2[[i]]) #spectral density matrix of the second population
  
  # if (i %% 5 == 0)
  #   print(i)
}

#filling the rest of the lists (except for the boundary of the spectrum)
for (i in (0.5 * freq.len + 1):(freq.len - 1)) {
  #filling up the second half of the list such that
  #the freq.len/2 + i is the same as the transpose freq.len/2 - i
  IF1[[i]] = t(IF1[[(freq.len - i)]])
  IF2[[i]] = t(IF2[[(freq.len - i)]])
  F1[[i]] = solve(IF1[[i]]) #spectral density matrix of the first population
  F2[[i]] = solve(IF2[[i]]) #spectral density matrix of the second population
}

#---- setting the parameters for the frequency 2pi----
IF1[[freq.len]] = sqrt(abs(IF1[[1]]))

F1[[freq.len]] = solve(IF1[[freq.len]]) #spectral density matrix of the first population

IF2[[freq.len]] = sqrt(abs(IF2[[1]])) 
F2[[freq.len]] = solve(IF2[[freq.len]]) #spectral density matrix of the second population

#The function generates a bloc matrix as follows.
# For a complex valued matrix A it generates the following block matrix
# Re(A) -Im(A)
# Im(A) Re(A)
unfold.IF1 = lapply(IF1, function(df1)
  rbind(cbind(Re(df1),-Im(df1)), cbind(Im(df1), Re(df1))))
unfold.IF2 = lapply(IF2, function(df1)
  rbind(cbind(Re(df1),-Im(df1)), cbind(Im(df1), Re(df1))))

#plotting the frobenius distance between the ith matrix of population 1
#and the ith matrix of population 2

plot(freqs / (2 * pi), sapply(1:freq.len, function(val)
  sqrt(sum(abs(IF1[[val]] - IF2[[val]])^2))), type = "l")
plot(freqs / (2 * pi), sapply(1:freq.len, function(val)
  sqrt(sum(abs(F1[[val]] - F2[[val]])^2))), type = "l")

barplot(height = sapply(c(freq.len,1:(freq.len-1)), function(val) #changed
  sqrt(sum(abs(IF1[[val]] - IF2[[val]])^2))))

gc()
#-----------------------------------------------------------------
cl = makeCluster(detectCores() , type = 'FORK')
registerDoParallel(cl)
clusterExport(cl, ls())

# data generation starts now.


#the data generation procedure is as developed in /Users/roys/Library/CloudStorage/Dropbox/KAUST/Research Proposal/Classification of Time Series/DataGenerationGivenSpecDensityMatrix.pdf
err = NULL
for (itr in 1:ITER) {
  X1 = rep(list(0), N1) #each entry of the list is a d-variate time series of lenght lenT
  X2 = rep(list(0), N2) #each entry is a matrix of dimension dxlenT which is one realization of the time series
  
  for (n1 in 1:N1) {
    Z.Tby2 = rmvnorm(n = 1,
                     mean = rep(0, d),
                     sigma = 2 * diag(d)) #d-variate normal zero mean, 2I covariance goes on the 25th
    Z.T = rmvnorm(n = 1,
                  mean = rep(0, d),
                  sigma = 2 * diag(d)) #d-variate normal zero mean, 2I covariance goes on the 50th
    Z.tmp = rmvnorm(n = lenT / 2 - 1,
                    mean = rep(0, 2 * d),
                    sigma = diag(2 * d)) #2d-variate normal zero mean, I covariance (lenT/2 - 1) observations
    Z.half1 = t(apply(Z.tmp, 1, function(vec) {
      complex(real = vec[1:d], imaginary = vec[d + (1:d)])
    })) #lenT/2 - 1 observations from a complex d-variate normal
    
    Z.half2 = t(apply(Z.tmp, 1, function(vec) {
      complex(real = vec[1:d], imaginary = -vec[d + (1:d)])
    })) #complex conjugate of Z.half1
    #vertically stacking Z.half1, Z.tby2, the reverse of the complex conjugate of Z.half1 and Z.T
    Z = rbind(Z.half1, Z.Tby2, Z.half2[nrow(Z.half2):1, ], Z.T)
    
    rm(list = c('Z.half1', 'Z.half2', 'Z.tmp', 'Z.T', 'Z.Tby2'))
    
    obj = lapply(1:freq.len, function(val) {
      if ((lenT / 2) + 1 <= val && val <= lenT - 1) {
        tmp = eigen(F1[[lenT - val]]) #if its part of the second half then the covariance matrix is conj(C)
        #if C = UDU^H then conj(C) = conj(U)DU^T
        #thus to simulate we take conj(U) D^1/2 conj(Z)
        Conj(tmp$vectors) %*% diag(sqrt(tmp$values)) %*% Conj(Z[lenT - val, ])
      } else{
        tmp = eigen(F1[[val]])
        #To have covariance matrix C = UDU^H we compute UD^1/2 Z
        tmp$vectors %*% diag(sqrt(tmp$values)) %*% Z[val, ]
      }
    })
    
    rm(Z)
    clusterExport(cl, 'obj')
    X1[[n1]] = do.call('cbind', parLapply(cl, 1:lenT, function(val) {
      # X1[[n1]] = do.call('cbind', lapply(1:lenT, function(val) {
      tmp = lapply(1:freq.len, function(Vind) {
        obj[[Vind]] * complex(
          length.out = 1,
          real = cos(val * freqs[Vind]),
          imaginary = sin(val * freqs[Vind])
        )
      }) #it is multiplying each spectral observation times exp(j*wk) for fixed j and all wk
      # Reduce('+', tmp) * sqrt(pi / lenT)
      Reduce('+', tmp) * sqrt(1 / (2 * lenT))
    }))
    X1[[n1]] = Re(X1[[n1]])
    rm(obj)
  }
  
  for (n2 in 1:N2) {
    Z.Tby2 = rmvnorm(n = 1,
                     mean = rep(0, d),
                     sigma = 2 * diag(d))
    Z.T = rmvnorm(n = 1,
                  mean = rep(0, d),
                  sigma = 2 * diag(d))
    Z.tmp = rmvnorm(n = lenT / 2 - 1,
                    mean = rep(0, 2 * d), 
                    # mean = c(rep(5, 10), rep(0, d-10), rep(5, 10) ,rep(0, d-10)), # for location model: change the mean here. the first d components corresponds to real part, next d correspond to imaginary part.
                    sigma = diag(2 * d))
    Z.half1 = t(apply(Z.tmp, 1, function(vec) {
      complex(real = vec[1:d], imaginary = vec[d + (1:d)])
    }))
    
    Z.half2 = t(apply(Z.tmp, 1, function(vec) {
      complex(real = vec[1:d], imaginary = -vec[d + (1:d)])
    }))
    Z = rbind(Z.half1, Z.Tby2, Z.half2[nrow(Z.half2):1, ], Z.T)
    
    rm(list = c('Z.half1', 'Z.half2', 'Z.tmp', 'Z.T', 'Z.Tby2'))
    
    obj = lapply(1:freq.len, function(val) {
      if ((lenT / 2) + 1 <= val && val <= lenT - 1) {
        tmp = eigen(F2[[lenT - val]])
        Conj(tmp$vectors) %*% diag(sqrt(tmp$values)) %*% Conj(Z[lenT - val, ])
      } else{
        tmp = eigen(F2[[val]])
        tmp$vectors %*% diag(sqrt(tmp$values)) %*% Z[val, ]
      }
    })
    rm(Z)
    clusterExport(cl, 'obj')
    X2[[n2]] = do.call('cbind', parLapply(cl, 1:lenT, function(val) {
      # X2[[n2]] = do.call('cbind', lapply(1:lenT, function(val) {
      tmp = lapply(1:freq.len, function(Vind) {
        obj[[Vind]] * complex(
          length.out = 1,
          real = cos(val * freqs[Vind]),
          imaginary = sin(val * freqs[Vind])
        )
      })
      # Reduce('+', tmp) * sqrt(pi / lenT) #[-c(lenT/2, lenT)]
      Reduce('+', tmp) * sqrt(1 / (2 * lenT))
    }))
    X2[[n2]] = Re(X2[[n2]])
    rm(obj)
  }
  #data generation is complete.
  if(itr==1){
    tmp = data.frame('Iteration' = rep(itr, (N1 + N2)*d), 'Label' = c(rep(1, N1*d), rep(2, N2*d)),'sampleID' = c(rep(1:N1, each = d), rep(1:N2, each = d)), rbind(do.call('rbind', X1), do.call('rbind', X2)))
    write.table(x = tmp, file = './TSdata_d100_V2.txt', sep = ',', row.names = F)
  }else{
    tmp = cbind(rep(itr, (N1 + N2)*d), c(rep(1, N1*d), rep(2, N2*d)),c(rep(1:N1, each = d), rep(1:N2, each = d)), rbind(do.call('rbind', X1), do.call('rbind', X2)))
    write.table(x = tmp, file = './TSdata_d100_V2.txt', sep = ',', row.names = F, append = T, col.names = F)
  }
  rm(tmp)
  if (itr %% 5 == 0)
    print(paste('iteration: ', itr, sep = ''))
}
print(paste("We good in", round(Sys.time() - init.time,4), "seconds"))
stopCluster(cl)
gc()
