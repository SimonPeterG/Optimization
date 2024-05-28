rm(list = ls())
set.seed(123)
PCKG = c(
  'Matrix',
  'mvtnorm',
  'doParallel',
  'foreach',
  # 'CVXR',
  'freqdom',
  'complexplus',
  # 'rlist',
  # 'dineR',
  # 'DiffGraph',
  'dplyr',
  'sets'
)

install.packages(setdiff(PCKG, rownames(installed.packages())))
lapply(PCKG, library, character.only = T)
##########################################################################
real.mat = function(rho, d) {
  grid = expand.grid(1:d, 1:d)
  x = apply(grid, 1, function(vec) {
    rho ^ abs(vec[1] - vec[2])
  })
  matrix(x, nrow = d, byrow = T)
}

im.mat = function(rho, d) {
  grid = expand.grid(1:d, 1:d)
  x = apply(grid, 1, function(vec) {
    rho ^ abs(vec[2] - vec[1]) * (-1) ^ (vec[2] > vec[1])
  })
  matrix(x, nrow = d, byrow = F) - diag(d)
}
##########################################################################
D_trace = function(A, B, lambda, rho, maxIter, tolerance) {
  #######Input:
  #A,B: covariance matrices of the data matrices
  #lambda: tuning parameter
  #rho: the coefficient of the augmented Lagrangian in D-trace,
  #usually set to 1
  #######Output:
  #D3 = inv(B) - inv(A) : result of differential matrix
  
  # eA = eigen(A)
  # eB = eigen(B)
  # a1 = matrix(eA$values, ncol = 1)
  # b1 = matrix(eB$values, ncol = 1)
  # a2 = eA$vectors
  # b2 = eB$vectors
  #
  k = ncol(A)
  
  D1 = solve(diag(diag(B)) + diag(k)) - solve(diag(diag(A)) + diag(k))
  # D1 = solve(b2 %*% (diag(c(b1)) + diag(k)) %*% t(b2)) - solve(a2 %*% (diag(c(a1)) +                                                                         diag(k)) %*% t(a2))
  D2 = D1
  D3 = D1
  
  oldD1 = D1 + 0.01 * diag(k)
  
  oldD2 = D2 + 0.01 * diag(k)
  
  oldD3 = D3 + 0.01 * diag(k)
  
  Z = matrix(0, k, k)
  
  L1 = Z
  
  L2 = Z
  
  L3 = Z
  
  iter = 0
  
  
  eA = eigen(A)
  eB = eigen(B)
  a1 = matrix(eA$values, ncol = 1)
  b1 = matrix(eB$values, ncol = 1)
  a2 = eA$vectors
  b2 = eB$vectors
  # a1=svd(A)[[1]]#S
  # a2=svd(A)[[2]]#U
  # b1=svd(B)[[1]]
  # b2=svd(B)[[2]]
  temp1 = (t(a1 %*% t(b1)) + 4 * rho)
  temp2 = (t(b1 %*% t(a1)) + 4 * rho)
  
  S.A.lamb = function(AA, lamb) {
    A1 = matrix(0, nrow(AA), ncol(AA))
    A1 = (AA - lamb) * (AA > lamb) + (AA + lamb) * (AA < (-lamb))
    return(A1)
  }
  while (iter < maxIter &
         {
           norm((D1 - oldD1), "F") / max(1, norm(D1, "F"), norm(oldD1, "F")) > tolerance |
             norm((D2 - oldD2), "F") / max(1, norm(D2, "F"), norm(oldD2, "F")) >
             tolerance |
             norm((D3 - oldD3), "F") / max(1, norm(D3, "F"), norm(oldD3, "F")) >
             tolerance
         }) {
    oldD1 = D1
    oldD2 = D2
    oldD3 = D3
    D1 = a2 %*% (t(a2) %*% ((2 * rho * D3 + 2 * rho * D2 + A - B + 2 * L1 -
                               2 * L3) / temp1) %*% t(b2)) %*% b2
    D2 = b2 %*% (t(b2) %*% ((2 * rho * D3 + 2 * rho * D1 + A - B + 2 * L3 -
                               2 * L2) / temp2) %*% t(a2)) %*% a2
    # D3=sign((rho*D1+rho*D2-L1+L2)/(2*rho))*pmax(abs((rho*D1+rho*D2-L1+L2)/(2*rho))-lambda/(2*rho),0)
    D3 = S.A.lamb(AA = (rho * D1 + rho * D2 - L1 + L2) / (2 * rho),
                  lamb = lambda / (2 * rho))
    L1 = L1 + rho * (D3 - D1)
    L2 = L2 + rho * (D2 - D3)
    L3 = L2 + rho * (D1 - D2)
    iter = iter + 1
    print(iter)
    # print(norm((D1 - oldD1), "F"))
  }
  # out = c(norm((D1 - oldD1), "F") / max(1, norm(D1, "F"), norm(oldD1, "F")),
  #         norm((D2 - oldD2), "F") / max(1, norm(D2, "F"), norm(oldD2, "F")),
  #         norm((D3 - oldD3), "F") / max(1, norm(D3, "F"), norm(oldD3, "F")))
  return(D3)
}
##########################################################################
dft.fun = function(ts.mat, f.val) {
  #ts.mat is the vector values time series d X T
  # f.val is the frequency lying in the set freqs.sipe
  #calculates DFT at a given frequency
  Ma = matrix(
    complex(
      length.out = lenT * d,
      real = cos(2 * pi * f.val * (1:lenT)),
      imaginary = -sin(2 * pi * f.val * (1:lenT))
    ),
    nrow = d,
    ncol = lenT,
    byrow = T
  )
  rowSums(ts.mat * Ma) / sqrt(lenT)
}

# dft.fun(t(T1.test[[1]][,-1]), rel.freq[1])
###########################################################################
simdOmega = result = NULL
ITER = 100
real.rho1 = 0.5
im.rho1 = 0.5
sig1 = 1 / 5
sparse.real.rho2 = 0.5
sparse.im.rho2 = 0.5

N1.train = N2.train = 20 # training sizes of the two class
N1.test = N2.test = 20 # test sizes of the two class
N1 = N1.train + N1.test
N2 = N2.train + N2.test

d = 100 # dimension of the time series
varNo = d * (d + 1) / 2 + d * (d - 1) / 2 # total number of unknown parameters
lenT = 200 #length of the time series

n1.train = n1.test = 20
n2.train = n2.test = 20
n1 = n1.train + n1.test
n2 = n2.train + n2.test

true.train.lab = rep(1:2, c(n1.train, n2.train))
true.test.lab = rep(1:2, c(n1.test, n2.test))

freqs = (1:lenT) * 2 * pi / lenT # let lenT be even
freq.len = length(freqs)

freqs.sipe = seq(1, floor(lenT / 2) - 1, 1) / lenT #seq(0, lenT / (2 * lenT), by = 1 / lenT)  #frequencies where the estimates will be obtained at
len.sipe = length(freqs.sipe)

##########################################################################
# true parameter values that were used to generate the data
F1 = F2 = IF1 = IF2 = D = replicate(lenT, array(dim = c(d, d)), simplify = FALSE)
for (i in 1:(freq.len * 0.5)) {
  IF1.R = real.mat(real.rho1, d)
  IF1.I = im.mat(im.rho1, d)
  # tmp.IF1 = tmp.IF2 = matrix(0,d,d)
  for (k in 1:d) {
    for (j in 1:d) {
      IF1[[i]][k, j] = complex(length.out = 1,
                               real = IF1.R[k, j],
                               imaginary = IF1.I[k, j])
      IF2[[i]][k, j] = IF1[[i]][k, j]
      if (abs(k - j) == 1 & (i <= round(0.1 * lenT))) {
        IF2[[i]][k, j] = -IF1[[i]][k, j]#complex(length.out = 1, real = sparse.real.rho2, imaginary = (sparse.im.rho2)* (-1)^(j>k)) #^abs(k-j)
        # print(k-j); print(i)
      }
    }
  }
  
  IF1[[i]] = sig1 * (IF1[[i]] + diag(d))
  IF2[[i]] = sig1 * (IF2[[i]] + diag(d))
  #-------------------------------------------------------------------
  D[[i]] = IF2[[i]] - IF1[[i]]
  
  
  if (freqs[i] == pi) {
    IF1[[i]] = Re(IF1[[i]])
    IF2[[i]] = Re(IF2[[i]])
  }
  
  F1[[i]] = solve(IF1[[i]]) #spectral density matrix
  F2[[i]] = solve(IF2[[i]]) #spectral density matrix of the second population
  
  if (i %% 5 == 0)
    print(i)
}

for (i in (0.5 * freq.len + 1):(freq.len - 1)) {
  IF1[[i]] = t(IF1[[(freq.len - i)]])
  IF2[[i]] = t(IF2[[(freq.len - i)]])
  F1[[i]] = solve(IF1[[i]]) #spectral density matrix
  F2[[i]] = solve(IF2[[i]]) #spectral density matrix of the second population
}

#---- setting the parameters for the frequency 2pi----
IF1[[freq.len]] = sqrt(abs(IF1[[1]]))

F1[[freq.len]] = solve(IF1[[freq.len]]) #spectral density matrix of the second population

IF2[[freq.len]] = sqrt(abs(IF2[[1]]))
F2[[freq.len]] = solve(IF2[[freq.len]]) #spectral density matrix of the second population


unfold.IF1 = lapply(IF1, function(df1)
  rbind(cbind(Re(df1), -Im(df1)), cbind(Im(df1), Re(df1))))
unfold.IF2 = lapply(IF2, function(df1)
  rbind(cbind(Re(df1), -Im(df1)), cbind(Im(df1), Re(df1))))

gc()
##########################################################################

cl = makeCluster(detectCores() , type = 'FORK')
registerDoParallel(cl)
clusterExport(cl, ls())
##########################################################################
for (itr in 1:ITER) {
  X = read.table(
    './TSdata_d100.txt',
    sep = ',',
    nrows = d * n1,
    skip = (itr - 1) * (n1 + n2) * d + 1,
    header = F
  )
  Y = read.table(
    './TSdata_d100.txt',
    sep = ',',
    nrows = d * n2,
    skip = (itr - 1) * (n1 + n2) * d + d * n1 + 1,
    header = F
  )
  
  X = X[, -c(1, 2)]
  Y = Y[, -c(1, 2)]
  X = split(x = X, f = factor(X$V3))
  Y = split(x = Y, f = factor(Y$V3))
  ############################################################################
  # processing of training data and estimation of spectral density matrices at fundamental Fourier frequencies
  tmp = parLapply(cl, X, function(df)
    spectral.density(t(df[, -1]), freq = freqs.sipe * 2 * pi)$operators) # estimation of spectral density matrices for the first class
  
  S1 = Reduce('+', tmp[1:n1.train]) / n1.train
  # S1 = S1[, , 101:201]
  S1.R = Re(S1)
  S1.R = lapply(1:dim(S1.R)[3], function(val)
    (S1.R[, , val] + t(S1.R[, , val])) / 2)
  S1.I = Im(S1)
  S1.I = lapply(1:dim(S1.I)[3], function(val)
    (S1.I[, , val] - t(S1.I[, , val])) / 2)
  
  # c(sapply(S1.R, isSymmetric),
  #   sapply(S1.I, is.skew.symmetric.matrix))
  S1.test = tmp[n1.train + (1:n1.test)]
  rm(tmp)
  
  # --------------------------------------------------------
  tmp = parLapply(cl, Y, function(df)
    spectral.density(t(df[, -1]), freq = freqs.sipe * 2 * pi)$operators) # estimation of spectral density matrices for the second class
  
  S2 = Reduce('+', tmp[1:n2.train]) / n2.train
  
  S2.R = Re(S2)
  S2.R = lapply(1:dim(S2.R)[3], function(val)
    (S2.R[, , val] + t(S2.R[, , val])) / 2)
  S2.I = Im(S2)
  S2.I = lapply(1:dim(S2.I)[3], function(val)
    (S2.I[, , val] - t(S2.I[, , val])) / 2)
  
  S2.test = tmp[n2.train + (1:n2.test)]
  rm(tmp)
  # plot(freqs.sipe, sapply(1:len.sipe, function(val)
  #   sqrt(sum(abs(
  #     S1[, , val] - S2[, , val]
  #   ) ^ 2))), ty= 'l')
  # 
  #  unfolding the complex coherency matrices -------------------------------------
  Big.S1 = Big.S2 = Big.S = list(NULL)
  for (i in 1:len.sipe) {
    Big.S1[[i]] = cbind(rbind(S1.R[[i]], S1.I[[i]]), rbind(-S1.I[[i]], S1.R[[i]]))
    Big.S2[[i]] = cbind(rbind(S2.R[[i]], S2.I[[i]]), rbind(-S2.I[[i]], S2.R[[i]]))
  }
  clusterExport(cl, c('Big.S1', 'Big.S2'))
  
  ##########################################################################
  # Estiamtion of Dk begins now
  est1.D = parLapply(cl, 1:length(Big.S1), function(val) { #estimates Dk for all k. returns a list of estimated matrices.
    nlambda = 10 # length of lambda sequence
    lambda_min_ratio = 0.05 
    lambda_max <- max(abs(Big.S2[[val]] - Big.S1[[val]])) # this might be some other choice, e.g., a fixed value.
    lambda_min <- lambda_min_ratio * lambda_max
    
    lambdas <-
      exp(seq(log(lambda_min), log(lambda_max), length = nlambda)) # could be some other reasonable sequence
    # lambdas = 2 ^ seq(-1, 5, length.out = nlambda)
    
    ESTs = list(NULL)
    bics = rep(0, nlambda)
    for (indx in 1:nlambda) {
      # .export = c('lambdas', 'Big.S1', 'Big.S2', 'D_trace', 'val')
      # ) %do% {
      obj = D_trace(
        A = Big.S1[[val]],
        B = Big.S2[[val]],
        lambda = lambdas[indx],
        rho = 1,
        maxIter = 5e3,
        tolerance = 1e-3
      )
      ESTs[[indx]] = obj
    
      if(sum(abs(obj))==0){
        break
      }
    }
    
    bics = sapply(1:length(ESTs), function(val1) { #modify the bic criterion
      bic1 = 0.25 * sum(diag((Big.S1[[val]] %*% ESTs[[val1]] %*% Big.S2[[val]] %*% t(ESTs[[val1]])) + (Big.S2[[val]] %*% ESTs[[val1]] %*% Big.S1[[val]] %*% t(ESTs[[val1]]))
      ))
      bic2 = sum(diag(ESTs[[val1]] %*% (Big.S1[[val]] - Big.S2[[val]]))) #+  (log(log(N1.train + N2.train)) * log(d^2) * sum(ESTs[[val1]] != 0))
      
      bic3 = (bic1 - bic2)
      # bic3 = (n1.train + n2.train) * (bic1 - bic2)
      
      bic4 = log(n1.train + n2.train) * sum(ESTs[[val1]] != 0)
      # bic4 =  (log(log(N1.train + N2.train)) * log(4 * (d ^ 2)) * sum(ESTs[[val1]] != 0))
      
      return(bic3 + bic4)
    })
    # bics
    # plot(lambdas, bics,)
    rm.lambda.indx = which(sapply(ESTs, sum) != 0)
    if (length(rm.lambda.indx) > 0) {
      bics.new = bics[rm.lambda.indx]
      ESTs.new = ESTs[rm.lambda.indx]
      i0 = which.min(bics.new)
      obj = ESTs.new[[i0]] #$path[[1]]
    } else{
      i0 = which.min(bics)
      obj = ESTs[[i0]] #$path[[1]]
    }
    
    #now we make the estimate symmetric
    Re.D = (obj[1:d, 1:d] + obj[d + (1:d), d + (1:d)]) / 2
    Im.D = (obj[d + (1:d), 1:d] - obj[(1:d), d + (1:d)]) / 2
    
    Re.D = Re.D * (Re.D <= t(Re.D)) + t(Re.D) * (Re.D > t(Re.D))
    Im.D = Im.D * (Im.D <= t(Im.D)) + t(Im.D) * (Im.D > t(Im.D))
    # return(obj)
    return(rbind(cbind(Re.D, -Im.D), cbind(Im.D, Re.D)))
  })
  # dtrace estimation is complete.
  #------------------------------------------------------------------------
}

# Frobenius norm of the estimated Dk-s
estD.Fnorm = sapply(est1.D, norm, type = 'F')
plot(freqs.sipe,estD.Fnorm, ty='l')

true.Fnorm = sapply(1:length(unfold.IF1), function(val) norm(unfold.IF1[[val]] - unfold.IF2[[val]], type = 'F'))
plot(freqs/(2*pi),true.Fnorm, ty='l')

library(plot.matrix)

plot(unfold.IF1[[1]] - unfold.IF2[[1]]) #true pattern of sparsity
plot(est1.D[[1]]) #estimated pattern

stopCluster(cl)
gc()
