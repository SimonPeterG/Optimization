library(Rcpp)
sourceCpp("./data_generationr.cpp")
a1 <- Sys.time()
gen()
print(Sys.time() - a1)



