arma::mat real_matrix(float rho, int d) {
    /*generates the real part of the precision matrix of the time series
    this matrix is symmetric*/
   arma::mat A(d, d);
    for (int i = 0; i < d; i++) {
        for (int j = i; j < d; j++) {
            A(i ,j) = pow(rho, abs(i - j));
            A(j ,i) = A(i, j);
        }
    }
    return A;
}

arma::mat im_matrix(float rho, int d) {
    /*generates the real part of the precision matrix of the time series
    this matrix is anti-symmetric*/
    arma::mat A(d, d);
    for (int i = 0; i < d; i++) {
        for (int j = i+1; j < d; j++) {
            A(i ,j) = -pow(rho, abs(i - j));
            A(j ,i) = -A(i, j);
        }
    }
    return A;
}

//note: the precision matrices are hermitian since the real part is symmetric and the imaginary part is anti-symmetric