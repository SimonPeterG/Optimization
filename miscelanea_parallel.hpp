arma::cx_mat gen_norm(int d, arma::mat sigma, int n = 1) {
    arma::mat Z = arma::randn(n, d);
    arma::mat R = arma::chol(sigma);
    arma::cx_mat X(n, d);
    X.set_real(Z * R);
    return X;
}

arma::cx_mat rbind(std::vector<arma::cx_mat> matrices) {
  arma::cx_mat X = matrices[0];
  for (int i = 1; i < matrices.size(); i++) {
    X = arma::join_cols(X, matrices[i]);
  }
  return X;
}

std::vector<arma::cx_vec> gen_spectral_obs(int lenT, std::vector<arma::cx_mat> F, arma::cx_mat Z) {
    //ask wheter freqlen and lenT can be different or not
    std::vector<arma::cx_vec> X(lenT);
    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::mat lambda;
    #pragma omp parallel for
        for (int i = 0; i < lenT; i++) {
            if(lenT/2 <= i & i < lenT - 1) {
                int idx = lenT - 2 - i;
                arma::eig_sym(eigval, eigvec, F[idx]);
                lambda = arma::diagmat(eigval);
                X[i] = arma::conj(eigvec) * arma::sqrt(lambda) * arma::conj(Z.row(idx).st());
            } else {
                arma::eig_sym(eigval, eigvec, F[i]);
                lambda = arma::diagmat(eigval);
                X[i] = eigvec * arma::sqrt(lambda) * Z.row(i).st();
            }
        }
    return X;
}


arma::cx_vec gen_ts_obs(int d, int lenT, int val, std::vector<arma::cx_vec> spectral_obs, double freqs[]) {
    arma::cx_vec X(d);
    
    #pragma omp parallel for reduction(+:X)
        for (int i = 0; i < lenT; i++) {
            std::complex<double> z(cos(val * freqs[i]), sin(val * freqs[i]));
            X += z * spectral_obs[i];
        }
    return X / sqrt(2 * lenT);
}

arma::mat gen_time_series(int d, int lenT, std::vector<arma::cx_mat> F, double freqs[]) {
    arma::cx_mat ZTby2  = gen_norm(d, 2 * arma::eye(d, d));
    arma::cx_mat ZT = gen_norm(d, 2 * arma::eye(d, d));
    arma::cx_mat Ztmp = gen_norm(2 * d, arma::eye(2 * d, 2 * d), lenT/2 - 1);
    arma::cx_mat Z_half1(lenT/2 - 1, d);
    Z_half1.set_real(arma::real(Ztmp.cols(0, d - 1)));
    Z_half1.set_imag(arma::real(Ztmp.cols(d, 2 * d - 1)));
    arma::uvec idx(lenT/2 - 1);
    for (int i = 0; i < lenT/2 - 1; i++) {
        idx(i) = lenT/2 - 1 - (i+1);
    }
    arma::cx_mat Z_half2 = arma::conj(Z_half1);
    arma::cx_mat Z = rbind(std::vector<arma::cx_mat> {Z_half1, ZTby2, Z_half2.rows(idx), ZT});
    
    std::vector<arma::cx_vec> spectral_obs = gen_spectral_obs(lenT, F, Z);
    arma::cx_mat ts_realization(d, lenT);
    for (int i = 0; i < lenT; i++) {
        ts_realization.col(i) = gen_ts_obs(d, lenT, i, spectral_obs, freqs);
    }
    return arma::real(ts_realization);
    
}

void write_data(std::string filename , int d, int iter, int label, int sampleid, arma::mat data) {
    arma::mat left_mat(d, 3);
    left_mat.col(0).fill(iter);
    left_mat.col(1).fill(label);
    left_mat.col(2).fill(sampleid);
    data = arma::join_rows(left_mat, data);


    std::stringstream buffer;
    buffer << std::fixed << std::setprecision(16);
    data.raw_print(buffer);

    std::ofstream file;
    file.open(filename, std::ios_base::app);
    file << buffer.str();
    file.close();
}