#include <iostream>
#include <armadillo>
#include <math.h>
#include "cov_matrices.h"
#include <complex>

using namespace std;
using namespace arma;

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

int main()
  {
    int d = 3, lenT = 10;

    // arma::cx_mat ZTby2  = gen_norm(d, 2 * arma::eye(d, d));
    // arma::cx_mat ZT = gen_norm(d, 2 * arma::eye(d, d));
    // arma::cx_mat Ztmp = gen_norm(2 * d, arma::eye(2 * d, 2 * d), lenT/2 - 1);
    // arma::cx_mat Z_half1(lenT/2 - 1, d);
    // Z_half1.set_real(arma::real(Ztmp.cols(0, d - 1)));
    // Z_half1.set_imag(arma::real(Ztmp.cols(d, 2 * d - 1)));
    // arma::uvec idx(lenT/2 - 1);
    // for (int i = 0; i < lenT/2 - 1; i++) {
    //     idx(i) = lenT/2 - 1 - (i+1);
    // }
    // arma::cx_mat Z_half2 = arma::conj(Z_half1);
    // arma::cx_mat Z = rbind(std::vector<arma::cx_mat> {Z_half1, ZTby2, Z_half2.rows(idx), ZT});
    // arma::cx_mat Zt = Z.st();
    // cout << Z << endl;
    // cout << Z.row(0).st() << endl;

    arma::mat sizas(d, d);
    arma::mat left_mat(d, 3);
    left_mat.col(0).fill(1);
    left_mat.col(1).fill(1);
    left_mat.col(2).fill(3);
    sizas = arma::join_rows(left_mat, sizas);
    cout << sizas << endl;
    return 0;
  }