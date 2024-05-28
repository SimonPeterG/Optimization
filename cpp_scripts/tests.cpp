#include <iostream>
#include <armadillo>
#include <math.h>
#include "cov_matrices.h"
#include <complex>

using namespace std;
using namespace arma;

int main() {
    int d = 3;
    cx_mat A(d, d, arma::fill::zeros);
    mat real_mat = real_matrix(0.5, d);
    mat comp_mat = im_matrix(0.5, d);
    A.set_real(real_mat), A.set_imag(comp_mat);
    cout << arma::sqrt(arma::abs(A)) << endl;

    return 0;
}
