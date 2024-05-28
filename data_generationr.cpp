#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <iomanip>
#include <math.h>
#include "cov_matrices.hpp"
#include "miscelanea.hpp"
#include <vector>
#include <complex>

using namespace std;
using namespace arma;

// [[Rcpp::export]]

void gen() {
    double real_rho1, im_rho1, sig1, sparse_real_rho2, sparse_im_rho2;
    int iter, n1_train, n2_train, n1_test, n2_test, n1, n2, d, varNo, lenT, len, freq_len, len_sipe;
    vector<cx_mat> F1, F2, IF1, IF2, D;

    double pi = M_PI;

    real_rho1 = im_rho1 = 0.5;
    sig1 = (double) 1/5;
    sparse_real_rho2 = sparse_im_rho2 = 0.5;

    iter = 6;
    n1_train = n2_train = 20;
    n1_test = n2_test = 20;
    n1 = n1_train + n1_test;
    n2 = n2_train + n2_test;

    int true_train_lab[n1_train + n2_train], true_test_lab[n1_train + n2_train];
    
    for (int i = 0; i < n1_train + n2_train; i++) {
        true_train_lab[i] = (i < n1_train) ? 1 : 2;
    }

    for (int i = 0; i < n1_test + n2_test; i++) {
        true_test_lab[i] = (i < n1_test) ? 1 : 2;
    }

    d = 5;
    varNo = d*(d+1)/2 + d*(d-1)/2;
    lenT = 10;

    double freqs[lenT], freqs_sipe[lenT + 1];
    cx_mat complex_zero(d, d);
    mat real_zero(d, d);

    for (int i = 0; i < lenT; i++) {
        freqs[i] = (i+1)*2*pi/lenT;
    }
    
    freq_len = lenT;

    for (int i = 0; i < lenT + 1; i++) {
        freqs_sipe[i] = -lenT/(2*lenT) + i/lenT;
    }

    len_sipe = lenT + 1;

    for (int i = 0; i < lenT; i++) {
        F1.push_back(complex_zero);
        F2.push_back(complex_zero);
        IF1.push_back(complex_zero);
        IF2.push_back(complex_zero);
        D.push_back(complex_zero);
    }

    for (int i = 0; i < (lenT/2); i++) {
        IF1[i].set_real(real_matrix(real_rho1, d));
        IF1[i].set_imag(im_matrix(real_rho1, d));
        IF2[i] = IF1[i];
        for (int k = 0; k < d; k++) {
            if (i <= 0.1*lenT) {
                if (k == 0) IF2[i](k, k+1) = -IF1[i](k, k+1);
                else if (k == d-1) IF2[i](k, k-1) = -IF1[i](k, k-1);
                else {
                    IF2[i](k, k-1) = -IF1[i](k, k-1);
                    IF2[i](k, k+1) = -IF1[i](k, k+1);
                }
            } 
        }


        IF1[i] = sig1 * (IF1[i] + arma::eye(d,d));

        IF2[i] = sig1 * (IF2[i] + arma::eye(d,d));

        D[i] = IF1[i] - IF2[i];

        if (i == lenT/2 - 1) { //pi will always be at the position lenT/2
            IF1[i].set_imag(real_zero);
            IF2[i].set_imag(real_zero);
        }

        F1[i] = arma::inv(IF1[i]);
        F2[i] = arma::inv(IF2[i]);
    }
    
    for (int i = (lenT/2 + 1); i < lenT - 1; i++) {
        IF1[i] = IF1[freq_len - i].st();
        IF2[i] = IF2[freq_len - i].st();
        F1[i] = arma::inv(IF1[i]);
        F2[i] = arma::inv(IF2[i]);
    }

    IF1[lenT - 1].set_real(arma::sqrt(arma::abs(IF1[1])));
    F1[lenT - 1] = arma::inv(IF1[lenT - 1]);

    IF2[lenT - 1].set_real(arma::sqrt(arma::abs(IF2[1])));
    F2[lenT - 1] = arma::inv(IF2[lenT - 1]);

    mat ts_obs;
    for (int i = 1; i <= iter; i++) {
        for (int sample_number = 1; sample_number <= n1; sample_number++) {
            ts_obs = gen_time_series(d, lenT, F1, freqs);
            write_data("aver.txt", d, i, 1, sample_number, ts_obs);    
        }
        for (int sample_number = 1; sample_number <= n1; sample_number++) {
            ts_obs = gen_time_series(d, lenT, F1, freqs);
            write_data("aver.txt", d, i, 2, sample_number, ts_obs);    
        }
        if (i % 5 == 0) {
            cout << "iteration " << i << endl;
        }
    }
    
    cout << "We good" << endl;
}