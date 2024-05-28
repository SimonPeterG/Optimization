#include<iostream>
#include<omp.h>
#include<armadillo>
#include<math.h>
#include"cov_matrices.h"

using namespace std;
using namespace arma;

int main() {
    int d = 100;
    double rho = 0.5;
    
    double start_time = omp_get_wtime();

    mat single = real_matrix(rho, d);

    cout << "It took " << omp_get_wtime() << " seconds with the sequential version" << endl;

    start_time = omp_get_wtime();

    mat parallel = paralell_real_matrix(rho, d);

    cout << "It took " << omp_get_wtime() << " seconds with the parallel version" << endl;
    
    return 0;
}