#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include "lib.h"
#include <armadillo>


using namespace std;

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

arma::imat initialize(arma::imat matrix, int n_spins, double& E, double& M)
{
    srand (time(NULL));
    for(int i = 0; i< n_spins; i++){
        for(int j = 0; j< n_spins; j++)
        {
            int val = rand() %2;
            if (val == 1)
            {
                matrix(i, j) = 1;
                M += 1;
            }
            else
            {
                matrix(i, j) = -1;
                M += -1;
            }
        }
    }

    for(int i = 0; i< n_spins; i++){
        for(int j = 0; j< n_spins; j++)
        {
            E -= (double) matrix(i, j)*
            (matrix(periodic(i, n_spins, -1), j) +
             matrix(i, periodic(j, n_spins, -1))
            );
        }
    }
    return matrix;
}

arma::imat Metropolis(int n_spins, arma::imat spin_matrix, double& E, double& M, double* w)
{
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
    // Find random position
            int ix = (int) rand() % n_spins;
            int iy = (int) rand() % n_spins;
            int deltaE = 2*spin_matrix(iy, ix)*
                                             (spin_matrix(iy, periodic(ix,n_spins,-1))+
                                              spin_matrix(periodic(iy,n_spins,-1), ix) +
                                              spin_matrix(iy, periodic(ix,n_spins,1)) +
                                              spin_matrix(periodic(iy,n_spins,1), ix));
            // Here we perform the Metropolis test

            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            if ( r <= w[deltaE+8] ) {
                spin_matrix(iy, ix) *= -1; // flip one spin and accept new spin config
                // update energy and magnetization
                M += (double) 2*spin_matrix(iy, ix);
                E += (double) deltaE;
            }
        }
    }

    return spin_matrix;
}

int main(int argc, char *argv[])
{
    string output_filename;
    int n_spins, max_carlos;
    double w[17], average[5];
    double init_temp = 1;
    double final_temp = 1;
    double step_temp = 0.1;
    double E, M  = 0;

    if (argc < 3){
        cout<< "wrong usage" << endl;
        exit(1);
    }
    else{
        n_spins = atoi(argv[1]);
        max_carlos = atoi(argv[2]);
        output_filename = string(argv[3]);
    }

    arma::imat spin_matrix(n_spins, n_spins);

    for(double temp = init_temp; temp <= final_temp ; temp += step_temp){
        E = 0;
        M = 0;

        for(int de = -8; de <= 8; de ++) w[de+8] = 0;
        for(int de = -8; de <= 8; de +=4) w[de+8] = exp(-de/temp);


        for(int i = 0; i<5 ; i++ ) average[i] =0;

        spin_matrix = initialize(spin_matrix, n_spins, E, M);
        spin_matrix.print();
        for(int i = 0 ; i<=max_carlos; i++)
        {
            spin_matrix = Metropolis(n_spins, spin_matrix, E, M, w);
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
    }

}
