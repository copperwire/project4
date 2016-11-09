#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include "lib.h"

using namespace std;

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

arma::imat initialize (arma::imat matrix, int n_spins, double& E, double& M)
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

arma::imat Metropolis(int n_spins, arma::imat spin_matrix, double& E, double& M, double* w, long& idum)
{
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (ran1(&idum) * (double) n_spins);
            int iy = (int) (ran1(&idum) * (double) n_spins);

            cout << ix << "  " << iy << endl;
            int deltaE = 2*spin_matrix(iy, ix)*
                                             (spin_matrix(iy, periodic(ix,n_spins,-1))+
                                              spin_matrix(periodic(iy,n_spins,-1), ix) +
                                              spin_matrix(iy, periodic(ix,n_spins,1)) +
                                              spin_matrix(periodic(iy,n_spins,1), ix));
            if(deltaE <= 0)
            {
                //cout << "-----------------interim--------------" << endl;
                spin_matrix(iy, ix) *= -1; // flip one spin and accept new spin config
                //spin_matrix.print();
                //cout << "-----------------exterim--------------" << endl;
                // update energy and magnetization
                M += (double) 2*spin_matrix(iy, ix);
                E += (double) deltaE;
            }
            else if(deltaE > 0)
            {
                double r = ran1(&idum);
                if ( r <= w[deltaE+8] )
                {
                    spin_matrix(iy, ix) *= -1; // flip one spin and accept new spin config
                    // update energy and magnetization
                    M += (double) 2*spin_matrix(iy, ix);
                    E += (double) deltaE;
                }
            }
        }
    }
    spin_matrix.print();
    cout << "--------------hello-------------" << endl;
    return spin_matrix;
}

void output(arma::imat spin_matrix, double* average, int n_spins, int m_carlos)
{
    double over_iter = 1/(double) m_carlos;
    double E = average[0] * over_iter;
    double EE = average[1] * over_iter;

    double M = average[2] * over_iter;
    double MM = average[3] * over_iter;
    double M_abs = average[4] * over_iter;

    double E_var = (EE - E*E )/n_spins/n_spins;
    double M_var = (MM - M_abs * M_abs)/n_spins/n_spins;

    cout << "Variance in energy: " << E_var << endl;
    cout << "Variance in magnetization: " << M_var << endl;
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
    spin_matrix.zeros();

    long idum = -1;

    for(double temp = init_temp; temp <= final_temp ; temp += step_temp){
        E = 0;
        M = 0;

        for(int de = -8; de <= 8; de ++) w[de+8] = 0;
        for(int de = -8; de <= 8; de +=4) w[de+8] = exp(-de/temp);
        for(int i = 0; i<5 ; i++ ) average[i] =0;
        spin_matrix = initialize (spin_matrix,  n_spins,  E,  M);

        for(int i = 0 ; i<=max_carlos; i++)
        {
            spin_matrix = Metropolis(n_spins, spin_matrix, E, M, w, idum);

            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
        output(spin_matrix, average, n_spins, max_carlos);
    }
}
