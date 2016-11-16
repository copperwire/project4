#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include "lib.h"

using namespace std;

ofstream ofile;

inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

arma::imat initialize (arma::imat matrix, int n_spins, double& E, double& M, string type)
{
    string alt_1 = "random";
    string alt_2 = "order";
    if(type == alt_1){
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
    else if(type == alt_2){

        for(int i = 0; i< n_spins; i++){
            for(int j = 0; j< n_spins; j++){
            matrix(i, j) = 1;
            M += 1;
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
}

arma::imat Metropolis(int n_spins, arma::imat spin_matrix, double& E, double& M, double* w, int& counter_accept)
{
    // loop over all spins
    //std::random_device rd;
    //std::mt19937_64 gen(rd());


    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            // Set up the uniform distribution for x \in [0, 1]

            //std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

            double r1 = ((double) rand() / (RAND_MAX)) ;
            double r2 = ((double) rand() / (RAND_MAX)) ;
            double r3 = ((double) rand() / (RAND_MAX)) ;

            int ix = (int) (r1 * (double) (n_spins -1));
            int iy = (int) (r2 * (double) (n_spins -1));

            int deltaE = 2*spin_matrix(iy, ix)*
                                             (spin_matrix(iy, periodic(ix,n_spins,-1))+
                                              spin_matrix(periodic(iy,n_spins,-1), ix) +
                                              spin_matrix(iy, periodic(ix,n_spins,1)) +
                                              spin_matrix(periodic(iy,n_spins,1), ix));
            //double r = RandomNumberGenerator(gen);
            if ( r3 <= w[deltaE+8] )
            {
                spin_matrix(iy, ix) *= -1; // flip one spin and accept new spin config
                // update energy and magnetization
                counter_accept += 1;
                M += (double) 2*spin_matrix(iy, ix);
                E += (double) deltaE;
            }

        }
    }
    return spin_matrix;
}

void output(double temperature, double* average, int n_spins, int counter_accept, int m_carlos)
{
    double over_iter = 1/(double) m_carlos;
    double E = average[0] * over_iter;
    double EE = average[1] * over_iter;

    double M = average[2] * over_iter;
    double MM = average[3] * over_iter;
    double M_abs = average[4] * over_iter;

    double E_var = (EE - E*E )/n_spins/n_spins;
    double M_var = (MM - M_abs * M_abs)/n_spins/n_spins;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << E/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << E_var/temperature/temperature;
    // ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << M_var/temperature;
    ofile << setw(15) << setprecision(8) << M_abs/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << (double) counter_accept << endl;
}

int main(int argc, char *argv[])
{
    string output_filename;
    int n_spins, max_carlos;
    double w[17], average[5];
    double init_temp = 2.4;
    double final_temp = 2.4;
    double step_temp = 0.1;
    double E, M  = 0;
    string initial_type = "random";

    if (argc < 3){
        cout<< "wrong usage" << endl;
        exit(1);
    }
    else{
        n_spins = atoi(argv[1]);
        max_carlos = atoi(argv[2]);
        output_filename = string(argv[3]);
        ofile.open(output_filename);
    }

    arma::imat spin_matrix(n_spins, n_spins);
    spin_matrix.zeros();

    int counter_accept = 0;


    for(double temp = init_temp; temp <= final_temp ; temp += step_temp){
        E = 0;
        M = 0;

        for(int de = -8; de <= 8; de ++) w[de+8] = 0;
        for(int de = -8; de <= 8; de +=4) w[de+8] = exp(-de/temp);
        for(int i = 0; i<5 ; i++ ) average[i] =0;
        spin_matrix = initialize (spin_matrix,  n_spins,  E,  M, initial_type);

        int writer_controller = 0;
        int write_every_n = 100;
        srand (time(NULL));
        for(int i = 0 ; i<=max_carlos; i++)
        {
            spin_matrix = Metropolis(n_spins, spin_matrix, E, M, w, counter_accept);

            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);


            if(writer_controller == write_every_n){
                output(temp, average, n_spins, counter_accept, i+1);
                writer_controller = 0;
                cout <<  (i / (double) max_carlos ) * 100 << "%\r";
                cout.flush();
            }
            writer_controller ++;


        }
        //output(temp, average, n_spins, counter_accept, max_carlos);
    }
    ofile.close();
}
