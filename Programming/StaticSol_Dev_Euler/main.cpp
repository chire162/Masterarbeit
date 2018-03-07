#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;


int main(){

    ofstream file("data.txt");
    double omega_D = 1.0;
    int len_energies = 1000;
    double lambda = 0.4;
    double Delta_x = 0.0, Delta_y = 0.0;
    double Delta_0 = omega_D/sinh(1/lambda);
    double dt = 1e-2;

    VectorXd energies(2*len_energies);

    for (int i = 0; i < len_energies; i++){

            energies(i) = -omega_D * (len_energies - i)/len_energies;
            energies(i + len_energies) = -energies(i); 
    }

    MatrixXd Spins(2 * len_energies, 3);

    for(int i = 0; i < 2 * len_energies; i++){

        Spins(i,0) = Delta_0 / sqrt(pow(Delta_0, 2) + pow(energies(i),2));
        Spins(i,1) = 0.0;
        Spins(i,2) = -energies(i) / sqrt(pow(Delta_0, 2) + pow(energies(i),2));
    }

    for (int i = 0; i < 2* len_energies; i++){

        Delta_x += Spins(i,0);
        Delta_y += Spins(i,1);
    }

    Delta_x *= lambda / 2.0 /len_energies;
    Delta_y *= lambda / 2.0/ len_energies;


    for (int i = 0; i <= 1e5; i++){

        for (int k = 0; k < 2 * len_energies; k++){

            Spins(k,0) += 2.0 * (- Delta_y * Spins(k, 2) - energies(k) * Spins(k,1)) * dt;
            Spins(k,1) += 2.0 * (energies(k) * Spins(k,0) + Delta_x * Spins(k,2)) * dt;
            Spins(k,2) += 2.0 * (-Delta_x * Spins(k,1) + Delta_y * Spins(k,0)) * dt;
        }

        for (int j = 0; j < 2 * len_energies; j++){

            Delta_x += Spins(j,0);
            Delta_y += Spins(j,1);
        }

        Delta_x *= lambda * 0.8 / 2.0/ len_energies;
        Delta_y *= lambda * 0.8 / 2.0/len_energies;

        file << i * dt << "\t" << Delta_x << "\t" << Delta_y << endl;

        if( i % 10000 == 0)
            cout << i << endl;
    }

    return 0;
}