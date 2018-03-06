#include <iostream>
#include <fstream>
#include <cmath>

#include "../header/rhs.h"


using namespace std;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;


int main(int ac, char** av){

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
                    ("help,h", "Print help messages")
                    ("lambda,l", po::value<double>(), "g BCS interaction strength (default=1)")
                    ("dlambda,d", po::value<double>(), "dg quench of the BCS interaction strength (default=0)")
                    ("tf", po::value<double>(), "tf finite time (default = 10)")
                    ("omega_D,w", po::value<double>(), "omega_D Debye-frequency")
                    ("Delta_0,D",po::value<double>(), "Delta_0 BCS-orderparameter")
                    ("file_name,f", po::value<string>(), "filename");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm); 

    ode_params params;
    string fname = "data.txt";
    fstream file;

    //Physikalisch relevante Größen und Simulationsgrößen setzen
    params.lambda = 0.3147;
    params.dlambda = 0.0;
    double Delta_0 =2.6;
    double omega_Debye =12.0 * Delta_0;
    params.len_energies = 3;
    params.energies = new double[params.len_energies];

    double t_i = 0.0, t_f = 10.0;
    int N_times = 10;
    double t_step;

    //Abfragen ob irgendwelche parameter übergeben worden sind und setze diese

    if(vm.count("help"))
    {
        std::cout << "Basic Command Line Parameter App" << std::endl 
                  << desc << std::endl;

        return 0;
    }

    if(vm.count("lambda"))
        params.lambda = vm["lambda"].as<double>();
 
    if(vm.count("dlambda"))
        params.dlambda = vm["dlambda"].as<double>();

    if(vm.count("tf"))
        t_f = vm["tf"].as<double>();

    if (vm.count("omega_D"))
        omega_Debye = vm["omega_D"].as<double>();

    if(vm.count("Delta_0"))
        Delta_0 = vm["Delta_0"].as<double>();

    if(vm.count("file_name"))
        fname = vm["file_name"].as<string>();
        

   
    //Diskretisierung der Energie-Shell
    for (int i = 0; i < params.len_energies; i++){

        params.energies[i] = omega_Debye * (-1.0 + 2.0*i/(params.len_energies -1));
    }

    //Anfangswert setzen
    double spins[3 * params.len_energies];

    for (int i = 0; i < params.len_energies; i++){

        spins[3 * i + 0] = Delta_0; ///sqrt(pow(Delta_0, 2) + pow(params.energies[i],2));
        spins[3 * i + 1] = 0.0;
        spins[3 * i + 2] = -params.energies[i]; ///sqrt(pow(Delta_0, 2) + pow(params.energies[i],2));
    }
    //############## SIMULATION ####################

    //step1 : Definiere das ODE-System
    gsl_odeiv2_system sys = {rhs, NULL, 3 * params.len_energies, &params};
    //step2 : Definiere den Solver
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6,0);


    //Dateiverwaltung
    file.open(fname.c_str(), ios::out);
    file << "#lambda \t dlambda \t omega_D \t Delta_0" << endl;
    file << "#" << "\t" << params.lambda << "\t" << params.dlambda << "\t" << omega_Debye << "\t" << Delta_0 << endl;
    file << "#t \t Delta_x \t Delta_y" << endl;


    //Iteration des Systems
    for (int i = 0; i <= N_times; i++){

        t_step = i * t_f/N_times;

        int status = gsl_odeiv2_driver_apply(driver, &t_i, t_step, spins);

        if(status != GSL_SUCCESS){

            cout << "error, return value=" << status << endl;
            break;
        }

    double Delta_x = 0.0, Delta_y = 0.0;
    for (int i = 0; i < params.len_energies; i++){

        Delta_x += spins[i * 3 + 0];
        Delta_y += spins[i * 3 + 1];
    }

    Delta_x *= (params.lambda + params.dlambda)/(2.0 * params.len_energies);
    Delta_y *= (params.lambda + params.dlambda)/(2.0 * params.len_energies);

    file << t_step << "\t" << Delta_x << "\t" << Delta_y << endl;
    }

    //clean up
    file.close();
    gsl_odeiv2_driver_free(driver);
    delete params.energies;
    return 0;
}