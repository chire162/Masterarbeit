#include "../header/rhs.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
using namespace std;

int rhs (double t, const double spin[], double f[], void * _params){

    ode_params params = *(ode_params*)_params;
    double Delta_x = 0.0, Delta_y = 0.0;

    for (int i = 0; i < params.len_energies; i++){

        Delta_x += spin[3 * i + 0];
        Delta_y += spin[3 * i + 1];
    }

    Delta_x *= (params.lambda + params.dlambda)/(2.0 * params.len_energies);
    Delta_y *= (params.lambda + params.dlambda)/(2.0 * params.len_energies);

    for (int i = 0; i < params.len_energies; i++){

        f[3 * i + 0] = 2.0 * ( -Delta_y * spin[3 * i + 2] -  params.energies[i]  * spin[3 * i + 1]);
        f[3 * i + 1] = 2.0 * (params.energies[i] * spin[3 * i + 0] + Delta_x * spin[3 * i +2]);
        f[3 * i + 2] = 2.0 * (-Delta_x * spin[3 * i + 1] + Delta_y * spin[3 * i + 0]);

        cout << i << endl;
        cout << f[3 * i + 0] << endl;
        cout << f[3 * i + 1] << endl;
        cout << f[3 * i + 2] << endl;
        cout << "_____________" << endl;
    }
    return GSL_SUCCESS;
}