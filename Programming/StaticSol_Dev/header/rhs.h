#ifndef RHS_H
#define RHS_H


struct ode_params{

    double lambda;
    double dlambda;
    double *energies;
    int len_energies;
};

int rhs (double t, const double spin[], double f[], void * _params);
#endif