#ifndef CALCULATOR_H
#define CALCULATOR_H


class Calculator {
public:
    virtual void update_x(double rho, const double *y, const double *z, double *x) = 0;

    virtual void update_y(double rho, const double *x, const double *z, double *y) = 0;

    virtual void update_z(int worker_num, double rho, const double *rho_x_plus_y, double *z) = 0;

    virtual void init(double *x, double *y, double *z) = 0;

    virtual ~Calculator() = default;
};


#endif //MCCADMM_CALCULATOR_H
