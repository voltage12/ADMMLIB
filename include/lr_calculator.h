#ifndef LR_CALCULATOR_H
#define LR_CALCULATOR_H

#include "calculator.h"
#include "utils/properties.h"
#include "optimizer/optimizer.h"
#include "optimizer/function.h"


class LRCalculator: public Calculator {
public:
    LRCalculator(int id_, Properties &properties);

    void update_x(double rho, const double *y, const double *z, double *x);

    void update_y(double rho, const double *x, const double *z, double *y);

    void update_z(int worker_num, double rho, const double *rho_x_plus_y, double *z);

    void init(double *x, double *y, double *z);

    ~LRCalculator();

private:
    int id;
    int dim;

    double l2reg;

    Optimizer *optimizer;
    Function *function;
};


#endif //MCCADMM_LRCALCULATOR_H
