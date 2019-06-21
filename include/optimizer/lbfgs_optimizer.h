#ifndef LBFGS_OPTIMIZER_H
#define LBFGS_OPTIMIZER_H

#include "optimizer/optimizer.h"
#include "optimizer/function.h"
#include "utils/properties.h"


class LBFGSOptimizer: public Optimizer {
public:
    LBFGSOptimizer(Properties &properties, Function *function_);

    void minimize(const double *y, const double *z, double *x);

private:
    int optimizer_max_iter_num;
    int dim;
    int m;

    double optimizer_epsilon;

    Function *function;
};


#endif
