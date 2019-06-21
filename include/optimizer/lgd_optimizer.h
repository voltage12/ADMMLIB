#ifndef LGD_OPTIMIZER_H
#define LGD_OPTIMIZER_H

#include "optimizer/optimizer.h"
#include "optimizer/function.h"
#include "utils/properties.h"


class LGDOptimizer: public Optimizer {
public:
    LGDOptimizer(Properties &properties, Function *function_);

    void minimize(const double *y, const double *z, double *x);

private:
    int optimizer_max_iter_num;
    int dim;

    double optimizer_epsilon;

    Function *function;
};


#endif
