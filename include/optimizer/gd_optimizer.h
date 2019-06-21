#ifndef GD_OPTIMIZER_H
#define GD_OPTIMIZER_H

#include "optimizer/optimizer.h"
#include "optimizer/function.h"
#include "utils/properties.h"


class GDOptimizer: public Optimizer {
public:
    GDOptimizer(Properties &properties, Function *function_);

    void minimize(const double *y, const double *z, double *x);

private:
    int optimizer_max_iter_num;
    int dim;

    double gamma;
    double learning_rate;
    double base_learning_rate;
    double end_learning_rate;
    double momentum;
    double optimizer_epsilon;

    Function *function;

    void learning_rate_decay(int step);
};


#endif
