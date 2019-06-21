#ifndef LR_FUNCTION_H
#define LR_FUNCTION_H

#include "optimizer/function.h"
#include "optimizer/sparse_sample_set.h"
#include "utils/properties.h"


/* 用ADMM解决Logistics Regression问题 */
class LRFunction: public Function {
public:
    LRFunction(Properties &properties, const std::string &data_file_path);

    double function_value(const double *x, const double *y, const double *z);

    void gradient(const double *x, const double *y, const double *z, double *g);

private:
    int dim;

    double rho;

    SparseSampleSet sample_set;
};


#endif
