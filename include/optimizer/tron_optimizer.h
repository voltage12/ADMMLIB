#ifndef TRONOPTIMIZER_H
#define TRONOPTIMIZER_H

#include "optimizer/optimizer.h"
#include "optimizer/sparse_sample_set.h"
#include "utils/properties.h"


class TronOptimizer : public Optimizer {
public:
    TronOptimizer(Properties &properties, const std::string &data_file_path);

    ~TronOptimizer();

    void minimize(const double *y, const double *z, double *x);

private:
    int dim;
    int optimizer_max_iter_num;

    double optimizer_epsilon;
    double optimizer_epsilon_cg;
    double rho;
    double *D;

    SparseSampleSet sample_set;

    int trpcg(double delta, double *g, double *M, double *s, double *r, bool *reach_boundary);

    double function_value(const double *x, const double *y, const double *z);

    void gradient(const double *x, const double *y, const double *z, double *g);

    void get_diag_preconditioner(double *M);

    void Hv(double *s, double *Hs);
};



#endif
