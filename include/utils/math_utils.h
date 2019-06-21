#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "optimizer/function.h"


double dot_product(const double *a, const double *b, int dim);

double l1_norm(const double *a, int dim);

double l2_norm(const double *a, int dim);

double sigmoid(double x);

void soft_threshold(double t, int dim, double *z);

bool line_search(const double *x, const double *y, const double *z, const double *g, const double *d, int dim,
                 Function *function, double &step, int max_search_time = 30, double r = 0.5, double c = 0.1);

#endif
