#include <cmath>
#include <iostream>

#include "utils/math_utils.h"


double dot_product(const double *a, const double *b, int dim) {
    double temp = 0.0;
    for (int i = 0; i < dim; ++i) {
        temp += a[i] * b[i];
    }
    return temp;
}

double l1_norm(const double *a, int dim) {
    double dmax = fabs(a[0]);
    for (int i = 1; i < dim; ++i) {
        if (fabs(a[i]) > dmax) {
            dmax = fabs(a[i]);
        }
    }
    return dmax;
}

double l2_norm(const double *a, int dim) {
    return std::sqrt(dot_product(a, a, dim));
}

double sigmoid(double x) {
    return 1.0 / (1 + exp(-x));
}

void soft_threshold(double t, int dim, double *z) {
    for (int i = 0; i < dim; ++i) {
        if (z[i] > t)
            z[i] -= t;
        else if (z[i] <= t && z[i] >= -t) {
            z[i] = 0.0;
        } else {
            z[i] += t;
        }
    }
}

bool line_search(const double *x, const double *y, const double *z, const double *g, const double *d, int dim,
                 Function *function, double &step, int max_search_time, double r, double c) {
    int k = 0;
    double func_val, func_val_next;
    double dginit = dot_product(g, d, dim);
    if (step < 0 || dginit > 0) {
        std::cerr << "step < 0 || dginit > 0" << std::endl;
        exit(-1);
    }

    func_val = function->function_value(x, y, z);
    double *x_new = new double[dim];
    while (k < max_search_time) {
        for (int i = 0; i < dim; ++i) {
            x_new[i] = x[i] + step * d[i];
        }
        func_val_next = function->function_value(x_new, y, z);
        if (func_val_next <= func_val + c * step * dginit) {
            delete[] x_new;
            return true;
        } else {
            step *= r;
            if (!std::isfinite(step)) {
                std::cerr << "line search failed: step become inf" << std::endl;
                exit(-1);
            }
        }
        ++k;
    }
    std::cout << "line search failed: max search time" << std::endl;
    delete[] x_new;
    return false;
}
