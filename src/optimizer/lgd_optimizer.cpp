#include <cmath>

#include "optimizer/lgd_optimizer.h"
#include "utils/math_utils.h"


LGDOptimizer::LGDOptimizer(Properties &properties, Function *function_)
        : function(function_) {
    optimizer_max_iter_num = properties.get_int("optimizer_max_iter_num");
    dim = properties.get_int("dim");
    optimizer_epsilon = properties.get_double("optimizer_epsilon");
}

void LGDOptimizer::minimize(const double *y, const double *z, double *x) {
    int k = 0;
    double g_norm, x_norm;
    double *g = new double[dim];
    function->gradient(x, y, z, g);
    g_norm = l2_norm(g, dim);
    x_norm = l2_norm(x, dim);
    x_norm = (x_norm > 1 ? x_norm : 1);
    if (l1_norm(g, dim) < optimizer_epsilon || g_norm < optimizer_epsilon * x_norm) {
        delete[] g;
        return;
    }

    double f_value_pre = function->function_value(x, y, z);
    double f_value;
    double *g_old = new double[dim];
    double *d = new double[dim];
    double *delta_x = new double[dim];
    double *delta_g = new double[dim];
    while (k < optimizer_max_iter_num) {
        for (int i = 0; i < dim; ++i) {
            d[i] = -g[i];
        }
        double step;
        if (k == 0) {
            step = 1;
        } else {
            step = dot_product(delta_x, delta_x, dim) / dot_product(delta_x, delta_g, dim);
            step = step > 0 ? step : 1;
        }
        line_search(x, y, z, g, d, dim, function, step);
        if (k % 10 == 0) {
            f_value_pre = function->function_value(x, y, z);
        }
        for (int i = 0; i < dim; ++i) {
            g_old[i] = g[i];
            delta_x[i] = d[i] * step;
            x[i] = x[i] + delta_x[i];
        }
        if (k % 10 == 0) {
            f_value = function->function_value(x, y, z);
            if (std::fabs((f_value - f_value_pre) / f_value_pre) * 10 < optimizer_epsilon) {
                break;
            }
        }
        function->gradient(x, y, z, g);
        g_norm = l2_norm(g, dim);
        x_norm = l2_norm(x, dim);
        x_norm = (x_norm > 1 ? x_norm : 1);
        if (l1_norm(g, dim) < optimizer_epsilon || g_norm < optimizer_epsilon * x_norm) {
            break;
        }
        for (int i = 0; i < dim; ++i) {
            delta_g[i] = g[i] - g_old[i];
        }
        ++k;
    }
    delete[] g;
    delete[] g_old;
    delete[] d;
    delete[] delta_g;
    delete[] delta_x;
}
