#include <cmath>

#include "utils/math_utils.h"
#include "optimizer/lbfgs_optimizer.h"


LBFGSOptimizer::LBFGSOptimizer(Properties &properties, Function *function_)
        : function(function_) {
    m = properties.get_int("m");
    dim = properties.get_int("dim");
    optimizer_max_iter_num = properties.get_int("optimizer_max_iter_num");
    optimizer_epsilon = properties.get_double("optimizer_epsilon");
}

void LBFGSOptimizer::minimize(const double *y, const double *z, double *x) {
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
    double beta;
    double *alpha = new double[m];
    double **delta_x = new double *[m];
    double **delta_g = new double *[m];
    for (int i = 0; i < m; ++i) {
        delta_x[i] = new double[dim];
        delta_g[i] = new double[dim];
    }
    double *g_old = new double[dim];
    double *d = new double[dim];
    while (k < optimizer_max_iter_num) {
        for (int i = 0; i < dim; ++i) {
            d[i] = -g[i];
        }
        int index;
        for (int i = k - 1; i >= 0 && i >= k - m; --i) {
            index = i % m;
            alpha[index] = dot_product(delta_x[index], d, dim) / dot_product(delta_x[index], delta_g[index], dim);
            for (int j = 0; j < dim; ++j) {
                d[j] -= alpha[index] * delta_g[index][j];
            }
        }
        if (k >= 1) {
            index = (k - 1) % m;
            double temp = dot_product(delta_x[index], delta_g[index], dim) /
                          dot_product(delta_g[index], delta_g[index], dim);
            for (int j = 0; j < dim; ++j) {
                d[j] *= temp;
            }
        }
        for (int i = (k - m >= 0 ? k - m : 0); i <= k - 1; ++i) {
            index = i % m;
            beta = dot_product(delta_g[index], d, dim) / dot_product(delta_x[index], delta_g[index], dim);
            for (int j = 0; j < dim; ++j) {
                d[j] += (alpha[index] - beta) * delta_x[index][j];
            }
        }
        double step = 1;
        line_search(x, y, z, g, d, dim, function, step);
        index = k % m;
        if (k % 10 == 0) {
            f_value_pre = function->function_value(x, y, z);
        }
        for (int j = 0; j < dim; ++j) {
            g_old[j] = g[j];
            delta_x[index][j] = d[j] * step;
            x[j] += delta_x[index][j];
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
        for (int j = 0; j < dim; ++j) {
            delta_g[index][j] = g[j] - g_old[j];
        }
        ++k;
    }
    for (int i = 0; i < m; ++i) {
        delete[] delta_g[i];
        delete[] delta_x[i];
    }
    delete[] g;
    delete[] g_old;
    delete[] d;
    delete[] delta_g;
    delete[] delta_x;
    delete[] alpha;
}
