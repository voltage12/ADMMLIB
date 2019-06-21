#include <cmath>

#include "optimizer/gd_optimizer.h"
#include "utils/math_utils.h"


GDOptimizer::GDOptimizer(Properties &properties, Function *function_)
        : function(function_) {
    optimizer_max_iter_num = properties.get_int("optimizer_max_iter_num");
    dim = properties.get_int("dim");
    gamma = properties.get_double("gamma");
    base_learning_rate = properties.get_double("base_learning_rate");
    end_learning_rate = properties.get_double("end_learning_rate");
    learning_rate = base_learning_rate;
    momentum = properties.get_double("momentum");
    optimizer_epsilon = properties.get_double("optimizer_epsilon");
}

void GDOptimizer::minimize(const double *y, const double *z, double *x) {
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
    double *v = NULL;
    if (momentum > 0) {
        v = new double[dim];
        for (int i = 0; i < dim; ++i) {
            v[i] = 0;
        }
    }
    while (k < optimizer_max_iter_num) {
        if (k % 10 == 0) {
            f_value_pre = function->function_value(x, y, z);
        }
        learning_rate_decay(k);
        if (momentum <= 0) { /* 使用最简单的梯度下降法 */
            for (int i = 0; i < dim; ++i) {
                x[i] -= learning_rate * g[i];
            }
        } else {
            for (int i = 0; i < dim; ++i) {
                v[i] = momentum * v[i] - learning_rate * g[i];
                x[i] += v[i];
            }
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
        ++k;
    }
    delete[] g;
    delete[] v;
}

void GDOptimizer::learning_rate_decay(int step) {
    learning_rate = end_learning_rate +
                    (base_learning_rate - end_learning_rate) * pow(1 - 1.0 * step / optimizer_max_iter_num, gamma);
}
