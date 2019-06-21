#include <cmath>

#include "optimizer/lr_function.h"
#include "utils/math_utils.h"


LRFunction::LRFunction(Properties &properties, const std::string &data_file_path)
        : sample_set(data_file_path) {
    dim = properties.get_int("dim");
    rho = properties.get_double("rho");
}

double LRFunction::function_value(const double *x, const double *y, const double *z) {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i) {
        double temp = x[i] - z[i] + y[i] / rho;
        sum += temp * temp;
    }
    sum *= (rho / 2.0);
    int sample_num = sample_set.get_sample_num();
    for (int i = 0; i < sample_num; ++i) {
        sum += (std::log(1 + exp(-sample_set.get_label(i) * sample_set.dot(i, x))));
    }
    return sum;
}

void LRFunction::gradient(const double *x, const double *y, const double *z, double *g) {
    for (int i = 0; i < dim; ++i) {
        g[i] = y[i] + rho * (x[i] - z[i]);
    }
    int sample_num = sample_set.get_sample_num();
    for (int i = 0; i < sample_num; ++i) {
        double d = (sigmoid(sample_set.get_label(i) * sample_set.dot(i, x)) - 1) * sample_set.get_label(i);
        const Feature *sample = sample_set.get_sample(i);
        while (sample->index != -1) {
            g[(sample->index - 1)] += (sample->value * d);
            ++sample;
        }
    }
}

