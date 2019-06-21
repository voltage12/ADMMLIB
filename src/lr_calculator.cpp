#include <string>
#include <sstream>
#include <iostream>

#include "lr_calculator.h"
#include "optimizer/lr_function.h"
#include "optimizer/tron_optimizer.h"
#include "optimizer/gd_optimizer.h"
#include "optimizer/lgd_optimizer.h"
#include "optimizer/lbfgs_optimizer.h"


LRCalculator::LRCalculator(int id_, Properties &properties) : id(id_) {
    dim = properties.get_int("dim");
    l2reg = properties.get_double("l2reg");
    std::string optimizer_name = properties.get_string("optimizer_name");
    std::string data_file_dir_path = properties.get_string("data_file_dir_path");
    std::stringstream string_builder;
    string_builder << data_file_dir_path << "data";
    if (id < 10) {
        string_builder << "00";
    } else if (id < 100) {
        string_builder << "0";
    }
    string_builder << id;
    if (optimizer_name == "gd") {
        function = new LRFunction(properties, string_builder.str());
        optimizer = new GDOptimizer(properties, function);
    } else if (optimizer_name == "tron") {
        optimizer = new TronOptimizer(properties, string_builder.str());
        function = NULL;
    } else if (optimizer_name == "lgd") {
        function = new LRFunction(properties, string_builder.str());
        optimizer = new LGDOptimizer(properties, function);
    } else if (optimizer_name == "lbfgs") {
        function = new LRFunction(properties, string_builder.str());
        optimizer = new LBFGSOptimizer(properties, function);
    }
}

void LRCalculator::update_x(double rho, const double *y, const double *z, double *x) {
    optimizer->minimize(y, z, x);
}

void LRCalculator::update_y(double rho, const double *x, const double *z, double *y) {
    for (int i = 0; i < dim; ++i) {
        y[i] += rho * (x[i] - z[i]);
    }
}

void LRCalculator::update_z(int worker_num, double rho, const double *rho_x_plus_y, double *z) {
    double temp = 2 * l2reg + worker_num * rho;
    for (int i = 0; i < dim; ++i) {
        z[i] = rho_x_plus_y[i] / temp;
    }
}

void LRCalculator::init(double *x, double *y, double *z) {
    /* 可以使用热启动，读取先前训练的结果，在错误恢复中有用 */
    for (int i = 0; i < dim; ++i) {
        x[i] = 0;
        y[i] = 0;
        z[i] = 0;
    }
}

LRCalculator::~LRCalculator() {
    delete optimizer;
    delete function;
}
