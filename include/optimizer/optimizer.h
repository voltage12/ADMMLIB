#ifndef OPTIMIZER_H
#define OPTIMIZER_H


class Optimizer {
public:
    virtual void minimize(const double *y, const double *z, double *x) = 0;

    virtual ~Optimizer() = default;
};


#endif
