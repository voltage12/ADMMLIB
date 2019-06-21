#ifndef FUNCTION_H
#define FUNCTION_H


class Function {
public:
    virtual double function_value(const double *x, const double *y, const double *z) = 0;

    virtual void gradient(const double *x, const double *y, const double *z, double *g) = 0;

    virtual ~Function() = default;
};


#endif
