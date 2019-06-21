#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "optimizer/tron_optimizer.h"
#include "utils/math_utils.h"


int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, ix, iy, nn, iincx, iincy;
    register double ssa;

    /* constant times a vector plus a vector.
       uses unrolled loop for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    ssa = *sa;
    iincx = *incx;
    iincy = *incy;

    if (nn > 0 && ssa != 0.0) {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn - 3;
            for (i = 0; i < m; i += 4) {
                sy[i] += ssa * sx[i];
                sy[i + 1] += ssa * sx[i + 1];
                sy[i + 2] += ssa * sx[i + 2];
                sy[i + 3] += ssa * sx[i + 3];
            }
            for (; i < nn; ++i) /* clean-up loop */
                sy[i] += ssa * sx[i];
        } else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
            iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
            for (i = 0; i < nn; i++) {
                sy[iy] += ssa * sx[ix];
                ix += iincx;
                iy += iincy;
            }
        }
    }

    return 0;
} /* daxpy_ */

double ddot_(int *n, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, nn, iincx, iincy;
    double stemp;
    long int ix, iy;

    /* forms the dot product of two vectors.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    iincy = *incy;

    stemp = 0.0;
    if (nn > 0) {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn - 4;
            for (i = 0; i < m; i += 5)
                stemp += sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * sy[i + 2] +
                         sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];

            for (; i < nn; i++)        /* clean-up loop */
                stemp += sx[i] * sy[i];
        } else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = 0;
            iy = 0;
            if (iincx < 0)
                ix = (1 - nn) * iincx;
            if (iincy < 0)
                iy = (1 - nn) * iincy;
            for (i = 0; i < nn; i++) {
                stemp += sx[ix] * sy[iy];
                ix += iincx;
                iy += iincy;
            }
        }
    }

    return stemp;
} /* ddot_ */

double dnrm2_(int *n, double *x, int *incx) {
    long int ix, nn, iincx;
    double norm, scale, absxi, ssq, temp;

/*  DNRM2 returns the euclidean norm of a vector via the function
    name, so that

       DNRM2 := sqrt( x'*x )

    -- This version written on 25-October-1982.
       Modified on 14-October-1993 to inline the call to SLASSQ.
       Sven Hammarling, Nag Ltd.   */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;

    if (nn > 0 && iincx > 0) {
        if (nn == 1) {
            norm = fabs(x[0]);
        } else {
            scale = 0.0;
            ssq = 1.0;

            /* The following loop is equivalent to this call to the LAPACK
               auxiliary routine:   CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

            for (ix = (nn - 1) * iincx; ix >= 0; ix -= iincx) {
                if (x[ix] != 0.0) {
                    absxi = fabs(x[ix]);
                    if (scale < absxi) {
                        temp = scale / absxi;
                        ssq = ssq * (temp * temp) + 1.0;
                        scale = absxi;
                    } else {
                        temp = absxi / scale;
                        ssq += temp * temp;
                    }
                }
            }
            norm = scale * sqrt(ssq);
        }
    } else
        norm = 0.0;

    return norm;

} /* dnrm2_ */

int dscal_(int *n, double *sa, double *sx, int *incx) {
    long int i, m, nincx, nn, iincx;
    double ssa;

    /* scales a vector by a constant.
       uses unrolled loops for increment equal to 1.
       jack dongarra, linpack, 3/11/78.
       modified 3/93 to return if incx .le. 0.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    ssa = *sa;

    if (nn > 0 && iincx > 0) {
        if (iincx == 1) /* code for increment equal to 1 */
        {
            m = nn - 4;
            for (i = 0; i < m; i += 5) {
                sx[i] = ssa * sx[i];
                sx[i + 1] = ssa * sx[i + 1];
                sx[i + 2] = ssa * sx[i + 2];
                sx[i + 3] = ssa * sx[i + 3];
                sx[i + 4] = ssa * sx[i + 4];
            }
            for (; i < nn; ++i) /* clean-up loop */
                sx[i] = ssa * sx[i];
        } else /* code for increment not equal to 1 */
        {
            nincx = nn * iincx;
            for (i = 0; i < nincx; i += iincx)
                sx[i] = ssa * sx[i];
        }
    }

    return 0;
} /* dscal_ */

double uTMv(int n, double *u, double *M, double *v) {
    const int m = n - 4;
    double res = 0;
    int i;
    for (i = 0; i < m; i += 5)
        res += u[i] * M[i] * v[i] + u[i + 1] * M[i + 1] * v[i + 1] + u[i + 2] * M[i + 2] * v[i + 2] +
               u[i + 3] * M[i + 3] * v[i + 3] + u[i + 4] * M[i + 4] * v[i + 4];
    for (; i < n; i++)
        res += u[i] * M[i] * v[i];
    return res;
}

TronOptimizer::TronOptimizer(Properties &properties, const std::string &data_file_path) : sample_set(data_file_path) {
    dim = properties.get_int("dim");
    rho = properties.get_double("rho");
    optimizer_max_iter_num = properties.get_int("optimizer_max_iter_num");
    optimizer_epsilon = properties.get_double("optimizer_epsilon");
    optimizer_epsilon_cg = properties.get_double("optimizer_epsilon_cg");
    int sample_num = sample_set.get_sample_num();
    D = new double[sample_num];
    int pos = 0;
    int neg = 0;
    for (int i = 0; i < sample_num; i++)
        if (sample_set.get_label(i) > 0)
            ++pos;
    neg = sample_num - pos;
    optimizer_epsilon = optimizer_epsilon * std::max(std::min(pos, neg), 1) / sample_num;
}

TronOptimizer::~TronOptimizer() {
    delete[] D;
}

void TronOptimizer::minimize(const double *y, const double *z, double *x) {
    // Parameters for updating the iterates.
    double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;
    // Parameters for updating the trust region size delta.
    double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

    int i;
    double delta = 0, sMnorm, one = 1.0;
    double alpha, f, fnew, prered, actred, gs;
    int search = 1, iter = 1, inc = 1;
    double *s = new double[dim];
    double *r = new double[dim];
    double *g = new double[dim];

    const double alpha_pcg = 0.01;
    double *M = new double[dim];

    // calculate gradient norm at x=0 for stopping condition.
    double *x0 = new double[dim];
    for (i = 0; i < dim; i++) {
        x0[i] = 0;
    }
    function_value(x0, y, z);
    gradient(x0, y, z, g);
    double gnorm0 = dnrm2_(&dim, g, &inc);
    delete[] x0;

    f = function_value(x, y, z);
    gradient(x, y, z, g);
    double gnorm = dnrm2_(&dim, g, &inc);

    if (gnorm <= optimizer_epsilon * gnorm0)
        search = 0;

    get_diag_preconditioner(M);
    for (i = 0; i < dim; i++)
        M[i] = (1 - alpha_pcg) + alpha_pcg * M[i];
    delta = sqrt(uTMv(dim, g, M, g));

    double *x_new = new double[dim];
    bool reach_boundary;
    bool delta_adjusted = false;
    while (iter <= optimizer_max_iter_num && search) {
        trpcg(delta, g, M, s, r, &reach_boundary);

        memcpy(x_new, x, sizeof(double) * dim);
        daxpy_(&dim, &one, s, &inc, x_new, &inc);

        gs = ddot_(&dim, g, &inc, s, &inc);
        prered = -0.5 * (gs - ddot_(&dim, s, &inc, r, &inc));
        fnew = function_value(x_new, y, z);

        // Compute the actual reduction.
        actred = f - fnew;

        // On the first iteration, adjust the initial step bound.
        sMnorm = sqrt(uTMv(dim, s, M, s));
        if (iter == 1 && !delta_adjusted) {
            delta = std::min(delta, sMnorm);
            delta_adjusted = true;
        }

        // Compute prediction alpha*sMnorm of the step.
        if (fnew - f - gs <= 0)
            alpha = sigma3;
        else
            alpha = std::max(sigma1, -0.5 * (gs / (fnew - f - gs)));

        // Update the trust region bound according to the ratio of actual to predicted reduction.
        if (actred < eta0 * prered)
            delta = std::min(alpha * sMnorm, sigma2 * delta);
        else if (actred < eta1 * prered)
            delta = std::max(sigma1 * delta, std::min(alpha * sMnorm, sigma2 * delta));
        else if (actred < eta2 * prered)
            delta = std::max(sigma1 * delta, std::min(alpha * sMnorm, sigma3 * delta));
        else {
            if (reach_boundary)
                delta = sigma3 * delta;
            else
                delta = std::max(delta, std::min(alpha * sMnorm, sigma3 * delta));
        }

        if (actred > eta0 * prered) {
            ++iter;
            memcpy(x, x_new, sizeof(double) * dim);
            f = fnew;
            gradient(x, y, z, g);
            get_diag_preconditioner(M);
            for (i = 0; i < dim; i++)
                M[i] = (1 - alpha_pcg) + alpha_pcg * M[i];

            gnorm = dnrm2_(&dim, g, &inc);
            if (gnorm <= optimizer_epsilon * gnorm0)
                break;
        }
        if (f < -1.0e+32) {
            std::cerr << "WARNING: f < -1.0e+32" << std::endl;
            break;
        }
        if (prered <= 0) {
            std::cerr << "WARNING: prered <= 0" << std::endl;
            break;
        }
        if (fabs(actred) <= 1.0e-12 * fabs(f) && fabs(prered) <= 1.0e-12 * fabs(f)) {
            std::cerr << "WARNING: actred and prered too small" << std::endl;
            break;
        }
    }

    delete[] g;
    delete[] r;
    delete[] x_new;
    delete[] s;
    delete[] M;
}

int TronOptimizer::trpcg(double delta, double *g, double *M, double *s, double *r, bool *reach_boundary) {
    int i, inc = 1;
    double one = 1;
    double *d = new double[dim];
    double *Hd = new double[dim];
    double zTr, znewTrnew, alpha, beta, cgtol;
    double *z = new double[dim];

    *reach_boundary = false;
    for (i = 0; i < dim; i++) {
        s[i] = 0;
        r[i] = -g[i];
        z[i] = r[i] / M[i];
        d[i] = z[i];
    }

    zTr = ddot_(&dim, z, &inc, r, &inc);
    cgtol = optimizer_epsilon_cg * sqrt(zTr);
    int cg_iter = 0;
    int max_cg_iter = std::max(dim, 5);

    while (cg_iter < max_cg_iter) {
        if (sqrt(zTr) <= cgtol)
            break;
        cg_iter++;
        Hv(d, Hd);

        alpha = zTr / ddot_(&dim, d, &inc, Hd, &inc);
        daxpy_(&dim, &alpha, d, &inc, s, &inc);

        double sMnorm = sqrt(uTMv(dim, s, M, s));
        if (sMnorm > delta) {
            //std::cerr << "cg reaches trust region boundary" << std::endl;
            *reach_boundary = true;
            alpha = -alpha;
            daxpy_(&dim, &alpha, d, &inc, s, &inc);

            double sTMd = uTMv(dim, s, M, d);
            double sTMs = uTMv(dim, s, M, s);
            double dTMd = uTMv(dim, d, M, d);
            double dsq = delta * delta;
            double rad = sqrt(sTMd * sTMd + dTMd * (dsq - sTMs));
            if (sTMd >= 0)
                alpha = (dsq - sTMs) / (sTMd + rad);
            else
                alpha = (rad - sTMd) / dTMd;
            daxpy_(&dim, &alpha, d, &inc, s, &inc);
            alpha = -alpha;
            daxpy_(&dim, &alpha, Hd, &inc, r, &inc);
            break;
        }
        alpha = -alpha;
        daxpy_(&dim, &alpha, Hd, &inc, r, &inc);

        for (i = 0; i < dim; i++)
            z[i] = r[i] / M[i];
        znewTrnew = ddot_(&dim, z, &inc, r, &inc);
        beta = znewTrnew / zTr;
        dscal_(&dim, &beta, d, &inc);
        daxpy_(&dim, &one, z, &inc, d, &inc);
        zTr = znewTrnew;
    }

    if (cg_iter == max_cg_iter)
        std::cerr << "WARNING: reaching maximal number of CG steps" << std::endl;

    delete[] d;
    delete[] Hd;
    delete[] z;

    return cg_iter;
}

double TronOptimizer::function_value(const double *x, const double *y, const double *z) {
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

void TronOptimizer::gradient(const double *x, const double *y, const double *z, double *g) {
    for (int i = 0; i < dim; ++i) {
        g[i] = y[i] + rho * (x[i] - z[i]);
    }
    int sample_num = sample_set.get_sample_num();
    for (int i = 0; i < sample_num; ++i) {
        double temp = sigmoid(sample_set.get_label(i) * sample_set.dot(i, x));
        D[i] = temp * (1 - temp);
        temp = (temp - 1) * sample_set.get_label(i);
        const Feature *sample = sample_set.get_sample(i);
        while (sample->index != -1) {
            g[(sample->index - 1)] += (sample->value * temp);
            ++sample;
        }
    }
}

void TronOptimizer::get_diag_preconditioner(double *M) {
    int sample_num = sample_set.get_sample_num();
    for (int i = 0; i < dim; i++) {
        M[i] = 1;
    }
    for (int i = 0; i < sample_num; i++) {
        const Feature *sample = sample_set.get_sample(i);
        while (sample->index != -1) {
            M[sample->index - 1] += sample->value * sample->value * D[i];
            ++sample;
        }
    }
}

void TronOptimizer::Hv(double *s, double *Hs) {
    int sample_num = sample_set.get_sample_num();
    for (int i = 0; i < dim; i++) {
        Hs[i] = 0;
    }
    for (int i = 0; i < sample_num; i++) {
        const Feature *sample = sample_set.get_sample(i);
        double xTs = sample_set.dot(i, s);
        xTs = D[i] * xTs;
        while (sample->index != -1) {
            Hs[sample->index - 1] += sample->value * xTs;
            ++sample;
        }
    }
    for (int i = 0; i < dim; i++)
        Hs[i] += rho * s[i];
}
