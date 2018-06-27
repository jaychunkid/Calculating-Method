#include"data_fitting.h"

void LeastSquares(long double &a, long double &b, unsigned size,
                  long double *x, long double *y) {
    long double s1 = 0;
    long double s2 = 0;
    long double s3 = 0;
    long double s4 = 0;
    for (unsigned i = 0; i < size; ++i) {
        s1 += x[i];
        s2 += y[i];
        s3 += x[i] * x[i];
        s4 += x[i] * y[i];
    }
    unsigned n = size - 1;
    a = (s2*s3 - s1 * s4) / (n*s3 - s1 * s1);
    b = (n*s4 - s1 * s2) / (n*s3 - s1 * s1);
    return;
}