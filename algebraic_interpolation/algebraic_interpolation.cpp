#include<vector>

#include"algebraic_interpolation.h"

using std::vector;

long double Lagrange(unsigned size, const long double *x, const long double *y, long double xx) {
    long double yy = 0.0;
    for (unsigned i = 0; i < size; ++i) {
        long double t = 1.0;
        for (unsigned k = 0; k < size; ++k) {
            if (k != i) {
                t *= (xx - x[k]) / (x[i] - x[k]);
            }
        }
        yy += t * y[i];
    }
    return yy;
}


long double Newton(unsigned size, const long double *x, const long double *y, long double xx){
	long double *V = new long double[size];
    long double w = 1.0;
    long double yy = y[0];
	V[0] = y[0];
	for(unsigned k = 1; k < size; ++k){
		V[k] = y[k];
		for(unsigned i = 0; i < k; ++i){
			V[k] = (V[i] - V[k]) / (x[i] - x[k]);
		}
		w *= xx - x[k - 1];
        long double tmp = w * V[k];
		yy += tmp;
	}
	delete[] V;
	return yy;
}


long double CubicSpline(unsigned size, const long double *x, const long double *y, long double xx){
	long double *h = new long double[size - 1];
	for(unsigned i = 0; i < size - 1; ++i){
		h[i] = x[i + 1] - x[i];
	}
	long double *alpha = new long double[size];
	long double *beta = new long double[size];
	alpha[0] = 1;
	alpha[size - 1] = 0;
	beta[0] = 3 / h[0] * (y[1] - y[0]);
	beta[size - 1] = 3 / h[size - 2] * (y[size - 1] - y[size - 2]);
	for(unsigned i = 1; i < size - 1; ++i){
		alpha[i] = h[i - 1] / (h[i -1] + h[i]);
		beta[i] = 3 * ((1 - alpha[i]) / h[i - 1] * (y[i] - y[i - 1]) + alpha[i] / h[i] * (y[i + 1] - y[i]));
	}
	long double *a = new long double[size];
	long double *b = new long double[size];
	a[0] = -alpha[0] / 2;
	b[0] = beta[0] / 2;
	for(unsigned i = 1; i < size; ++i){
		a[i] = -alpha[i] / (2 + (1 - alpha[i]) * a[i -1]);
		b[i] = (beta[i] - (1 - alpha[i]) * b[i - 1]) / (2 + (1 - alpha[i]) * a[i - 1]);
	}
	long double *m = new long double[size + 1];
	m[size] = 0;
	for(unsigned i = size; i > 0; --i){
		m[i - 1] = a[i - 1] * m[i] + b[i - 1];
	}
	for(unsigned i = 1; i < size; ++i){
		if(xx <= x[i]){
			long double tmp1 = (xx - x[i - 1]) / (x[i] - x[i - 1]);
			long double tmp2 = (xx - x[i]) / (x[i - 1] - x[i]);
			long double yy = (1 + 2 * tmp1) * tmp2 * tmp2 * y[i - 1] + 
				(1 + 2 * tmp2) * tmp1 * tmp1 * y[i] + 
				(xx - x[i - 1]) * tmp2 * tmp2 * m[i - 1] + 
				(xx - x[i]) * tmp1 * tmp1 * m[i];
			return yy;
		}
	}
	delete[] alpha;
	delete[] beta;
	delete[] a;
	delete[] b;
	delete[] m;
	return 0.0;
}


long double CubicSpline(unsigned size, const long double *x, const long double *y, long double m0, long double mn, long double xx){
	long double *h = new long double[size - 1];
	for(unsigned i = 0; i < size - 1; ++i){
		h[i] = x[i + 1] - x[i];
	}
	long double *alpha = new long double[size];
	long double *beta = new long double[size];
	alpha[0] = 0;
	alpha[size - 1] = 1;
	beta[0] = 2 * m0;
	beta[size - 1] = 2 * mn;
	for(unsigned i = 1; i < size - 1; ++i){
		alpha[i] = h[i - 1] / (h[i -1] + h[i]);
		beta[i] = 3 * ((1 - alpha[i]) / h[i - 1] * (y[i] - y[i - 1]) + alpha[i] / h[i] * (y[i + 1] - y[i]));
	}
	long double *a = new long double[size];
	long double *b = new long double[size];
	a[0] = -alpha[0] / 2;
	b[0] = beta[0] / 2;
	for(unsigned i = 1; i < size; ++i){
		a[i] = -alpha[i] / (2 + (1 - alpha[i]) * a[i -1]);
		b[i] = (beta[i] - (1 - alpha[i]) * b[i - 1]) / (2 + (1 - alpha[i]) * a[i - 1]);
	}
	long double *m = new long double[size + 1];
	m[size] = 0;
	for(unsigned i = size; i > 0; --i){
		m[i - 1] = a[i - 1] * m[i] + b[i - 1];
	}
	for(unsigned i = 1; i < size; ++i){
		if(xx <= x[i]){
			long double tmp1 = (xx - x[i - 1]) / (x[i] - x[i - 1]);
			long double tmp2 = (xx - x[i]) / (x[i - 1] - x[i]);
			long double yy = (1 + 2 * tmp1) * tmp2 * tmp2 * y[i - 1] + 
				(1 + 2 * tmp2) * tmp1 * tmp1 * y[i] + 
				(xx - x[i - 1]) * tmp2 * tmp2 * m[i - 1] + 
				(xx - x[i]) * tmp1 * tmp1 * m[i];
			return yy;
		}
	}
	delete[] alpha;
	delete[] beta;
	delete[] a;
	delete[] b;
	delete[] m;
	return 0.0;
}


