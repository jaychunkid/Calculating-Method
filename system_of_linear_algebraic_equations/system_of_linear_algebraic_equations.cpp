#include<cmath>

#include"system_of_linear_algebraic_equations.h"

bool Gauss1(double** A, double* b, double e, unsigned n, double* x) {
    for (unsigned k = 0; k < n - 1; ++k) {
        if (std::abs(A[k][k]) <= e) {
            return false;
        }
        for (unsigned i = k + 1; i < n; ++i) {
            double T = A[i][k] / A[k][k];
            b[i] -= T * b[k];
            for (unsigned j = k + 1; j < n; ++j) {
                A[i][j] -= T * A[k][j];
            }
        }
    }
    if (std::abs(A[n - 1][n - 1]) <= e) {
        return false;
    }
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (unsigned i = n - 1; i > 0; --i) {
        double sum = A[i - 1][i] * x[i];
        for (unsigned j = i + 1; j < n; ++j) {
            sum += A[i - 1][j] * x[j];
        }
        x[i - 1] = (b[i - 1] - sum) / A[i - 1][i - 1];
    }
    return true;
}


bool Gauss2(double** A, double* b, double e, unsigned n, double* x) {
    for (unsigned k = 0; k < n - 1; ++k) {
        unsigned mark = k;
        for (unsigned i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > std::abs(A[mark][k])) {
                mark = i;
            }
        }
        if (std::abs(A[mark][k]) < e) {
            return false;
        }
        if (mark != k) {
            for (unsigned i = 0; i < n; ++i) {
                double tmp = A[k][i];
                A[k][i] = A[mark][i];
                A[mark][i] = tmp;
            }
            double tmp = b[k];
            b[k] = b[mark];
            b[mark] = tmp;
        }
        for (unsigned i = k + 1; i < n; ++i) {
            double T = A[i][k] / A[k][k];
            b[i] -= T * b[k];
            for (unsigned j = k + 1; j < n; ++j) {
                A[i][j] -= T * A[k][j];
            }
        }
    }
    if (std::abs(A[n - 1][n - 1]) <= e) {
        return false;
    }
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (unsigned i = n - 1; i > 0; --i) {
        double sum = A[i - 1][i] * x[i];
        for (unsigned j = i + 1; j < n; ++j) {
            sum += A[i - 1][j] * x[j];
        }
        x[i - 1] = (b[i - 1] - sum) / A[i - 1][i - 1];
    }
    return true;
}


bool Gauss3(double** A, double* b, double e, unsigned n, double* x){
    unsigned* d = new unsigned[n];
    for(unsigned i = 0; i < n; ++i){
        d[i] = i;
    }
    for(unsigned k = 0; k < n - 1; ++k){
        unsigned iMark = k;
        unsigned jMark = k;
        for(unsigned i = k + 1; i < n; ++i){
            for(unsigned j = k + 1; j < n;++j){
                if(A[i][j] > A[iMark][jMark]){
                    iMark = i;
                    jMark = j;
                }
            }
        }
        if(std::abs(A[iMark][jMark]) < e){
            delete[] d;
            return false;
        }
        if(iMark != k){
            for(unsigned i = 0; i < n; ++i){
                double tmp = A[k][i];
                A[k][i] = A[iMark][i];
                A[iMark][i] = tmp;
            }
            double tmp = b[k];
            b[k] = b[iMark];
            b[iMark] = tmp;
        }
        if(jMark != k){
            for(unsigned i = 0; i < n; ++i){
                double tmp = A[i][k];
                A[i][k] = A[i][jMark];
                A[i][jMark] = tmp;
            }
            unsigned tmp = d[jMark];
            d[jMark] = d[k];
            d[k] = tmp;
        }
        for(unsigned i = k + 1; i < n; ++i){
            double T = A[i][k] / A[k][k];
            b[i] -= T * b[k];
            for(unsigned j = k + 1; j < n; ++j){
                A[i][j] -= T * A[k][j];
            }
        }
    }
    double* Z = new double[n];
    Z[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for(unsigned i = n - 1; i > 0; --i){
        double sum = A[i - 1][i] * Z[i];
        for(unsigned j = i + 1; j < n; ++j){
            sum += A[i - 1][j] * Z[j];
        }
        Z[i - 1] = (b[i - 1] - sum) / A[i - 1][i - 1];
    }
    for(unsigned j = 0; j < n; ++j){
        x[d[j]] = Z[j];
    }
    delete[] d;
    delete[] Z;
    return true;
}

double** generateMatrix(unsigned n){
    double** M = new double*[n];
    for(unsigned i = 0; i < n; ++i){
        M[i] = new double[n];
    }
    return M;
}

void deposeMatrix(double** M, unsigned n){
    for(unsigned i = 0; i < n; ++i){
        delete[] M[i];
    }
    delete[] M;
}


bool LU(double** A, double* b, double e, unsigned n, double* x){
    double** u = generateMatrix(n);
    double** l = generateMatrix(n);
    for(unsigned k = 0; k < n; ++k){
        for(unsigned j = k; j < n; ++j){
            u[k][j] = A[k][j];
            for(unsigned m = 0; m < k; ++m){
                u[k][j] -= l[k][m] * u[m][j];
            }
            if(std::abs(u[k][k]) < e){
                deposeMatrix(u, n);
                deposeMatrix(l, n);
                return false;
            }
            for(unsigned i = k + 1; i < n; ++i){
                double sum = 0;
                for(unsigned m = 0; m < k; ++m){
                    sum += l[i][m] * u[m][k];
                }
                l[i][k] = (A[i][k] - sum) / u[k][k];
            }
        }
    }
    double* y = new double[n];
    y[0] = b[0];
    for(unsigned i = 1; i < n; ++i){
        y[i] = b[i];
        for(unsigned j = 0; j < i; ++j){
            y[i] -= l[i][j] * y[j];
        }
    }
    x[n - 1] = y[n - 1] / u[n - 1][n - 1];
    for(unsigned i = n - 1; i > 0; --i){
        double sum = 0;
        for(unsigned j = i; j < n; ++j){
            sum += u[i - 1][j] * x[j];
        }
        x[i - 1] = (y[i - 1] - sum) / u[i - 1][i - 1];
    }
    deposeMatrix(l, n);
    deposeMatrix(u, n);
    delete[] y;
    return true;
}


double maxDelta(double* x, double* y, unsigned n){
    double delta = std::abs(x[0] - y[0]);
    for(unsigned i = 1; i < n; ++i){
        double tmp = std::abs(x[i] - y[i]);
        if(tmp > delta){
            delta = tmp;
        }
    }
    return delta;
}


int Jacobi(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x){
    unsigned k = 1;
    double* g = new double[n];
    for(unsigned i = 0; i < n; ++i){
        if(std::abs(A[i][i]) < e){
            delete[] g;
            return -1;
        }
        double T = A[i][i];
        for(unsigned j = 0; j < n; ++j){
            A[i][j] = -A[i][j] / T;
        }
        A[i][i] = 0;
        g[i] = b[i] / T;
    }
    while(true){
        for(unsigned i = 0; i < n; ++i){
            x[i] = g[i];
            for(unsigned j = 0; j < n; ++j){
                if(j != i){
                    x[i] += A[i][j] * y[j];
                }
            }
        }
        double delta = maxDelta(x, y, n);
        if(delta < e){
            return k; 
        }
        if(k < M){
            ++k;
            for(unsigned i = 0; i < n; ++i){
                y[i] = x[i];
            }
        } else {
            return -1;
        }
    }
}


int Seidel(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x){
    unsigned k = 1;
    double* g = new double[n];
    for(unsigned i = 0; i < n; ++i){
        x[i] = y[i];
    }
    for(unsigned i = 0; i < n; ++i){
        if(std::abs(A[i][i]) < e){
            delete[] g;
            return -1;
        }
        double T = A[i][i];
        for(unsigned j = 0; j < n; ++j){
            A[i][j] = -A[i][j] / T;
        }
        A[i][i] = 0;
        g[i] = b[i] / T;
    }
    while(true){
        for(unsigned i = 0; i < n; ++i){
            x[i] = g[i];
            for(unsigned j = 0; j < n; ++j){
                if(j != i){
                    x[i] += A[i][j] * x[j];
                }
            }
        }
        double delta = maxDelta(x, y, n);
        if(delta < e){
            return k; 
        }
        if(k < M){
            ++k;
            for(unsigned i = 0; i < n; ++i){
                y[i] = x[i];
            }
        } else {
            return -1;
        }
    }
}


int SOR(double** A, double* b, double* y, double w, double e, unsigned n, 
unsigned M, double* x){
    unsigned k = 1;
    double* g = new double[n];
    for(unsigned i = 0; i < n; ++i){
        x[i] = y[i];
    }
    for(unsigned i = 0; i < n; ++i){
        if(std::abs(A[i][i]) < e){
            delete[] g;
            return -1;
        }
        double T = A[i][i];
        for(unsigned j = 0; j < n; ++j){
            A[i][j] = -w * A[i][j] / T;
        }
        A[i][i] = 1 - w;
        g[i] = w * b[i] / T;
    }
    while(true){
        for(unsigned i = 0; i < n; ++i){
            x[i] = g[i] + A[i][i] * x[i];
            for(unsigned j = 0; j < n; ++j){
                if (j != i) {
                    x[i] += A[i][j] * x[j];
                }
            }
        }   
        double delta = maxDelta(x, y, n);
        if(delta < e){
            return k; 
        }
        if(k < M){
            ++k;
            for(unsigned i = 0; i < n; ++i){
                y[i] = x[i];
            }
        } else {
            return -1;
        }
    }
}
