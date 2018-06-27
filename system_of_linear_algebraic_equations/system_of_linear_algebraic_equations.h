/*解线性代数方程组算法*/
#ifndef SYSTEM_OF_LINEAR_ALGEBRAIC_EQUATIONS_H
#define SYSTEM_OF_LINEAR_ALGEBRAIC_EQUATIONS_H

//顺序高斯消去法
bool Gauss1(double** A, double* b, double e, unsigned n, double* x);

//列主元高斯消去法
bool Gauss2(double** A, double* B, double e, unsigned n, double* x);

//全主元高斯消去法
bool Gauss3(double** A, double* b, double e, unsigned n, double* x);

//直接LU分解法
bool LU(double** A, double* b, double e, unsigned n, double* x);

//Jacobi迭代
int Jacobi(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x);

//Seidel迭代
int Seidel(double** A, double* b, double* y, double e, unsigned n, 
unsigned M, double* x);

//松弛迭代
int SOR(double** A, double* b, double* y, double w, double e, unsigned n, 
unsigned M, double* x);

#endif