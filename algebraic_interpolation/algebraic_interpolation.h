/*代数插值算法*/
#ifndef ALGEBRAIC_INTERPOLATION_H
#define ALGEBRAIC_INTERPOLATION_H

//Lagrange型n次插值算法
long double Lagrange(unsigned n, const long double *x, const long double *y, 
                     long double xx);

//Newton型n次插值算法
long double Newton(unsigned n, const long double *x, const long double *y, 
                   long double xx);

//三次样条插值算法，给定边界条件为f的二次导函数在x0、xn处的值均为0
long double CubicSpline(unsigned size, const long double *x, 
                        const long double *y, long double xx);

//三次样条插值算法，给定边界条件为f的一次导函数在x0、xn处的值分别为m0、mn
long double CubicSpline(unsigned size, const long double *x, 
                        const long double *y, long double m0, 
                        long double mn, long double xx);

#endif // !ALGEBRAIC_INTERPOLATION_H

