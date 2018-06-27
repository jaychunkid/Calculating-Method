/*������ֵ�㷨*/
#ifndef ALGEBRAIC_INTERPOLATION_H
#define ALGEBRAIC_INTERPOLATION_H

//Lagrange��n�β�ֵ�㷨
long double Lagrange(unsigned n, const long double *x, const long double *y, 
                     long double xx);

//Newton��n�β�ֵ�㷨
long double Newton(unsigned n, const long double *x, const long double *y, 
                   long double xx);

//����������ֵ�㷨�������߽�����Ϊf�Ķ��ε�������x0��xn����ֵ��Ϊ0
long double CubicSpline(unsigned size, const long double *x, 
                        const long double *y, long double xx);

//����������ֵ�㷨�������߽�����Ϊf��һ�ε�������x0��xn����ֵ�ֱ�Ϊm0��mn
long double CubicSpline(unsigned size, const long double *x, 
                        const long double *y, long double m0, 
                        long double mn, long double xx);

#endif // !ALGEBRAIC_INTERPOLATION_H

