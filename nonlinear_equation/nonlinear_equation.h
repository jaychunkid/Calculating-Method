/*解非线性方程算法*/
#ifndef NONLINEAR_EQUATION_H
#define NONLINEAR_EQUATION_H

#include<functional>

//对分法
double Halving(std::function<double(double)> f, double a, double b,
               double e1, double e2);

//松弛法
double Relaxation(std::function<double(double)> f, std::function<double(double)> ff,
                  double x, double e);

//牛顿法
double Newton(std::function<double(double)> f, std::function<double(double)> ff,
              double x, double e);

#endif