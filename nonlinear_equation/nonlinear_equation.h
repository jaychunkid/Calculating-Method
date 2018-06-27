/*������Է����㷨*/
#ifndef NONLINEAR_EQUATION_H
#define NONLINEAR_EQUATION_H

#include<functional>

//�Էַ�
double Halving(std::function<double(double)> f, double a, double b,
               double e1, double e2);

//�ɳڷ�
double Relaxation(std::function<double(double)> f, std::function<double(double)> ff,
                  double x, double e);

//ţ�ٷ�
double Newton(std::function<double(double)> f, std::function<double(double)> ff,
              double x, double e);

#endif