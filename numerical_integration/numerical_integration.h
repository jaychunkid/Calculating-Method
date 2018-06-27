/*数值积分算法*/
#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include<functional>

//自动选取步长梯形法
long double AutoStepTrapezoidalMethod(long double a, long double b, 
                                      std::function<long double(long double)> f, 
                                      long double deta);

//Romberg求积法
long double Romberg(long double a, long double b, 
                    std::function<long double(long double)> f, 
                    long double deta);

#endif // !NUMERICAL_INTEGRATION_H