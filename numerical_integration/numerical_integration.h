/*��ֵ�����㷨*/
#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include<functional>

//�Զ�ѡȡ�������η�
long double AutoStepTrapezoidalMethod(long double a, long double b, 
                                      std::function<long double(long double)> f, 
                                      long double deta);

//Romberg�����
long double Romberg(long double a, long double b, 
                    std::function<long double(long double)> f, 
                    long double deta);

#endif // !NUMERICAL_INTEGRATION_H