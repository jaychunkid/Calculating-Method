#include<functional>
#include<vector>
#include<cmath>

#include"nonlinear_equation.h"

double Halving(std::function<double(double)> f, double a, double b, 
               double e1, double e2){
    unsigned k = 0;
    double x = 0;
    while(true){
        x = (a + b) / 2;
        double fx = f(x);
        if(std::abs(fx) < e1){
            return x;
        }
        if(f(a) * f(x) < 0){
            b = x;
        } else {
            a = x;
        }
        if(b - a <= e2){
            return (a + b) / 2;
        } else {
            ++k;
        }
    }
}


double Relaxation(std::function<double(double)> f, std::function<double(double)> ff, 
                  double x, double e){
    unsigned k = 0;
    double w = 0;
    double nx = 0;
    while(true){
        w = 1 / (1 - ff(x));
        nx = (1 - w) * x + w * f(x);
        if(std::abs(nx - x) < e){
            return nx;
        } else {
            ++k;
            x = nx;
        }
    }
}


double Newton(std::function<double(double)> f, std::function<double(double)> ff, 
              double x, double e){
    unsigned k = 0;
    double nx = 0;
    while(true){
        nx = x - f(x) / ff(x);
        if(std::abs(nx - x) < e){
            return nx;
        } else {
            ++k;
            x = nx;
        }
    }
}