#include<vector>
#include<cmath>

#include"numerical_integration.h"

long double AutoStepTrapezoidalMethod(long double a, long double b,
                                      std::function<long double(long double)> f,
                                      long double deta){
	long double h = (b - a) / 2;
	long double T1 = (f(a) + f(b)) * h;
	long n = 1;
	while(1){
		long double T0 = T1;
		long double S = 0;
		for(long k = 1; k < n + 1; ++k){
			S += f(a + (2 * k - 1) * h / n);
		}
		T1 = T0 / 2 + S * h / n;
		if(std::abs(T1 - T0) < 3 * deta){
			return T1;
		} else {
			n *= 2;
		}
	}
}	


long double Romberg(long double a, long double b,
                    std::function<long double(long double)> f,
                    long double deta){
	std::vector<std::vector<long double>> T;
	T.push_back(std::vector<long double>());
	T[0].push_back((b - a) / 2 * (f(a) - f(b)));
	unsigned k = 1;
	while(1){
		unsigned long n1 = static_cast<unsigned long>(std::pow(2, k));
		long double m = 0;
        for (unsigned long i = 1; i < static_cast<unsigned long>(n1 / 2); ++i) {
            m += f(a + (2 * i - 1) * (b - a) / n1);
        }
		T[0].push_back(0.5 * (T[0][k - 1] + (b - a) / (n1 / 2) * m));
		for(unsigned m = 1; m < k + 1; ++m){
			if(m == k){
				T.push_back(std::vector<long double>());
			}
			long double n2 = std::pow(4.0L, m);
			T[m].push_back((n2 * T[m - 1][k - m + 1] - T[m - 1][k - m]) / (n2 - 1));
		}
		if(std::abs(T[k][0] - T[k - 1][0]) < deta){
			return T[k][0];
        } else {
            k += 1;
        }
	}
}
		
