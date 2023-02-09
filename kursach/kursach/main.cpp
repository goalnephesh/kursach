#include<iostream>
#include<cmath>
#include"mathmethods.h"

double f(double x, double y){
	return -y * std::cos(x) + std::sin(x) * std::cos(x);
}

int main (int argc, char *argv[]) {
	auto Values {methodRungeKutta(-1, 0, 0.5, 0.1, f)};
//	std::cout << "Таблица значений рунге-кутты:" << std::endl;
//	for(auto &val: Values) std::cout << val << ' ';
	std::cout << std::endl;
	methodAdams(Values, 0, 0.5, 20, 0.1, f);
/*	for(auto i {0}; i < Values.size() - 1; ++i) printf("Y[%d] = %.10f\n", i, Values[i]);*/
	auto x{0.};
	for(auto i {0}; i < Values.size() - 1; ++i){
		std::cout << "|Y[" << i << "] - phi(x)| = " << std::abs(Values[i] -(std::sin(x) - 1) ) << std::endl;
		x += 0.1;
	}
	/*for(auto i {0}; i < Values.size() - 1; ++i) printf("DeltaY[%d] = %.6f\n", i, Values[i + 1] - Values[i]);*/
	return 0;
}	
	
