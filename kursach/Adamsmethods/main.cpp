#include<iostream>
#include<cmath>
#include"MathMetods.h"

double f(double x, double y){
	return -y * std::cos(x) + std::sin(x) * std::cos(x);
}

double f1(double x, int k){
	return x * x;
}

int main (int argc, char *argv[]) {
//	auto Values {methodRungeKutta(-1, 0, 1.5, 0.1, f)};
//	methodAdams(Values, 0, 1.5, 2, 0.1, f);
//	for(auto &val: Values) std::cout << val << ' ';
	std::cout << integrate(0, 3, 0.001, f1);
	return 0;
}

