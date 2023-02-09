#include"MathMetods.h"
#include<cmath>
#include<iostream>


double integrate(double a, double b, double eps, int k, func1 f){
	auto sum1 {0.0};
	auto N1 {(b - a) / eps};
	auto fin {false};
	auto xi {a + 1 / N1};
	while(xi + 1 / N1 <= b){
		sum1 += f(xi - 1 / N1, k) + 4 * f(xi, k) + f(xi + 1 / N1, k);
		xi += 2 / N1;
	}
	sum1 *= 1 / (3 * N1);
	while(!fin){
		auto N2 {2 * N1};
		auto sum2 {0.0};
		xi = a + 1 / N2;
		while(xi + 1 / N2 < b){
			sum2 += f(xi - 1 / N2, k) + 4 * f(xi, k) + f(xi + 1 / N2, k);
			xi += 2 / N2;
		}
		sum2 *= 1 / (3 * N2);
		fin = std::abs(sum1 - sum2) < eps;
		sum1 = sum2;
		N1 = N2;
	}
	return sum1;
}

double coef(double u, int k){
	auto fact {1.0};
	auto num {1.0};
	auto i {0};
	while(i < k){
		num *= u + i;
		++i;
		fact *= i;
	}
	return num / fact;
}

//Метод Рюнге-Кутта четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 на отрезке [x0, xn] с шагом h
std::vector<double> methodRungeKutta(double y0, double x0, double xn, double h, func f){
	auto yi {y0},
		 xi {x0};
	std::vector<double> Values;
	Values.push_back(yi);
	while(xi <= xn){
		auto K1 {h * f(xi, yi)};
		auto K2 {h * f(xi + h / 2, yi + K1 / 2)};
		auto K3 {h * f(xi + h / 2, yi + K2 / 2)};
		auto K4 {h * f(xi + h, yi + K3)};
		auto DeltaY {(K1 + 2 * K2 + 2 * K3 + K4) / 6};
		yi += DeltaY;
		xi += h;
		Values.push_back(yi);
	}
	return Values;
}

//Интерполяционный метод Адамса четвертого порядка точности для вычисления коррекционного значения yn с заданной точностью eps
void correctionMethodAdams(std::vector<double> &Values, const std::vector<double> &X, double h, double eps, func f){
	std::vector<double> Eta(5);
	std::vector<double> Delta1(4);
	std::vector<double> Delta2(3);
	std::vector<double> Delta3(2);
	auto EtaSize {Eta.size()},
		Delta1Size {Delta1.size()},
		Delta2Size {Delta2.size()},
		Delta3Size {Delta3.size()};
	auto N {Values.size()};
	auto flag {false};
	double Delta4 {0.0};
	while(!flag){
		for(auto i {0}; i < EtaSize; ++i){
			Eta[i] = h * f(X[N - i - 1], Values[N - i - 1]);
		}
		for(auto i {0}; i < Delta1Size; ++i){
			Delta1[i] = Eta[i] - Eta[i + 1];
		}
		for(auto i {0}; i < Delta2Size; ++i){
			Delta2[i] = Delta1[i] - Delta1[i + 1];
		}
		for(auto i {0}; i < Delta3Size; ++i){
			Delta3[i] = Delta2[i] - Delta2[i + 1];
		}
		Delta4 = Delta3[0] - Delta3[1];
		auto yn {Values[N - 2] + Eta[0] - Delta1[0] / 2 - Delta2[0] / 12 - Delta3[0] / 24 - 19 * Delta4 / 720};
		flag = std::abs(Values[N - 1] - yn) < eps;
		Values[N - 1] = yn;
	}
}

//Экстраполяционный метод Адамса четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 и найденными значениями yi (i = 1, ..., n) по методу Рюнге-Кутты
void methodAdams(std::vector<double> &Values, double x0, double xn, double b, double h, double eps, func f){
	auto N {Values.size()};
	auto flag {true};
	std::vector<double> X(N);
	
	std::vector<std::vector<double>> Delta(N);
	for(auto i {0}; i < N; ++i){
		Delta[i] = std::vector<double>(N - i);
	}
//	std::vector<double> Eta(5);
//	std::vector<double> Delta1(4);
//	std::vector<double> Delta2(3);
//	std::vector<double> Delta3(2);
//	auto EtaSize {Eta.size()},
//		  Delta1Size {Delta1.size()},
//		  Delta2Size {Delta2.size()},
//		  Delta3Size {Delta3.size()};
//	double Delta4 {0.};
	auto xi {x0};
	for(auto i {0}; i < N; ++i){
		X[i] = xi;
		xi += h;
	}
	xi -= h;
	for(auto i {0}; i < N; ++i)
		Delta[0][i] = h * f(X[N - i - 1], Values[N - i - 1]);
	auto k {0};
	for(auto i {1}; flag && i < N; ++i){
		for(auto j {0}; flag && j < N - i; ++j){
			Delta[i][j] = Delta[i - 1][j] - Delta[i - 1][j + 1];
			flag = Delta[i][j] < eps;
		}
		k++;
	}
	std::vector<double> C(k + 2);
	for(auto i {0}; i < k + 2; ++i)
		C[i] = integrate(0, 1, eps, i, coef);
	while(xi <= b){
		for(auto i {0}; i < EtaSize; ++i)
			Eta[i] = h * f(X[N - i - 1], Values[N - i - 1]);
		for(auto i {0}; i < Delta1Size; ++i)
			Delta1[i] = Eta[i] - Eta[i + 1];
		for(auto i {0}; i < Delta2Size; ++i)
			Delta2[i] = Delta1[i] - Delta1[i + 1];
		for(auto i {0}; i < Delta3Size; ++i)
			Delta3[i] = Delta2[i] - Delta2[i + 1];
		Delta4 = Delta3[0] - Delta3[1];
		Values.push_back(Values[N - 1] + Eta[0] + Delta1[0] / 2 + 5 * Delta2[0] / 12 + 3 * Delta3[0] / 8  + 251 * Delta4 / 720);
		xi += h;
		X.push_back(xi);
		N++;
		correctionMethodAdams(Values, X, h, h * h * h * h * h * h * 95 / 283, f);
	}
}
