#include"mathmethods.h"
#include<cmath>
#include<iostream>
//Метод Рунге-Кутта четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 на отрезке [x0, X] с шагом h
std::vector<double> methodRungeKutta(double y0, double x0, double X, double h, func f){
	auto yi {y0},
		xi {x0};
	std::vector<double> Values;
	Values.push_back(yi);
	while(xi < X){
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
//			std::cout << "Cor Eta[" << N - i - 1 << "] = " << Eta[i] << std::endl;
		}
		for(auto i {0}; i < Delta1Size; ++i){
			Delta1[i] = Eta[i] - Eta[i + 1];
//			std::cout << "Cor Delta1[" << Delta1Size - i - 1 << "] = " << Delta1[i] << std::endl;
		}
		for(auto i {0}; i < Delta2Size; ++i){
			Delta2[i] = Delta1[i] - Delta1[i + 1];
//			std::cout << "Cor Delta2[" << Delta2Size - i - 1 << "] = " << Delta2[i] << std::endl;
		}
		for(auto i {0}; i < Delta3Size; ++i){
			Delta3[i] = Delta2[i] - Delta2[i + 1];
//			std::cout << "Cor Delta3[" << Delta3Size - i - 1 << "] = " << Delta3[i] << std::endl;
		}
		Delta4 = Delta3[0] - Delta3[1];
//		std::cout << "Cor Delta4 = " << Delta4 << std::endl;
		auto yn {Values[N - 2] + Eta[0] - Delta1[0] / 2 - Delta2[0] / 12 - Delta3[0] / 24 - 19 * Delta4 / 720};
		flag = std::abs(Values[N - 1] - yn) < eps;
		Values[N - 1] = yn;
//		printf("Corrected Values[%d] = %f\n", N - 1, Values[N - 1]);
	}
}

//Экстраполяционный метод Адамса четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 и найденными значениями yi (i = 1, ..., n) по методу Рюнге-Кутты
void methodAdams(std::vector<double> &Values, double x0, double xn, double b, double h, func f){
	std::vector<double> X(Values.size());
	std::vector<double> Eta(5);
	std::vector<double> Delta1(4);
	std::vector<double> Delta2(3);
	std::vector<double> Delta3(2);
	auto EtaSize {Eta.size()},
		Delta1Size {Delta1.size()},
		Delta2Size {Delta2.size()},
		Delta3Size {Delta3.size()};
	double Delta4 {0.};
	auto N {Values.size()};
	auto xi {x0};
	for(auto i {0}; i < N; ++i){
		X[i] = xi;
		xi += h;
	}
	xi -= h;
	while(xi < b){
		for(auto i {0}; i < EtaSize; ++i){
			Eta[i] = h * f(X[N - i - 1], Values[N - i - 1]);
//			std::cout << "Eta[" << N - i - 1 << "] = " << Eta[i] << std::endl;
		}
		for(auto i {0}; i < Delta1Size; ++i){
			Delta1[i] = Eta[i] - Eta[i + 1];
//			std::cout << "Delta1[" << Delta1Size - i - 1 << "] = " << Delta1[i] << std::endl;
		}
		for(auto i {0}; i < Delta2Size; ++i){
			Delta2[i] = Delta1[i] - Delta1[i + 1];
//			std::cout << "Delta2[" << Delta2Size - i - 1 << "] = " << Delta2[i] << std::endl;
		}
		for(auto i {0}; i < Delta3Size; ++i){
			Delta3[i] = Delta2[i] - Delta2[i + 1];
//			std::cout << "Delta3[" << Delta3Size - i - 1 << "] = " << Delta3[i] << std::endl;
		}
		Delta4 = Delta3[0] - Delta3[1];
//		std::cout << "Delta4 = " << Delta4 << std::endl;
		Values.push_back(Values[N - 1] + Eta[0] + Delta1[0] / 2 + 5 * Delta2[0] / 12 + 3 * Delta3[0] / 8  + 251 * Delta4 / 720);
		xi += h;
		X.push_back(xi);
		N++;
//		printf("Values[%d] = %.6f\n", N - 1,Values[N - 1]);
		correctionMethodAdams(Values, X, h, h * h * h * h * h * h * 95 / 283, f);
	}
}
