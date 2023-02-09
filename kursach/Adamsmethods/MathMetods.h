#ifndef MATHMETODS_H
#define MATHMETODS_H

#include<vector>

//тип функций возвращающий (double, double) -> (double)
using func = double (*)(double, double);

//тип функций возвращающий (double, int) -> (double)
using func1 = double (*)(double, int);

//Метод Симпсона для интегрирования функции f на отрезке [a, b] с точностью eps
double integrate(double a, double b, double eps, int k, func1 f);

//Метод Рюнге-Кутта четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 на отрезке [x0, xn] с шагом h
std::vector<double> methodRungeKutta(double y0, double x0, double xn, double h, func f);

//Интерполяционный метод Адамса четвертого порядка точности для вычисления коррекционного значения yn с заданной точностью eps
void correctionMethodAdams(std::vector<double> &Values, const std::vector<double> &X, double h, double eps, func f);

//Экстраполяционный метод Адамса четвертого порядка точности для нахождения значений приближенного
//решения задачи Коши с начальным значением y0 и найденными значениями yi (i = 1, ..., n) по методу Рюнге-Кутты
void methodAdams(std::vector<double> &Values, double x0, double xn, double b, double h, double eps, func f);

#endif
