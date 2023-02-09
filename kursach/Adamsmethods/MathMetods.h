#ifndef MATHMETODS_H
#define MATHMETODS_H

#include<vector>

//��� ������� ������������ (double, double) -> (double)
using func = double (*)(double, double);

//��� ������� ������������ (double, int) -> (double)
using func1 = double (*)(double, int);

//����� �������� ��� �������������� ������� f �� ������� [a, b] � ��������� eps
double integrate(double a, double b, double eps, int k, func1 f);

//����� �����-����� ���������� ������� �������� ��� ���������� �������� �������������
//������� ������ ���� � ��������� ��������� y0 �� ������� [x0, xn] � ����� h
std::vector<double> methodRungeKutta(double y0, double x0, double xn, double h, func f);

//���������������� ����� ������ ���������� ������� �������� ��� ���������� �������������� �������� yn � �������� ��������� eps
void correctionMethodAdams(std::vector<double> &Values, const std::vector<double> &X, double h, double eps, func f);

//����������������� ����� ������ ���������� ������� �������� ��� ���������� �������� �������������
//������� ������ ���� � ��������� ��������� y0 � ���������� ���������� yi (i = 1, ..., n) �� ������ �����-�����
void methodAdams(std::vector<double> &Values, double x0, double xn, double b, double h, double eps, func f);

#endif
