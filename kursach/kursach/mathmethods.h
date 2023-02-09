#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#include<vector>

//��� ������� ������������ (double) -> (double)
using func = double (*)(double, double);

//����� �����-����� ���������� ������� �������� ��� ���������� �������� �������������
//������� ������ ���� � ��������� ��������� y0 �� ������� [x0, X] � ����� h
std::vector<double> methodRungeKutta(double y0, double x0, double xn, double h, func f);

//���������������� ����� ������ ���������� ������� �������� ��� ���������� �������������� �������� yn � �������� ��������� eps
void correctionMethodAdams(std::vector<double> &Values, const std::vector<double> &X, double h, double eps, func f);

//����������������� ����� ������ ���������� ������� �������� ��� ���������� �������� �������������
//������� ������ ���� � ��������� ��������� y0 � ���������� ���������� yi (i = 1, ..., n) �� ������ �����-�����
void methodAdams(std::vector<double> &Values, double x0, double xn, double b, double h, func f);

#endif
