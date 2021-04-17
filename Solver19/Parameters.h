#pragma once
#include <functional>
#include <vector>

enum kind_of_solution {
	FIRST,  // ���������� ������ �����
	SECOND, // ���������� ���������� �������, ���������� ���� �� ������ ��������
	THIRD   // ���������� ��� ������������� ��������
};

class Parameters
{
	static double get_value(double t, size_t param_num);
public:	
	static kind_of_solution kind;
	static std::vector<std::function<double(double)>> smooth_params;
	static std::vector<double> const_params;
	static std::vector<double> points;
	static std::vector<std::vector<double>> piecewise_linear_params;
	static double evaluate(double t, size_t param_num);
//	static double kappa;
};


