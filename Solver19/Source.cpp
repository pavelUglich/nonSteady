#include <vector>
#include <cstdlib>
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include "OdeSolver.h"
#include <ctime>

#include "integrals.h"


using namespace std;
std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;








int main()
{
	setlocale(0, "");
	Parameters::kind = FIRST;
	const auto t = 1.0;
	const double b = 20.0;
	const size_t length = 50;
	const double h = b / length;
	//const std::complex<double> alpha = { 1,0 };
	std::vector<std::complex<double>> vector(length);
	const unsigned int start_time = clock(); // начальное время
	#pragma omp parallel for
	for (int i = 1; i <= length; i++)
	{
		const double eps = 0.001;
		const std::complex<double> alpha = { i*h,0 };
		vector[i - 1] = back_integrand(alpha, t, eps);
		cout << i << " \n";
	}
	const unsigned int end_time = clock(); // начальное время
	cout << end_time - start_time << endl;
	plotTheWaveField(h, { {"black", vector}/*, {"red", vector1}*/}, "xxx.txt");
	system("pause");
}



//Аналитическое выражение



//Условные вычисления




