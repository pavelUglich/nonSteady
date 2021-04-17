#include <vector>
#include <cstdlib>
#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include "OdeSolver.h"

using namespace std;

std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;
const double pi = 3.1415926535897932384626433832795;


const double nodes[] = { -0.989400934991650 , -0.944575023073233,
-0.865631202387832, -0.755404408355003 , -0.617876244402644 ,
-0.458016777657227 , -0.281603550779259, -0.950125098376374e-1,
0.950125098376374e-1 ,0.281603550779259, 0.458016777657227, 0.617876244402644,
0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650 };
const double weights[] = { 0.0271524594117540, 0.0622535239386475,
0.0951585116824939, 0.124628971255534, 0.149595988816577, 0.169156519395002,
0.182603415044922, 0.189450610455068, 0.189450610455068, 0.18260341504492,
0.16915651939500, 0.14959598881657, 0.124628971255534, 0.0951585116824939,
0.0622535239386475, 0.0271524594117540 };

std::complex<double> evaluate_integral(
	const std::function<std::complex<double>(double)> integrand, double a = -1,
	double b = 1);
std::complex<double> integrand(std::complex<double> alpha,
	std::complex<double> s, double t);
std::complex<double> contour(double s, double zeta = 0);
std::complex<double> new_integrand(std::complex<double> alpha, double s, double t);
std::complex<double> back_integrand(std::complex<double> alpha, /*std::complex<double> s,*/ double t, double eps);
std::complex<double> back_integrands(std::complex<double> alpha, double t, double eps);

/* std::complex<double> Solve1(double s, double t); */


void showVector(double step,
	const std::vector<std::complex<double>>& complexValuedVector,
	const std::function<double(std::complex<double>)>& mapping, std::ostream& stream, double a = 0)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		const auto x = a + (i + 0.5) * step;
		stream << "(" << x << ", " << mapping(complexValuedVector[i]) << ") ";
	}
}

void addTheCurve(double step, const std::vector<std::complex<double>>& fieldHomogeneous,
	std::ofstream& stream, const std::string& color, double a = 0)
{
	stream << " \\addplot[line width = 0.25mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.real(); }, stream, a);
	stream << " };\n";
	stream << " \\addplot[smooth, dashed,";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.imag(); }, stream, a);
	stream << " };\n";
}


void plotTheWaveField(double step, const std::map<std::string,
	std::vector<std::complex<double>>>& vectors, const std::string& fileName, double a = 0)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	for (auto item : vectors)
	{
		addTheCurve(step, item.second, stream, item.first, a);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}


int main()
{
	setlocale(0, "");
	Parameters::kind = FIRST;
	auto t = 0.01;
	double a = 1.0;
	double b = 10.0;
	double h = 0.1;
	double eps = 0.001;

	std::vector<std::complex<double>> vector;

	for (; a < b; a += h)
	{
		std::complex<double> alpha = { a,0 };
		//std::complex<double> alpha = { a,1 };
		//const auto value = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, 0, 2 * pi / t);		
		vector.push_back(back_integrands(alpha, t, 0.0001));
		//vector.push_back(new_integrand(alpha, a, t));

	}

	plotTheWaveField(h, { { "black", vector } }, "1.txt");


	system("pause");
}

// квадратурная формула Гаусса
std::complex<double> evaluate_integral(const std::function<std::complex<double>(double)> integrand, double a, double b)
{
	std::complex<double> value = 0;
	const size_t size = sizeof(nodes) / sizeof(double);
	for (size_t i = 0; i < size; i++)
	{
		double x = (a + b) / 2 + (b - a) / 2 * nodes[i];
		value += weights[i] * integrand(x);
	}
	return (b - a) * value / 2.0;
}

// подынтегральное выражение

std::complex<double> integrand(std::complex<double> alpha, std::complex<double> s, double t)
{
	const std::function<double(double)> mu = [](double x) {return 1; };
	const std::function<double(double)> rho = [](double x) {return 1; };
	const boundary_value_problem<std::complex<double>> bvp = { // создаётся объект типа краевая задача
		{
			[=](double x, const std::vector<std::complex<double>>& v) {return 1 / mu(x) * v[1]; }, // правые части уравнений системы
			[=](double x, const std::vector<std::complex<double>>& v) {return (rho(x) * s * s + alpha * alpha * mu(x)) * v[0]; }

		},{ { 0,0 } },{ { 1, 1.0 } } }; // краевые условия
	return bvp.solve()[0] / s * exp(s * t); // набор точек
}

//Аналитическое выражение

std::complex<double> analytical(std::complex<double> alpha, std::complex<double> s, double t)
{
	const auto gamma = sqrt(alpha * alpha + s * s);
	const auto e = exp(gamma);
	const auto sh = (e - 1.0 / e) / 2.0;
	const auto ch = (e + 1.0 / e) / 2.0;
	return  exp(s * t) * (1.0 / s) * sh / (gamma * ch);
}

//Условные вычисления

std::complex<double> back_integrands(std::complex<double> alpha, double t, double eps)
{
	std::complex<double> sm = { 0.0, 0.0 };
	int i = 0;

	std::complex<double> value = { 1, 0 };

	while (value.real() * value.real() + value.imag() * value.imag() > eps * eps)
	{
		if (i >= 1000) break;
		const auto tmp = (2.0 * i - 1) * pi / 2.0;
		const auto sqr = sqrt(alpha * alpha + tmp * tmp);
		value = 2.0 * sin(tmp) * sin(t * sqr) / sqr;
		if (i % 2 != 0) value *= -1;
		i++;
		sm += value;
	}
	return sm;
}


//std::complex<double> back_integrands(std::complex<double> alpha, /*std::complex<double> s,*/ double t, double eps)
/*
{
std::complex<double> sum = { 0.0, 0.0 };
const auto sqr = sqrt((alpha * alpha) + (((2 - 1) * pi) / 2) * (((2 - 1) * pi) / 2));
return 2 * sin(((2 - 1) / 2) * (pi)) * ((sin(t * sqr)) / (sqr));
}
*/



// контур интегрирования
std::complex<double> contour(double s, double zeta)
{
	return{ zeta, s };
}

// интеграл по контуру
std::complex<double> new_integrand(std::complex<double> alpha, double s, double t)
{
	std::complex<double>point = contour(s, 0.1);
	return integrand(alpha, point, t);
}

// Обратное преобразование Лапласа
std::complex<double> back_integrand(std::complex<double> alpha, /*std::complex<double> s,*/ double t, double eps)
{
	std::complex<double> sum = { 0.0, 0.0 };
	const std::complex<double> denumerator = { 0, 2 * pi };
	int i = 0;
	std::complex<double> value = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, 0, 2 * pi / t);
	while (value.real() * value.real() + value.imag() * value.imag() > eps * eps)
	{
		sum += value;
		i++;
		if (i >= 1000) break;
		value = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, i * 2 * pi / t, (i + 1) * 2 * pi / t);
	}
	return sum / denumerator;
}