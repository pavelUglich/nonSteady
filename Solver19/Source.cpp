#include <vector>
#include <cstdlib>
#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include "OdeSolver.h"
#include <omp.h>
#include <ctime>



using namespace std;

std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;
const double pi = 3.1415926535897932384626433832795;


/**
 * \brief узловые точки 
 */
const double nodes[] = { -0.989400934991650 , -0.944575023073233,
-0.865631202387832, -0.755404408355003 , -0.617876244402644 ,
-0.458016777657227 , -0.281603550779259, -0.950125098376374e-1,
0.950125098376374e-1 ,0.281603550779259, 0.458016777657227, 0.617876244402644,
0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650 };


/**
 * \brief весовые коэффициенты
 */
const double weights[] = { 0.0271524594117540, 0.0622535239386475,
0.0951585116824939, 0.124628971255534, 0.149595988816577, 0.169156519395002,
0.182603415044922, 0.189450610455068, 0.189450610455068, 0.18260341504492,
0.16915651939500, 0.14959598881657, 0.124628971255534, 0.0951585116824939,
0.0622535239386475, 0.0271524594117540 };


/**
 * \brief отыскание интеграла
 * \param integrand подынтегральное выржение
 * \param a лева€ граница отрезка 
 * \param b права€ граница отрезка
 * \return значение интеграла
 */
std::complex<double> evaluate_integral(
	const std::function<std::complex<double>(double)>& integrand, double a = -1,
	double b = 1);


/**
 * \brief подынтегральное выражение
 * \param alpha параметр преобразовани€ Ћапласа
 * \param s параметр преобразовани€ ‘урье
 * \param t врем€
 * \return подынтегральное выражение
 */
std::complex<double> integrand(std::complex<double> alpha,
	std::complex<double> s, double t);

/**
 * \brief контур интегрировани€
 * \param s параметр преобразовани€ Ћапласа
 * \param zeta абсцисса контура интегрировани€
 * \return точка контура
 */
std::complex<double> contour(double s, double zeta = 0);

/**
 * \brief подынтегральное выражение вдольконтура интегрировани€
 * \param alpha параметр преобразовани€ ‘урье
 * \param s параметр преобразовани€ Ћапласа
 * \param t врем€
 * \return значение
 */
std::complex<double> new_integrand(std::complex<double> alpha, double s, double t);

/**
 * \brief обратное преобразование Ћапласа
 * \param alpha параметр преобразовани€ ‘урье
 * \param t врем€
 * \param eps точность вычислений
 * \return значение обратной трансформанты
 */
std::complex<double> back_integrand(std::complex<double> alpha, double t, 
	double eps);

/**
 * \brief обратное преобразование Ћапласа в однородном случае (аналитическое
 * выражение)
 * \param alpha параметр преобразовани€ ‘урье
 * \param t врем€
 * \param eps точность вычислений
 * \return значение
 */
std::complex<double> back_integrands(std::complex<double> alpha, double t, 
	double eps);


/**
 * \brief добавить кривую на график
 * \param step шаг
 * \param complexValuedVector вектор
 * \param mapping отображение
 * \param stream файловый поток
 * \param a лева€ граница отрезка
 */
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

/**
 * \brief  добавить кривую на график
 * \param step шаг
 * \param fieldHomogeneous вектор
 * \param stream файловый поток
 * \param color цвет
 * \param a лева€ граница
 */
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

/**
 * \brief построить диаграмму
 * \param step шаг
 * \param vectors набор векторов
 * \param fileName им€ файла
 * \param a лева€ граница
 */
void plotTheWaveField(double step, const std::map<std::string,
	std::vector<std::complex<double>>>& vectors, const std::string& fileName, double a = 0)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	for (const auto& item : vectors)
	{
		addTheCurve(step, item.second, stream, item.first, a);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

std::complex<double> analytical(std::complex<double> alpha, std::complex<double> s, double t);

int main()
{
	setlocale(0, "");
	Parameters::kind = FIRST;
	const auto t = 1.0;
	const double b = 10.0;
	double eps = 0.001;
	const size_t length = 50;
	const double h = b / length;
	const std::complex<double> alpha = { 1,0 };
	std::vector<std::complex<double>> vector(length);
	unsigned int start_time = clock(); // начальное врем€
	#pragma omp parallel for
	for (size_t i = 1; i <= length; i++)
	{
		std::complex<double> alpha = { i*h,0 };
		vector[i - 1] = back_integrand(alpha, t, eps);
		cout << i << " \n";
	}
	unsigned int end_time = clock(); // начальное врем€
	cout << end_time - start_time << endl;
	plotTheWaveField(h, { {"black", vector}/*, {"red", vector1}*/}, "xxx.txt");


	system("pause");
}

// квадратурна€ формула √аусса
std::complex<double> evaluate_integral(const std::function<std::complex<double>(double)>& integrand, double a, double b)
{
	std::complex<double> value = 0;
	const size_t size = sizeof(nodes) / sizeof(double);
	#pragma omp parallel for
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
	const boundary_value_problem<std::complex<double>> bvp = { // создаЄтс€ объект типа краева€ задача
		{
			[=](double x, const std::vector<std::complex<double>>& v) {return 1 / mu(x) * v[1]; }, // правые части уравнений системы
			[=](double x, const std::vector<std::complex<double>>& v) {return (rho(x) * s * s + alpha * alpha * mu(x)) * v[0]; }

		},{ { 0,0 } },{ { 1, 1.0 } } }; // краевые услови€*/
	return bvp.solve()[0] / s; // набор точек
}

//јналитическое выражение

std::complex<double> analytical(std::complex<double> alpha, std::complex<double> s, double t)
{
	const auto gamma = sqrt(alpha * alpha + s * s);
	const auto e = exp(gamma);
	const auto sh = (e - 1.0 / e) / 2.0;
	const auto ch = (e + 1.0 / e) / 2.0;
	return  /*exp(s * t) **/ (1.0 / s) * sh / (gamma * ch);
}

//”словные вычислени€

std::complex<double> back_integrands(std::complex<double> alpha, double t, double eps)
{
	int i = 0;
	const auto e = exp(alpha);
	const auto sh = (e - 1.0 / e) / 2.0;
	const auto ch = (e + 1.0 / e) / 2.0;
	std::complex<double> sm = sh / alpha / ch;	
	std::complex<double> value;
	do {
		if (i >= 1000) break;
		const auto tmp = (2.0 * i + 1) * pi / 2.0;
		const auto sqr = alpha * alpha + tmp * tmp;
		value = -2.0 * cos(sqrt(sqr) * t) / sqr;
		i++;
		sm += value;
	} while (value.real() * value.real() + value.imag() * value.imag() > eps * eps);
	return sm;
}

std::complex<double> contour(double s, double zeta)
{
	return{ zeta, s };
}

std::complex<double> new_integrand(std::complex<double> alpha, double s, double t)
{
	const std::complex<double> I = { 0,1 };
	const std::complex<double> point = contour(s, 0.5);
	return integrand(alpha, point, t) * exp(0.5) * (cos(s) + I * sin(s)) * I;
}

std::complex<double> back_integrand(std::complex<double> alpha, double t, double eps)
{
	std::complex<double> sum = { 0.0, 0.0 };
	const std::complex<double> denumerator = { 0, 2 * pi };
	int i = 1;
	std::complex<double> value = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, -pi / t, pi / t);
	while (value.real() * value.real() + value.imag() * value.imag() > eps * eps)
	{
		sum += value;
		if (i >= 1000) break;
		auto x = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, pi * i / t, pi * (i + 2) / t);
		auto y = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, -pi * (i + 2) / t, -pi * i / t);
		value = x + y;
		i+=2;
	}
	return sum / denumerator;
}