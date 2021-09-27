#pragma once

#include<complex>
#include <functional>
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
 * \param integrand подынтегральное выражение
 * \param a левая граница отрезка
 * \param b правая граница отрезка
 * \return значение интеграла
 */
std::complex<double> evaluate_integral(
	const std::function<std::complex<double>(double)>& integrand, double a = -1,
	double b = 1);

/**
 * \brief подынтегральное выражение
 * \param alpha параметр преобразования Лапласа
 * \param s параметр преобразования Фурье
 * \param t время
 * \return подынтегральное выражение
 */
std::complex<double> integrand(std::complex<double> alpha,
	std::complex<double> s, double t);

/**
 * \brief контур интегрирования
 * \param s параметр преобразования Лапласа
 * \param zeta абсцисса контура интегрирования
 * \return точка контура
 */
std::complex<double> contour(double s, double zeta = 0);

/**
 * \brief подынтегральное выражение вдольконтура интегрирования
 * \param alpha параметр преобразования Фурье
 * \param s параметр преобразования Лапласа
 * \param t время
 * \return значение
 */
std::complex<double> new_integrand(std::complex<double> alpha, double s, double t);

/**
 * \brief обратное преобразование Лапласа
 * \param alpha параметр преобразования Фурье
 * \param t время
 * \param eps точность вычислений
 * \return значение обратной трансформанты
 */
std::complex<double> back_integrand(std::complex<double> alpha, double t,
	double eps);

/**
 * \brief обратное преобразование Лапласа в однородном случае (аналитическое
 * выражение)
 * \param alpha параметр преобразования Фурье
 * \param t время
 * \param eps точность вычислений
 * \return значение
 */
std::complex<double> back_integrands(std::complex<double> alpha, double t,
	double eps);

std::complex<double> analytical(std::complex<double> alpha, std::complex<double> s, double t);
