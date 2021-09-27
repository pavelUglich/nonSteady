#include "integrals.h"
#include <complex>
#include "boundary_value_problem.h"

std::complex<double> evaluate_integral(const std::function<std::complex<double>(double)>& integrand, double a, double b)
{
	std::complex<double> value = 0;
	const size_t size = sizeof(nodes) / sizeof(double);
    #pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		const double x = (a + b) / 2 + (b - a) / 2 * nodes[i];
		value += weights[i] * integrand(x);
	}
	return (b - a) * value / 2.0;
}


std::complex<double> integrand(std::complex<double> alpha, std::complex<double> s, double t)
{
	const std::function<double(double)> mu = [](double x) {return 1.0; };
	const std::function<double(double)> rho = [](double x) {return 1.0; };
	const boundary_value_problem<std::complex<double>> bvp = { // создаётся объект типа краевая задача
		{
			[=](double x, const std::vector<std::complex<double>>& v) {return 1 / mu(x) * v[1]; }, // правые части уравнений системы
			[=](double x, const std::vector<std::complex<double>>& v) {return (rho(x) * s * s + alpha * alpha * mu(x)) * v[0]; }

		},{ { 0,0 } },{ { 1, 1.0 } } }; // краевые условия
	return bvp.solve()[0] / s; // набор точек
}

std::complex<double> contour(double s, double zeta)
{
		return{ zeta, s };
}

std::complex<double> new_integrand(std::complex<double> alpha, double s, double t)
{
	const std::complex<double> I = { 0,1 };
	const std::complex<double> point = contour(s, 0.5);
	return integrand(alpha, point, t) * exp(0.5*t) * (cos(s*t) + I * sin(s*t)) * I;
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
		if (i >= 10000) break;
		auto x = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, pi * i / t, pi * (i + 2) / t);
		auto y = evaluate_integral([=](double s) {return new_integrand(alpha, s, t); }, -pi * (i + 2) / t, -pi * i / t);
		value = x + y;
		i += 2;
	}
	return sum / denumerator;
}

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

std::complex<double> analytical(std::complex<double> alpha, std::complex<double> s, double t)
{
	const auto gamma = sqrt(alpha * alpha + s * s);
	const auto e = exp(gamma);
	const auto sh = (e - 1.0 / e) / 2.0;
	const auto ch = (e + 1.0 / e) / 2.0;
	return  /*exp(s * t) **/ (1.0 / s) * sh / (gamma * ch);
}

