#pragma once
#include <complex>
#include <fstream>
#include <functional>
#include <vector>
#include <map>

void showVector(std::vector<double>& complexValuedVector,
                std::ostream& stream);


void showVector(const std::vector<double>& complexValuedVector,
                std::ostream& stream, double step, double min = 0);


/**
 * \brief �������� ������ �� ������
 * \param step ���
 * \param complexValuedVector ������
 * \param mapping �����������
 * \param stream �������� �����
 * \param a ����� ������� �������
 */
void showVector(double step,
    const std::vector<std::complex<double>>& complexValuedVector,
    const std::function<double(std::complex<double>)>& mapping, std::ostream& stream, double a = 0);


void addTheCurve(const std::vector<double>& fieldHomogeneous, std::ofstream& stream,
                 const std::string& color, double step, double min = 0);



void plotTheWaveField(const std::map<std::string, std::vector<double>>& vectors, const std::string& fileName,
                      double step, double min = 0);


/**
 * \brief  �������� ������ �� ������
 * \param step ���
 * \param fieldHomogeneous ������
 * \param stream �������� �����
 * \param color ����
 * \param a ����� �������
 */
void addTheCurve(double step, const std::vector<std::complex<double>>& fieldHomogeneous,
                 std::ofstream& stream, const std::string& color, double a = 0);

/**
 * \brief ��������� ���������
 * \param step ���
 * \param vectors ����� ��������
 * \param fileName ��� �����
 * \param a ����� �������
 */
void plotTheWaveField(double step, const std::map<std::string,   std::vector<std::complex<double>>>& vectors,
                      const std::string& fileName, double a = 0);

