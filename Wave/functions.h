#pragma once
#include <iostream>
#include <functional>
#include <fstream>

using namespace std;

class Wave2D
{
public:
	Wave2D(int _I, int _J, double _T, double _sigma_max, double _C, function<double(double, double)> _p_0, function<double(double, double, double)> _f);
	~Wave2D();

private:
	int I, J, N;
	double T, C;

	double *x, *y;

	int PML_Size = 20;
	double sigma_max;

	double** u_curr;
	double** u_next;

	double** v_curr;
	double** v_next;

	double** p_x_curr;
	double** p_x_next;

	double** p_y_curr;
	double** p_y_next;

	double** p_curr;
	double** p_next;

	function<double(double, double)> p_0;
	function<double(double, double, double)> f;

	void writeFile(int i);
	void compute(double, double, int);
};