#include "functions.h"

Wave2D::Wave2D(int _I, int _J, double _T, double _sigma_max, double _C, function<double(double, double)> _p_0, function<double(double, double, double)> _f)
{
	I = _I; J = _J; T = _T; sigma_max = _sigma_max; C = _C; N = 0; p_0 = _p_0; f = _f;

	double hx = 1. / (I - 1.), hy = 1. / (J - 1.);
	double t = 0.0001;

	if (hx != hy) {
		throw std::exception("Parameters erros");
	}

	if (sigma_max > 2. / t) {
		throw exception("Absorbing coef error");
	}

	if (t <= 2 * hx / C) {
		N = int(T / t);
	}

	x = new double[I + 2 * PML_Size];
	y = new double[I + 2 * PML_Size];

	u_curr = new double* [I + 2 * PML_Size];
	u_next = new double* [I + 2 * PML_Size];

	v_curr = new double* [J + 2 * PML_Size];
	v_next = new double* [J + 2 * PML_Size];

	p_x_curr = new double* [I + 2 * PML_Size];
	p_x_next = new double* [I + 2 * PML_Size];

	p_y_curr = new double* [J + 2 * PML_Size];
	p_y_next = new double* [J + 2 * PML_Size];

	p_curr = new double* [I + 2 * PML_Size];
	p_next = new double* [I + 2 * PML_Size];

	for (size_t i = 0; i < I + 2 * PML_Size; i++)
	{
		u_curr[i] = new double[I + 2 * PML_Size];
		u_next[i] = new double[I + 2 * PML_Size];

		v_curr[i] = new double[J + 2 * PML_Size];
		v_next[i] = new double[J + 2 * PML_Size];

		p_x_curr[i] = new double[I + 2 * PML_Size];
		p_x_next[i] = new double[I + 2 * PML_Size];

		p_y_curr[i] = new double[J + 2 * PML_Size];
		p_y_next[i] = new double[J + 2 * PML_Size];

		p_curr[i] = new double[I + 2 * PML_Size];
		p_next[i] = new double[I + 2 * PML_Size];
	}

	for (int i = 0; i < I + 2 * PML_Size; i++)
	{
		x[i] = ((double)i - (double)PML_Size + 1.) * hx;
		y[i] = ((double)i - (double)PML_Size + 1.) * hy;
	}

	for (size_t i = 0; i < I + 2 * PML_Size; i++)
	{
		for (size_t j = 0; j < J + 2 * PML_Size; j++)
		{
			u_curr[i][j] = 0;
			u_next[i][j] = 0;

			v_curr[i][j] = 0;
			v_next[i][j] = 0;

			p_x_curr[i][j] = 0;
			p_x_next[i][j] = 0;

			p_y_curr[i][j] = 0;
			p_y_next[i][j] = 0;

			p_curr[i][j] = 0;
			p_next[i][j] = 0;
		}
	}

	for (size_t i = PML_Size - 1; i < I + PML_Size - 1; i++)
	{
		for (size_t j = PML_Size - 1; j < J + PML_Size - 1; j++)
		{
			p_curr[i][j] = p_0(x[i], y[j]);

			p_x_curr[i][j] = p_curr[i][j] / 2;
			p_y_curr[i][j] = p_curr[i][j] / 2;
		}
	}

	// Reflecting condition

	for (size_t j = PML_Size - 1; j < I + PML_Size - 1; j++)
	{
		p_curr[PML_Size - 1][j] = 0;
		p_curr[I + PML_Size - 2][j] = 0;
	}

	for (size_t i = PML_Size - 1; i < I + PML_Size - 1; i++)
	{
		p_curr[i][PML_Size - 1] = 0;
		p_curr[i][I + PML_Size - 2] = 0;
	}


	std::ofstream param("param.txt");
	param << "N = " << N << "\n";
	param << "I = " << I << "\n";
	param << "J = " << J << "\n";
	param << "T = " << T << "\tt = " << t << "\n";
	param.close();

	writeFile(0);

	compute(hx, t, 1);
}

Wave2D::~Wave2D()
{
	delete[] x;
	delete[] y;

	for (size_t i = 0; i < I + 2 * PML_Size; i++)
	{
		delete[] u_curr[i];
		delete[] u_next[i];

		delete[] v_curr[i];
		delete[] v_next[i];

		delete[] p_x_curr[i];
		delete[] p_x_next[i];

		delete[] p_y_curr[i];
		delete[] p_y_next[i];

		delete[] p_curr[i];
		delete[] p_next[i];
	}

	delete[] u_curr;
	delete[] u_next;

	delete[] v_curr;
	delete[] v_next;

	delete[] p_x_curr;
	delete[] p_x_next;

	delete[] p_y_curr;
	delete[] p_y_next;

	delete[] p_curr;
	delete[] p_next;
}

void Wave2D::writeFile(int k)
{
	char fn[30];
	sprintf_s(fn, ".u.dat.%03d", k);
	ofstream datafile(fn);
	for (int i = PML_Size - 1; i < I + PML_Size - 1; i++) {
		for (int j = PML_Size - 1; j < J + PML_Size - 1; j++) {
			datafile << x[i] << "\t" << y[j] << "\t" << p_curr[i][j] << endl;
		}
		datafile << endl;
	}
	datafile.close();
}

void Wave2D::compute(double h, double t, int k)
{
	double C_1 = C * t / (2 * h);
	if (sigma_max == 0) {
		for (size_t n = 0; n < N; n++)
		{
			for (size_t j = PML_Size - 1; j < J + PML_Size - 1; j++)
			{
				//u_next[PML_Size - 1][j] = 0;
				for (size_t i = PML_Size; i < I + PML_Size - 2; i++)
				{
					u_next[i][j] = C_1 * (p_curr[i + 1][j] - p_curr[i - 1][j]) + u_curr[i][j];
				}
				//u_next[I + PML_Size - 2][j] = 0;
			}

			for (size_t i = PML_Size - 1; i < I + PML_Size - 1; i++)
			{
				//v_next[i][PML_Size - 1] = 0;
				for (size_t j = PML_Size; j < J + PML_Size - 2; j++)
				{
					v_next[i][j] = C_1 * (p_curr[i][j + 1] - p_curr[i][j - 1]) + v_curr[i][j];
				}
				//v_next[i][I + PML_Size - 2] = 0;
			}

			for (size_t j = PML_Size - 1; j < J + PML_Size - 1; j++)
			{
				p_x_next[PML_Size - 1][j] = 0;
				for (size_t i = PML_Size; i < I + PML_Size - 2; i++)
				{
					p_x_next[i][j] = C_1 * (u_curr[i + 1][j] - u_curr[i - 1][j]) + p_x_curr[i][j];
				}
				p_x_next[I + PML_Size - 2][j] = 0;
			}

			for (size_t i = PML_Size - 1; i < I + PML_Size - 1; i++)
			{
				p_y_next[i][PML_Size - 1];
				for (size_t j = PML_Size; j < J + PML_Size - 2; j++)
				{
					p_y_next[i][j] = C_1 * (v_curr[i][j + 1] - v_curr[i][j - 1]) + p_y_curr[i][j];
				}
				p_y_next[i][I + PML_Size - 2];
			}

			///////////////// Adding P
			for (int i = PML_Size - 1; i < I + PML_Size - 1; i++)
			{
				for (int j = PML_Size - 1; j < J + PML_Size - 1; j++)
				{
					p_next[i][j] = p_x_next[i][j] + p_y_next[i][j];
				}
			}

			// Reflecting condition

			for (size_t j = PML_Size - 1; j < I + PML_Size - 1; j++)
			{
				p_next[PML_Size - 1][j] = 0;
				p_next[ I + PML_Size - 2][j] = 0;
			}
			
			for (size_t i = PML_Size - 1; i < I + PML_Size - 1; i++)
			{
				p_next[i][PML_Size - 1] = 0;
				p_next[i][I + PML_Size - 2] = 0;
			}

			if (n % 100 == 0) {
				writeFile(k);
				k++;
			}

			for (int i = PML_Size - 1; i < I + PML_Size - 1; i++)
			{
				for (int j = PML_Size - 1; j < J + PML_Size - 1; j++)
				{
					p_curr[i][j] = p_next[i][j];

					p_x_curr[i][j] = p_x_next[i][j];
					p_y_curr[i][j] = p_y_next[i][j];

					u_curr[i][j] = u_next[i][j];
					v_curr[i][j] = v_next[i][j];
				}
			}
		}
	}
}
