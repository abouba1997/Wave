// Wave.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "functions.h"

double p_0(double x, double y) {
	return exp(-(pow(x * 10., 2) /*+ pow(y * 10. - 5, 2)*/) / 0.8); // Source in center
	//return exp(-(pow(x * 10. - 5, 2) /*+ pow(y * 10. - 5, 2)*/) / 0.8); // Source in center
	//return exp(-(pow(x * 10. - 2, 2) + pow(y * 10. - 2, 2)) / 0.8) + exp(-(pow(x * 10. - 8, 2) + pow(y * 10. - 8, 2)) / 0.8); // Two sources
}

double f(double t, double x, double y) {
	return 0;
}

int main()
{
    int I = 100, J = 100, T = 1;
    double sigma_max = 0, C = 2;
    Wave2D wave(I, J, T, sigma_max, C, p_0, f);

	cout << endl << endl << endl;

	cout << "\t\tProgram completed!" << endl;

	cout << endl << endl << endl;
}
