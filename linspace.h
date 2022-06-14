#ifndef LINSPACE_H
#define LINSPACE_H

#include "eigenIncludes.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Defining a function linspace which creates linearly spaced vectors
VectorXd linspace(double a, double b, double n)
{
	std::vector<double> v(n);

	double step;
	step = (b - a) / (n - 1);

	// assign values to the vector
	for (int i = 0; i < n; ++i)
	{
		v[i] = a + (i * step);
	}

	int N;
	N = int(n);
	VectorXd vec(N);

	vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());

	return vec;
}

#endif
