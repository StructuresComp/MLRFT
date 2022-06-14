#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity, std::vector<double> m_Ct, std::vector<double> m_Cp, std::vector<double> m_Cz);
	~dampingForce();
	void computeFd();
	void computeJd();

	VectorXd ForceVec; // added

private:
	elasticRod *rod;
	timeStepper *stepper;
	double viscosity, eta_per, eta_par, coeffa, x, y, z, xm1, ym1, zm1, xp1, yp1, zp1, dl, ux, uy, uz;
	std::vector<double> Ct, Cp, Cz;
	Vector3d t, u, f, normvec, parvec, lastvec, p1, p2;
	double tmp_Ct, tmp_Cp, tmp_Cz;
	int ind, indx, indy;
	Matrix3d Id3, jac;
	double dt;
};

#endif
