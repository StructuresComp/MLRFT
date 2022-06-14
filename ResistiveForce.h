#ifndef RESISTIVEFORCE_H
#define RESISTIVEFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class resistiveForce
{
public:
	resistiveForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity, double m_eta_per, double m_eta_par);
	~resistiveForce();
	void computeFrft();
	void computeJrft();

	VectorXd ForceVec; // added

private:
	elasticRod *rod;
	timeStepper *stepper;
	double viscosity, eta_per, eta_par;
	Vector3d t, u, f;
	int ind, indx, indy;
	Matrix3d Id3, jac;
	double dt;
};

#endif
