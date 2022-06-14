#include "ResistiveForce.h"
#include <iostream>

resistiveForce::resistiveForce(elasticRod &m_rod, timeStepper &m_stepper,
							   double m_viscosity, double m_eta_per, double m_eta_par)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;
	eta_per = m_eta_per;
	eta_par = m_eta_par;

	Id3 << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	f.setZero(3);
	ForceVec = VectorXd::Zero(rod->ndof); // added
}

resistiveForce::~resistiveForce()
{
	;
}

void resistiveForce::computeFrft()
{
	ForceVec = VectorXd::Zero(rod->ndof); // added

	for (int i = 2; i < rod->ne; i++)
	{
		u = rod->getVelocity(i);
		t = rod->getTangent(i);
		f = -((eta_par - eta_per) * (u.dot(t)) * t + eta_per * u) * rod->voronoiLen(i);

		for (int k = 0; k < 3; k++)
		{
			ind = 4 * i + k;

			ForceVec(ind) = ForceVec(ind) + f[k]; // added
			stepper->addForce(ind, -f[k]);		  // subtracting external force
		}
	}
}

void resistiveForce::computeJrft()
{
	// Remember that dF/dx = 1/dt * dF/dv

	for (int i = 2; i < rod->ne; i++)
	{
		u = rod->getVelocity(i);
		t = rod->getTangent(i);

		{
			jac = -((eta_par - eta_per) * t * t.transpose() + eta_per * Id3) * rod->voronoiLen(i) / dt;
		}

		for (int kx = 0; kx < 3; kx++)
		{
			indx = 4 * i + kx;
			for (int ky = 0; ky < 3; ky++)
			{
				indy = 4 * i + ky;
				stepper->addJacobian(indx, indy, -jac(kx, ky)); // subtracting external force
			}
		}
	}
}
