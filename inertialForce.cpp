#include "inertialForce.h"

inertialForce::inertialForce(elasticRod &m_rod, timeStepper &m_stepper)
{
	rod = &m_rod;
	stepper = &m_stepper;

	ForceVec = VectorXd::Zero(rod->ndof);
}

inertialForce::~inertialForce()
{
	;
}

void inertialForce::computeFi()
{

	for (int i = 0; i < rod->ndof; i++)
	{
		f = rod->massArray[i] * (rod->x[i] - rod->x0[i]) / ((rod->dt) * (rod->dt)) - (rod->massArray[i] * rod->u[i]) / (rod->dt);

		ForceVec(i) = f;

		stepper->addForce(i, f);
	}
}

void inertialForce::computeJi()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		jac = rod->massArray(i) / ((rod->dt) * (rod->dt));
		stepper->addJacobian(i, i, jac);
	}
}
