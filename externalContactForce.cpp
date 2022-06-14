#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticRod &m_rod, timeStepper &m_stepper)
{
	rod = &m_rod;
	stepper = &m_stepper;

	Force.setZero(rod->ndof);
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::setZeroForce()
{
	Force.setZero(rod->ndof);
}

void externalContactForce::getContactForce(int i, Vector3d localForce)
{
	Force.segment(4 * i, 3) = Force.segment(4 * i, 3) + localForce;
}

void externalContactForce::computeFc()
{
	for (int i = 0; i < rod->ndof; i++)
	{
		stepper->addForce(i, -Force[i]);
	}
}

void externalContactForce::computeJc()
{
	;
}
