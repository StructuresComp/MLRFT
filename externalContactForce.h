#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticRod &m_rod, timeStepper &m_stepper);
	~externalContactForce();

	void computeFc();
	void computeJc();
	void getContactForce(int i, Vector3d localForce);
	void setZeroForce();

private:
	elasticRod *rod;
	timeStepper *stepper;

	VectorXd Force;
};

#endif
