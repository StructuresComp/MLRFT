#include "RegularizedStokeslet.h"

RegularizedStokeslet::RegularizedStokeslet(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity, double m_epsilon)
{
	rod = &m_rod;
	stepper = &m_stepper;

	viscosity = m_viscosity;

	epsilon = m_epsilon;

	A = MatrixXd::Zero(3 * rod->nv, 3 * rod->nv);

	VelocityVec = VectorXd::Zero(3 * rod->nv);
	ViscousForce = VectorXd::Zero(3 * rod->nv);

	Id3 << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	ForceVec = VectorXd::Zero(rod->ndof);
}

RegularizedStokeslet::~RegularizedStokeslet()
{
	;
}

void RegularizedStokeslet::prepareForViscousForce()
{
	A = MatrixXd::Zero(3 * rod->nv, 3 * rod->nv);

	VelocityVec = VectorXd::Zero(3 * rod->nv);
	ViscousForce = VectorXd::Zero(3 * rod->nv);

	// loop for first rod
	for (int i = 2; i < rod->nv; i++)
	{
		uPos = rod->getVertexOld(i);
		uVelocity = rod->getVelocityOld(i);

		VelocityVec(3 * i) = uVelocity(0);
		VelocityVec(3 * i + 1) = uVelocity(1);
		VelocityVec(3 * i + 2) = uVelocity(2);

		for (int j = 2; j < rod->ne; j++)
		{

			y_0 = rod->getVertexOld(j);
			y_1 = rod->getVertexOld(j + 1);

			x_0 = uPos - y_0;
			x_1 = uPos - y_1;

			vDirection = y_0 - y_1;
			edgeLength = vDirection.norm();

			R_0 = sqrt(x_0.norm() * x_0.norm() + epsilon * epsilon);
			R_1 = sqrt(x_1.norm() * x_1.norm() + epsilon * epsilon);
			computeTnumber();

			M2 = ((T1_1 + epsilon * epsilon * T1_3) * Id3 + T1_3 * (x_0 * x_0.transpose()) + T2_3 * (x_0 * vDirection.transpose() + vDirection * x_0.transpose()) + T3_3 * (vDirection * vDirection.transpose()));
			M1 = ((T0_1 + epsilon * epsilon * T0_3) * Id3 + T0_3 * (x_0 * x_0.transpose()) + T1_3 * (x_0 * vDirection.transpose() + vDirection * x_0.transpose()) + T2_3 * (vDirection * vDirection.transpose())) - M2;

			A.block(3 * i, 3 * j, 3, 3) = A.block(3 * i, 3 * j, 3, 3) + M1;
			A.block(3 * i, 3 * (j + 1), 3, 3) = A.block(3 * i, 3 * (j + 1), 3, 3) + M2;
		}
	}
	
	ViscousForce = -A.ldlt().solve(8 * M_PI * viscosity * VelocityVec);

}

void RegularizedStokeslet::computeFrs()
{
	ForceVec = VectorXd::Zero(rod->ndof);

	for (int i = 2; i < rod->nv; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			stepper->addForce(4 * i + j, -ViscousForce(3 * i + j));

			ForceVec(4 * i + j) = ForceVec(4 * i + j) + ViscousForce(3 * i + j);
		}
	}
}

void RegularizedStokeslet::computeJrs()
{
	;
}

void RegularizedStokeslet::computeTnumber()
{
	T0_1 = (log(edgeLength * R_1 + (x_0.dot(vDirection) + edgeLength * edgeLength)) - log(edgeLength * R_0 + x_0.dot(vDirection))) / edgeLength;
	T0_3 = (-1 / (R_1 * (edgeLength * R_1 + (x_0.dot(vDirection) + edgeLength * edgeLength)))) - (-1 / (R_0 * (edgeLength * R_0 + (x_0.dot(vDirection)))));
	T1_1 = (R_1 / (edgeLength * edgeLength) - R_0 / (edgeLength * edgeLength)) - T0_1 * x_0.dot(vDirection) / (edgeLength * edgeLength);
	T1_3 = -(1 / (R_1 * edgeLength * edgeLength) - 1 / (R_0 * edgeLength * edgeLength)) - T0_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
	T2_3 = -(1 / (R_1 * edgeLength * edgeLength) - 0 / (R_0 * edgeLength * edgeLength)) + T0_1 / (edgeLength * edgeLength) - T1_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
	T3_3 = -(1 / (R_1 * edgeLength * edgeLength) - 0 / (R_0 * edgeLength * edgeLength)) + 2 * T1_1 / (edgeLength * edgeLength) - T2_3 * x_0.dot(vDirection) / (edgeLength * edgeLength);
}
