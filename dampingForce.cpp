#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper,
						   double m_viscosity, std::vector<double> m_Ct, std::vector<double> m_Cp, std::vector<double> m_Cz)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;
	Ct = m_Ct;
	Cp = m_Cp;
	Cz = m_Cz;
	Id3 << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	f.setZero(3);
	ForceVec = VectorXd::Zero(rod->ndof); // added
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	ForceVec = VectorXd::Zero(rod->ndof); // added

	for (int i = 2; i < rod->nv; i++)
	{

		u = rod->getVelocity(i);

		if (i > 2 && i < rod->nv - 1)
		{

			p1 = rod->getVertex(i - 1);
			p2 = rod->getVertex(i + 1);
			t = p2 - p1;
			t = t / t.norm();
		}
		else if (i == 2)
		{

			p1 = rod->getVertex(i);
			p2 = rod->getVertex(i + 1);
			t = p2 - p1;
			t = t / t.norm();
		}
		else
		{

			p1 = rod->getVertex(i - 1);
			p2 = rod->getVertex(i);
			t = p2 - p1;
			t = t / t.norm();
		}

		tmp_Ct = Ct[i - 2]; // tangential
		tmp_Cp = Cp[i - 2]; // normal
		tmp_Cz = Cz[i - 2]; // z-dir tan.crossnorm

		parvec = (u.dot(t)) * t;
		normvec = (u - (u.dot(t)) * t);
		lastvec = normvec.cross(t);

		if ((normvec.norm()) == 0)
		{
			if (i == 2)
			{
				f = -(tmp_Ct * t + tmp_Cp * normvec + tmp_Cz * lastvec) * u.norm() * viscosity * rod->voronoiLen(3) / 2;
			}
			else
			{
				f = -(tmp_Ct * t + tmp_Cp * normvec + tmp_Cz * lastvec) * u.norm() * viscosity * rod->voronoiLen(i);
			}
		}
		else
		{
			if (i == 2)
			{
				f = -(tmp_Ct * t + tmp_Cp * normvec / normvec.norm() + tmp_Cz * lastvec / lastvec.norm()) * u.norm() * viscosity * rod->voronoiLen(3) / 2;
			}
			else
			{
				f = -(tmp_Ct * t + tmp_Cp * normvec / normvec.norm() + tmp_Cz * lastvec / lastvec.norm()) * u.norm() * viscosity * rod->voronoiLen(i);
			}
		}
		// cout << f;
		for (int k = 0; k < 3; k++)
		{
			ind = 4 * i + k;

			ForceVec(ind) = ForceVec(ind) + f(k);
			stepper->addForce(ind, -f[k]); // subtracting external force
		}
	}
}

void dampingForce::computeJd()
{
	// Remember that dF/dx = 1/dt * dF/dv

	for (int i = 2; i < rod->nv; i++)
	{
		u = rod->getVelocity(i);
		ux = u(0);
		uy = u(1);
		uz = u(2);
		if (i > 2 && i < rod->nv - 1)
		{

			p1 = rod->getVertex(i - 1);
			xm1 = p1(0);
			ym1 = p1(1);
			zm1 = p1(2);
			p2 = rod->getVertex(i + 1);
			xp1 = p2(0);
			yp1 = p2(1);
			zp1 = p2(2);

			t = p2 - p1;
			t = t / t.norm();
			dl = t.norm();
		}
		else if (i == 2)
		{

			p1 = rod->getVertex(i);
			x = p1(0);
			y = p1(1);
			z = p1(2);
			p2 = rod->getVertex(i + 1);
			xp1 = p2(0);
			yp1 = p2(1);
			zp1 = p2(2);
			t = p2 - p1;
			t = t / t.norm();
			dl = t.norm();
		}
		else
		{

			p1 = rod->getVertex(i - 1);
			xm1 = p1(0);
			ym1 = p1(1);
			zm1 = p1(2);
			p2 = rod->getVertex(i);
			x = p2(0);
			y = p2(1);
			z = p2(2);
			t = p2 - p1;
			t = t / t.norm();
			dl = t.norm();
		}

		parvec = (u.dot(t)) * t;
		normvec = (u - (u.dot(t)) * t);
		lastvec = normvec.cross(t);

		tmp_Ct = Ct[i - 2]; // tangential
		tmp_Cp = Cp[i - 2]; // normal
		tmp_Cz = Cz[i - 2]; // z-dir tan.crossnorm

		if (i == 2)
		{
			jac(0, 0) = tmp_Cp * (1 - (x - xp1) * (x - xp1) / dl * dl) + (tmp_Ct * (x - xp1) * (x - xp1)) / dl * dl + (tmp_Cz * (x - xp1) * (uz * (-y + yp1) + uy * (z - zp1))) / dl * dl;
			jac(1, 0) = -((tmp_Cp * (x - xp1) * (y - yp1)) / dl * dl) + (tmp_Ct * (x - xp1) * (y - yp1)) / dl * dl + (tmp_Cz * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1)) * (-z + zp1)) / dl * dl + (tmp_Cz * (x - xp1) * (uz * (x - xp1) + ux * (-z + zp1))) / dl * dl;
			jac(2, 0) = (tmp_Cz * (x - xp1) * (uy * (-x + xp1) + ux * (y - yp1))) / dl * dl + (tmp_Cz * (y - yp1) * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1))) / dl * dl - (tmp_Cp * (x - xp1) * (z - zp1)) / dl * dl + (tmp_Ct * (x - xp1) * (z - zp1)) / dl * dl;
			jac(0, 1) = -((tmp_Cp * (x - xp1) * (y - yp1)) / dl * dl) + (tmp_Ct * (x - xp1) * (y - yp1)) / dl * dl + (tmp_Cz * (y - yp1) * (uz * (-y + yp1) + uy * (z - zp1))) / dl * dl + (tmp_Cz * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1)) * (z - zp1)) / dl * dl;
			jac(1, 1) = tmp_Cp * (1 - (y - yp1) * (y - yp1) / dl * dl) + (tmp_Ct * (y - yp1) * (y - yp1)) / dl * dl + (tmp_Cz * (y - yp1) * (uz * (x - xp1) + ux * (-z + zp1))) / dl * dl;
			jac(2, 1) = (tmp_Cz * (uy * (-x + xp1) + ux * (y - yp1)) * (y - yp1)) / dl * dl + (tmp_Cz * (-x + xp1) * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1))) / dl * dl - (tmp_Cp * (y - yp1) * (z - zp1)) / dl * dl + (tmp_Ct * (y - yp1) * (z - zp1)) / dl * dl;
			jac(0, 2) = (tmp_Cz * (-y + yp1) * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1))) / dl * dl - (tmp_Cp * (x - xp1) * (z - zp1)) / dl * dl + (tmp_Ct * (x - xp1) * (z - zp1)) / dl * dl + (tmp_Cz * (uz * (-y + yp1) + uy * (z - zp1)) * (z - zp1)) / dl * dl;
			jac(1, 2) = (tmp_Cz * (x - xp1) * (ux * (x - xp1) + uy * (y - yp1) + uz * (z - zp1))) / dl * dl - (tmp_Cp * (y - yp1) * (z - zp1)) / dl * dl + (tmp_Ct * (y - yp1) * (z - zp1)) / dl * dl + (tmp_Cz * (z - zp1) * (uz * (x - xp1) + ux * (-z + zp1))) / dl * dl;
			jac(2, 2) = tmp_Cp * (1 - (z - zp1) * (z - zp1) / dl * dl) + (tmp_Cz * (uy * (-x + xp1) + ux * (y - yp1)) * (z - zp1)) / dl * dl + (tmp_Ct * (z - zp1) * (z - zp1)) / dl * dl;
		}
		else if (i > 2 && i < rod->nv - 1)
		{
			jac(0, 0) = tmp_Cp * (1 - (xm1 - xp1) * (xm1 - xp1) / (4 * dl * dl)) + (tmp_Ct * (xm1 - xp1) * (xm1 - xp1)) / (4 * dl * dl) + (tmp_Cz * (xm1 - xp1) * (uz * (-ym1 + yp1) + uy * (zm1 - zp1))) / (4 * dl * dl);
			jac(1, 0) = -1 / 4 * (tmp_Cp * (xm1 - xp1) * (ym1 - yp1)) / dl * dl + (tmp_Ct * (xm1 - xp1) * (ym1 - yp1)) / (4 * dl * dl) + (tmp_Cz * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1)) * (-zm1 + zp1)) / (4 * dl * dl) + (tmp_Cz * (xm1 - xp1) * (uz * (xm1 - xp1) + ux * (-zm1 + zp1))) / (4 * dl * dl);
			jac(2, 0) = (tmp_Cz * (xm1 - xp1) * (uy * (-xm1 + xp1) + ux * (ym1 - yp1))) / (4 * dl * dl) + (tmp_Cz * (ym1 - yp1) * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1))) / (4 * dl * dl) - (tmp_Cp * (xm1 - xp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Ct * (xm1 - xp1) * (zm1 - zp1)) / (4 * dl * dl);
			jac(0, 1) = -1 / 4 * (tmp_Cp * (xm1 - xp1) * (ym1 - yp1)) / dl * dl + (tmp_Ct * (xm1 - xp1) * (ym1 - yp1)) / (4 * dl * dl) + (tmp_Cz * (ym1 - yp1) * (uz * (-ym1 + yp1) + uy * (zm1 - zp1))) / (4 * dl * dl) + (tmp_Cz * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1)) * (zm1 - zp1)) / (4 * dl * dl);
			jac(1, 1) = tmp_Cp * (1 - (ym1 - yp1) * (ym1 - yp1) / (4 * dl * dl)) + (tmp_Ct * (ym1 - yp1) * (ym1 - yp1)) / (4 * dl * dl) + (tmp_Cz * (ym1 - yp1) * (uz * (xm1 - xp1) + ux * (-zm1 + zp1))) / (4 * dl * dl);
			jac(2, 1) = (tmp_Cz * (uy * (-xm1 + xp1) + ux * (ym1 - yp1)) * (ym1 - yp1)) / (4 * dl * dl) + (tmp_Cz * (-xm1 + xp1) * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1))) / (4 * dl * dl) - (tmp_Cp * (ym1 - yp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Ct * (ym1 - yp1) * (zm1 - zp1)) / (4 * dl * dl);
			jac(0, 2) = (tmp_Cz * (-ym1 + yp1) * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1))) / (4 * dl * dl) - (tmp_Cp * (xm1 - xp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Ct * (xm1 - xp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Cz * (uz * (-ym1 + yp1) + uy * (zm1 - zp1)) * (zm1 - zp1)) / (4 * dl * dl);
			jac(1, 2) = (tmp_Cz * (xm1 - xp1) * (ux * (xm1 - xp1) + uy * (ym1 - yp1) + uz * (zm1 - zp1))) / (4 * dl * dl) - (tmp_Cp * (ym1 - yp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Ct * (ym1 - yp1) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Cz * (zm1 - zp1) * (uz * (xm1 - xp1) + ux * (-zm1 + zp1))) / (4 * dl * dl);
			jac(2, 2) = tmp_Cp * (1 - (zm1 - zp1) * (zm1 - zp1) / (4 * dl * dl)) + (tmp_Cz * (uy * (-xm1 + xp1) + ux * (ym1 - yp1)) * (zm1 - zp1)) / (4 * dl * dl) + (tmp_Ct * (zm1 - zp1) * (zm1 - zp1)) / (4 * dl * dl);
		}
		else
		{
			jac(0, 0) = (tmp_Cp * (dl * dl - (x - xm1) * (x - xm1)) + tmp_Ct * (x - xm1) * (x - xm1) + tmp_Cz * (x - xm1) * (uz * (-y + ym1) + uy * (z - zm1))) / dl * dl;
			jac(1, 0) = (-(tmp_Cp * (x - xm1) * (y - ym1)) + tmp_Ct * (x - xm1) * (y - ym1) + tmp_Cz * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) * (-z + zm1) + tmp_Cz * (x - xm1) * (uz * (x - xm1) + ux * (-z + zm1))) / dl * dl;
			jac(2, 0) = (tmp_Cz * (x - xm1) * (uy * (-x + xm1) + ux * (y - ym1)) + tmp_Cz * (y - ym1) * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) - tmp_Cp * (x - xm1) * (z - zm1) + tmp_Ct * (x - xm1) * (z - zm1)) / dl * dl;
			jac(0, 1) = (-(tmp_Cp * (x - xm1) * (y - ym1)) + tmp_Ct * (x - xm1) * (y - ym1) + tmp_Cz * (y - ym1) * (uz * (-y + ym1) + uy * (z - zm1)) + tmp_Cz * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) * (z - zm1)) / dl * dl;
			jac(1, 1) = (tmp_Cp * (dl * dl - (y - ym1) * (y - ym1)) + tmp_Ct * (y - ym1) * (y - ym1) + tmp_Cz * (y - ym1) * (uz * (x - xm1) + ux * (-z + zm1))) / dl * dl;
			jac(2, 1) = (tmp_Cz * (uy * (-x + xm1) + ux * (y - ym1)) * (y - ym1) + tmp_Cz * (-x + xm1) * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) - tmp_Cp * (y - ym1) * (z - zm1) + tmp_Ct * (y - ym1) * (z - zm1)) / dl * dl;
			jac(0, 2) = (tmp_Cz * (-y + ym1) * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) - tmp_Cp * (x - xm1) * (z - zm1) + tmp_Ct * (x - xm1) * (z - zm1) + tmp_Cz * (uz * (-y + ym1) + uy * (z - zm1)) * (z - zm1)) / dl * dl;
			jac(1, 2) = (tmp_Cz * (x - xm1) * (ux * (x - xm1) + uy * (y - ym1) + uz * (z - zm1)) - tmp_Cp * (y - ym1) * (z - zm1) + tmp_Ct * (y - ym1) * (z - zm1) + tmp_Cz * (z - zm1) * (uz * (x - xm1) + ux * (-z + zm1))) / dl * dl;
			jac(2, 2) = (tmp_Cp * (dl * dl - (z - zm1) * (z - zm1)) + tmp_Cz * (uy * (-x + xm1) + ux * (y - ym1)) * (z - zm1) + tmp_Ct * (z - zm1) * (z - zm1)) / dl * dl;
		}

		// Explicit setup

		//  jac(0, 0) = 0;
		//  jac(1, 0) = 0;
		//  jac(2, 0) = 0;
		//  jac(0, 1) = 0;
		//  jac(1, 1) = 0;
		//  jac(0, 2) = 0;
		//  jac(1, 2) = 0;
		//  jac(2, 2) = 0;

		for (int kx = 0; kx < 3; kx++)
		{
			indx = 4 * i + kx;
			for (int ky = 0; ky < 3; ky++)
			{
				indy = 4 * i + ky;
				stepper->addJacobian(indx, indy, +jac(kx, ky)); // subtracting external force
			}
		}
	}
}
