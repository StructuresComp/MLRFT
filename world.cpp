#include "world.h"
#include "linspace.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");			// boolean
	saveData = m_inputData.GetBoolOpt("saveData");		// boolean
	inputName = m_inputData.GetStringOpt("input-file"); // string

	// Physical parameters
	helixpitchRatio = m_inputData.GetScalarOpt("helixpitchRatio"); // multiplied by helix radius m
	rodRadiusRatio = m_inputData.GetScalarOpt("rodRadiusRatio");   // multiplied by helix radius m
	contourRatio = m_inputData.GetScalarOpt("contourRatio");	   // multiplied by helix radius m
	helixradius = m_inputData.GetScalarOpt("helixradius");		   // meter
	gVector = m_inputData.GetVecOpt("gVector");					   // m/s^2
	maxIter = m_inputData.GetIntOpt("maxIter");					   // maximum number of iterations
	maxIterContact = m_inputData.GetIntOpt("maxIterContact");	   // maximum number of iterations

	youngM = m_inputData.GetScalarOpt("youngM");		 // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");		 // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");	 // seconds
	totalTime = m_inputData.GetScalarOpt("totalTime");	 // seconds
	tol = m_inputData.GetScalarOpt("tol");				 // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");			 // small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");		 // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");	 // viscosity in Pa-s
	translation = m_inputData.GetBoolOpt("translation"); // MLRFT translation model? True means translation/false means rotation.
	headSize = m_inputData.GetScalarOpt("headSize");

	useRSS = m_inputData.GetBoolOpt("use-RSS"); // If true, RSS will be used. RFT needs to be false
	useRFT = m_inputData.GetBoolOpt("use-RFT"); // If true, RFT will be used. RSS needs to be false, when both false, MLRFT active.
	includeContact = m_inputData.GetBoolOpt("include-contact");

	// geometry of helix

	helixpitch = helixpitchRatio * helixradius; // meter
	rodRadius = helixradius / rodRadiusRatio;
	flagellaLength = contourRatio * helixradius;
	deltaLengthInput = 5 * rodRadius;
	axisLengthInput = flagellaLength / sqrt(((2 * M_PI * helixradius) * (2 * M_PI * helixradius)) + (helixpitch * helixpitch)) * helixpitch;

	double newRodLength = (flagellaLength * 1) + helixradius + helixradius; // RODLENGTH CORRECTED
	RodLength = newRodLength;

	// RFT/RSS coefficients and regularizer
	epsilon = 1.031 * rodRadius;
	eta_per = 4.0 * M_PI * viscosity / (log(2.0 * helixpitch / rodRadius) + 0.5);
	eta_par = 2.0 * M_PI * viscosity / (log(2.0 * helixpitch / rodRadius) - 0.5);

	int newNe = round(flagellaLength / deltaLengthInput);

	numVertices = newNe;
	int totVertices = numVertices + 2;

	shearM = youngM / (2.0 * (1.0 + Poisson)); // shear modulus

	// Read input file to get angular/translational velocity, if translation == 1 input-file needs to include the speed in m/s. if translation = 0 rpm.
	ReadOmegaData();

	smallDist2 = pow(rodRadius / 1000.0, 2);

	minDistace = 1000.0;

	contactNum = 0; // Initialize number of contacts to 0
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

// ML-RFT coefficients
void world::ReadCoefficientData(elasticRod &m_rod)
{
	// Loading NN model
	if (translation == 1)
	{
		modelname = "/home/sciws001/Documents/mlhydro/MLRFT_distribution/Translation/after_server_datafull_translation_case1";
		// Set the model name as the full path to the trained translation machine learning model.
	}
	if (translation == 0)
	{
		modelname = "/home/sciws001/Documents/mlhydro/MLRFT_distribution/Rotation/after_server_datafull_datafull_case1";
		// Set the model name as the full path to the trained rotation machine learning model.
	}
	cppflow::model model(modelname);
	rod = &m_rod;
	// Preparing the input to the NN
	numVert = rod->nv - 2;

	s1_NN = linspace(0, deltaLengthInput * (numVert - 1), numVert);
	s2_NN = linspace(deltaLengthInput * (numVert - 1), 0, numVert);

	// Defining tangent and perp. components
	Vnorm_NN = VectorXd::Zero(numVert);
	VelocityT_NN = VectorXd::Zero(numVert);
	VelocityP_NN = VectorXd::Zero(numVert);
	p1_NN.setZero(3);
	p2_NN.setZero(3);
	p_NN.setZero(3);

	for (int i = 2; i < rod->nv; i++)
	{
		// std::cout << i << std::endl;
		u_NN(0) = 0.0;
		u_NN(1) = 0.0;
		u_NN(2) = 1.0;
		p_NN(0) = rod->getVertex(i)(0);
		p_NN(1) = rod->getVertex(i)(1);
		p_NN(2) = rod->getVertex(i)(2);

		if (translation == 1)
		{
			Vel_NN = u_NN;
		}
		if (translation == 0)
		{
			Vel_NN = u_NN.cross(p_NN);
		}
		Vnorm_NN(i - 2) = Vel_NN.norm();

		if (i > 2 && i < rod->nv - 1)
		{

			p1_NN = rod->getVertex(i - 1);
			p2_NN = rod->getVertex(i + 1);
			t_NN = p2_NN - p1_NN;
			t_NN = t_NN / t_NN.norm();
		}
		else if (i == 2)
		{

			p1_NN = rod->getVertex(i);
			p2_NN = rod->getVertex(i + 1);

			t_NN = p2_NN - p1_NN;
			t_NN = t_NN / t_NN.norm();
		}
		else
		{

			p1_NN = rod->getVertex(i - 1);
			p2_NN = rod->getVertex(i);
			t_NN = p2_NN - p1_NN;
			t_NN = t_NN / t_NN.norm();
		}

		VelocityT_NN(i - 2) = Vel_NN.dot(t_NN);

		vP_NN = Vel_NN - ((Vel_NN.dot(t_NN)) * t_NN);
		VelocityP_NN(i - 2) = vP_NN.norm();
	}

	// Obtain the coefficients

	inputNN = vector<float>(10, 0);
	for (int c = 0; c < numVert; ++c)
	{
		inputNN[0] = helixpitchRatio * rodRadiusRatio;
		inputNN[1] = 1 / contourRatio;
		inputNN[2] = VelocityT_NN(c) / Vnorm_NN(c);
		inputNN[3] = VelocityP_NN(c) / Vnorm_NN(c);
		inputNN[4] = s1_NN(c) / flagellaLength;
		inputNN[5] = s2_NN(c) / rodRadius;
		inputNN[6] = 1 / (contourRatio * rodRadiusRatio);
		inputNN[7] = 0;
		inputNN[8] = 0;
		inputNN[9] = 1; // Converted RPM to rad/s

		auto input_to_NN = cppflow::tensor({inputNN[0], inputNN[1], inputNN[2], inputNN[3], inputNN[4], inputNN[5], inputNN[6], inputNN[7], inputNN[8], inputNN[9]});
		input_to_NN = cppflow::cast(input_to_NN, TF_DOUBLE, TF_FLOAT);
		input_to_NN = cppflow::expand_dims(input_to_NN, 0);
		auto output_NN = model({{"serving_default_dense_input", input_to_NN}}, {"StatefulPartitionedCall:0"}); // Run forward pass through the NN

		output_coef = output_NN[0].get_data<float>();
		Ct.push_back(output_coef[0]);
		Cp.push_back(output_coef[1]);
		Cz.push_back(output_coef[2]);
	}

	Ct_e = Eigen::VectorXd::Zero(numVert);
	Cp_e = Eigen::VectorXd::Zero(numVert);
	Cz_e = Eigen::VectorXd::Zero(numVert);

	Ct_norm = Eigen::VectorXd::Zero(numVert);
	Cp_norm = Eigen::VectorXd::Zero(numVert);
	Cz_norm = Eigen::VectorXd::Zero(numVert);

	for (int i = 0; i < numVert; ++i)
	{
		Ct_e(i) = Ct[i];
		Cp_e(i) = Cp[i];
		Cz_e(i) = Cz[i];
		Ct_norm(i) = Ct[i];
		Cp_norm(i) = Cp[i];
		Cz_norm(i) = Cz[i];
	}
	if (translation == 1)
	{

		// Output_coef * ystd --- ystd = [0.005623743550026,0.00371828939389,0.197145491977754]
		Ct_e = Ct_e * 0.00222720644663503;
		Cp_e = Cp_e * 0.00263120201277278;
		Cz_e = Cz_e * 0.0798620080306157;

		// mean + Output_coef * ystd --- mean = [7.63807726394991,7.32731453922228,0.756079435979494]
		for (int i = 0; i < numVert; ++i)
		{
			Ct_e(i) = Ct_e(i) + 7.57259616888332;
			Cp_e(i) = Cp_e(i) + 7.41220021815048;
			Cz_e(i) = Cz_e(i) + 0.702476212430043;
		}

		// e^(output_coef*ystd + ymean)  + ymin - 2*Vector3d::Ones() --- ymin = [-2073.7,-1517.8,9.1077E-18]
		for (int i = 0; i < numVert; ++i)
		{
			Ct_e(i) = std::exp(Ct_e(i)) - 1944.1;
			Cp_e(i) = std::exp(Cp_e(i)) - 1655.0;
			Cz_e(i) = std::exp(Cz_e(i)) - 2.0;
		}
	}

	if (translation == 0)
	{

		// Output_coef * ystd --- ystd = [0.005623743550026,0.00371828939389,0.197145491977754]
		Ct_e = Ct_e * 0.005623743550026;
		Cp_e = Cp_e * 0.00371828939389;
		Cz_e = Cz_e * 0.197145491977754;

		// mean + Output_coef * ystd --- mean = [7.63807726394991,7.32731453922228,0.756079435979494]
		for (int i = 0; i < numVert; ++i)
		{
			Ct_e(i) = Ct_e(i) + 7.63807726394991;
			Cp_e(i) = Cp_e(i) + 7.32731453922228;
			Cz_e(i) = Cz_e(i) + 0.756079435979494;
		}

		// e^(output_coef*ystd + ymean)  + ymin - 2*Vector3d::Ones() --- ymin = [-2073.7,-1517.8,9.1077E-18]
		for (int i = 0; i < numVert; ++i)
		{
			Ct_e(i) = std::exp(Ct_e(i)) - 2075.7;
			Cp_e(i) = std::exp(Cp_e(i)) - 1519.8;
			Cz_e(i) = std::exp(Cz_e(i)) + 9.1077E-18 - 2;
		}
	}

	// Storing denormalized coef values back into std::vector
	for (int i = 0; i < numVert; ++i)
	{
		Ct[i] = -Ct_e(i);
		Cp[i] = Cp_e(i);
		Cz[i] = Cz_e(i);
		std::cout << Ct[i] << "," << Cp[i] << "," << Cz[i] << std::endl;
	}
}

void world::ReadOmegaData()
{
	ifstream infile;

	infile.open(inputName.c_str());
	if (!infile.is_open())
	{
		cout << "Unable to open file to read omega";
		timeStep = Nstep; // we are exiting
	}

	numOmegaPoints = 0;
	double a, b;
	while (infile >> a >> b)
	{
		numOmegaPoints++;
		timeSeries.push_back(a);
		omegaSeries.push_back(b);
	}
	infile.close();

	currentOmegaIndex = 0; // keeps track of current angular velocity
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData == false)
		return;

	int systemRet = system("mkdir datafiles"); // make the directory
	if (systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	ReadOmegaData();

	ostringstream name;
	name.precision(4);
	name << fixed;
	name << "datafiles/simDER";
	name << "_numvertex_" << rod->nv;
	// name << "_axisLength_" << axisLengthInput;
	name << "_rodradiusRatio_" << rodRadiusRatio;
	name << "_contourRatio_" << contourRatio;
	name << "_helixPitchRatio_" << helixpitchRatio;
	name << "_helixRadius_" << helixradius;
	name << "_omega_" << omegaSeries[currentOmegaIndex];
	name << "_useRFT" << useRFT;
	name << "_useRSS" << useRSS;
	name << "_totalTime_render_" << totalTime;
	name << "simDER.txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData == false)
		return;
	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData == false)
	{
		return;
	}

	//  data output every 0.01 seconds.

	if (fmod(timeStep, 10) == 0)
	{
		double sumFx = 0;
		double sumFy = 0;
		double sumFz = 0;

		for (int i = 2; i < rod->nv; i++)
		{
			sumFx += reactionForce(4 * i);
			sumFy += reactionForce(4 * i + 1);
			sumFz += reactionForce(4 * i + 2);
		}

		Vector3d rxn = Vector3d::Zero();
		Vector3d Torque = Vector3d::Zero();
		double sumtorque = 0.0;
		double omega = 0.0;

		if (translation == 1)
		{
			omega = omegaSeries[currentOmegaIndex];
		}
		else
		{
			omega = omegaSeries[currentOmegaIndex] * (2.0 * M_PI / 60.0);
		}

		for (int i = 2; i < rod->nv; i++)
		{
			Vector3d xCurrent = rod->getVertex(i);
			rxn(0) = reactionForce(4 * i);
			rxn(1) = reactionForce(4 * i + 1);
			rxn(2) = reactionForce(4 * i + 2);
			Torque += xCurrent.cross(rxn);
			sumtorque = Torque(2);
		}

		if (translation == 1)
		{
			outfile << currentTime << ", " << sumFx / (viscosity * omega * helixradius) << ", " << sumFy / (viscosity * omega * helixradius) << ", " << sumFz / (viscosity * omega * helixradius) << endl;
		}
		else
		{
			outfile << currentTime << ", 0, " << sumFx / (viscosity * omega * helixradius * helixradius) << ", " << sumFy / (viscosity * omega * helixradius * helixradius) << ", " << sumFz / (viscosity * omega * helixradius * helixradius) << ", " << sumtorque / (viscosity * omega * helixradius * helixradius * helixradius) << endl;
		}
	}
}

void world::setRodStepper()
{
	rodGeometry();

	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
						 youngM, shearM, RodLength);

	rodBoundaryCondition();

	rod->setup();

	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);

	if (includeContact == true) // If contact should be included, declare that force
	{
		m_externalContactForce = new externalContactForce(*rod, *stepper);
	}

	// dampingForce added
	if (useRSS == true && useRFT == false)
	{
		m_RegularizedStokeslet = new RegularizedStokeslet(*rod, *stepper, viscosity, epsilon);
	}
	else if (useRSS == false && useRFT == false)
	{
		ReadCoefficientData(*rod);
		m_dampingForce = new dampingForce(*rod, *stepper, viscosity, Ct, Cp, Cz);
	}
	else if (useRSS == false && useRFT == true)
	{
		m_resistiveForce = new resistiveForce(*rod, *stepper, viscosity, eta_per, eta_par);
	}
	else if (useRSS == true && useRFT == true)
	{
		std::cout << "In correct option, please check useRSS useRFT" << endl;
	}
	rod->updateTimeStep();

	timeStep = 0;
	currentTime = 0.0;

	Nstep = totalTime / deltaTime;

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius, 4) / 4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	int totVertices = 1 * (numVertices) + 2;
	ne = totVertices - 1;
	nv = totVertices;
}

// Setup geometry
void world::rodGeometry()
{
	double helixA = helixradius;
	double helixB = helixpitch / (2.0 * M_PI);

	int newNe = round(flagellaLength / deltaLengthInput);

	int numVertices = newNe;
	int totVertices = 1 * (numVertices) + 2; // NOTE THIS CHANGE
	vertices = MatrixXd(totVertices, 3);

	double helixT = flagellaLength / sqrt(helixA * helixA + helixB * helixB);
	double delta_t = helixT / (numVertices - 3); // step for t->[0, T]

	// geometry of helix
	double delta_l = flagellaLength / (numVertices - 1);
	double RodLength = flagellaLength * 1 + helixradius + helixA; // RODLENGTH CORRECTED

	vertices(0, 0) = 0.0;
	vertices(0, 1) = 0.0;
	vertices(0, 2) = helixA;
	;

	vertices(1, 0) = 0.0;
	vertices(1, 1) = 0.0;
	vertices(1, 2) = 0.0;

	int j = 0;
	for (double tt = 0.0; j < numVertices; tt += delta_t)
	{
		vertices(j + 2, 0) = helixA * cos(tt);
		vertices(j + 2, 1) = helixA * sin(tt);
		vertices(j + 2, 2) = -(helixB * tt);
		j++;
	}
}

void world::rodBoundaryCondition()
{

	rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	rod->setThetaBoundaryCondition(rod->getTheta(0), 0);
	rod->setVertexBoundaryCondition(rod->getVertex(1), 1);
}

void world::updateTimeStep()
{
	if (currentOmegaIndex < numOmegaPoints - 1 && timeSeries[currentOmegaIndex + 1] <= currentTime)
	{
		currentOmegaIndex++;
	}

	if (translation == false)
	{
		deltaTwist = omegaSeries[currentOmegaIndex] * (2.0 * M_PI / 60.0) * deltaTime;
		rod->setThetaBoundaryCondition(rod->getTheta(0) + deltaTwist, 0);
		rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	}
	if (translation == true)
	{
		nextpos = Vector3d::Zero(3);
		deltaTwist = omegaSeries[currentOmegaIndex] * deltaTime;
		nextpos[2] = nextpos[2] + deltaTwist;
		rod->setVertexBoundaryCondition(rod->getVertex(0) + nextpos, 0);
		rod->setThetaBoundaryCondition(rod->getTheta(0), 0);
		rod->setVertexBoundaryCondition(rod->getVertex(1) + nextpos, 1);
	}

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt

	if (includeContact == true)
		m_externalContactForce->setZeroForce();

	// compute hydrodynamic force - viscous force appllied. make change on here. Write your own. Just add viscous force.
	if (useRSS == true)
		m_RegularizedStokeslet->prepareForViscousForce();

	// solve the EOM
	updateEachRod();

	// compute contact
	if (includeContact == true)
	{
		prepareForContact();
		int contiter = 0;
		// resolve contact
		while (contactNum != 0 && contiter <= maxIterContact)
		{
			if (render == 1)
			{
				cout << "Contact detected, contact number = " << contactNum << endl;
			}

			rod->updateGuess();
			updateEachRod();

			contiter = contiter + 1;
		}
	}

	computeReactionForce();

	// update time step
	rod->updateTimeStep();

	currentTime += deltaTime;

	timeStep++;
}

void world::updateEachRod()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;

	bool solved = false;

	iter = 0;

	while (solved == false)
	{
		rod->prepareForIteration();

		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();

		m_stretchForce->computeFs();
		m_stretchForce->computeJs();

		m_bendingForce->computeFb();
		m_bendingForce->computeJb();

		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		if (includeContact == true)
			m_externalContactForce->computeFc();

		if (useRSS == true && useRFT == false)
			m_RegularizedStokeslet->computeFrs();
		else if (useRSS == false && useRFT == false)
		{
			m_dampingForce->computeFd(); // DampingForce added
			m_dampingForce->computeJd(); // DampingForce added
		}
		else if (useRSS == false && useRFT == true)
		{
			m_resistiveForce->computeFrft(); // DampingForce added
			m_resistiveForce->computeJrft(); // DampingForce added
		}
		// Compute norm of the force equations.
		normf = 0.0;
		for (int i = 0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}

		normf = sqrt(normf);

		if (iter == 0)
			normf0 = normf;

		if (normf <= forceTol)
		{
			solved = true;
		}
		else if (iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce);

			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	if (render)
	{
		cout << "Time: " << currentTime << " iter=" << iter << endl;
	}

	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep)
		return 1;
	else
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int j)
{
	return rod->x[j] / RodLength * 1.5;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

void world::prepareForContact()
{
	contactNum = 0;

	for (int j = 0; j < ne; j++)
	{
		const Vector3d x_1 = rod->getVertex(j);
		const Vector3d x_2 = rod->getVertex(j + 1);

		for (int l = 0; l < ne; l++)
		{
			// Two edges side by side are always going to contact. We will ignore it.
			if (abs(l - j) <= 1)
				continue;

			const Vector3d x_3 = rod->getVertex(l);
			const Vector3d x_4 = rod->getVertex(l + 1);

			// compute min length of two segements
			Vector3d c1;
			Vector3d c2;
			double s, t;

			double sqrdist = ClosestPtSegmentSegment(x_1, x_2, x_3, x_4, s, t, c1, c2);

			if (sqrdist < (2.0 * rodRadius) * (2.0 * rodRadius))
			{
				double pen = rodRadius + rodRadius - sqrt(sqrdist);

				Vector3d n = c2 - c1; // contact normal
				n.normalize();

				double wi = s;
				double wj = t;
				double d = 2.0 * rodRadius;
				double mdij = d - pen;

				double del_r_i = 0.5 * (mdij - d) * wi;
				double del_r_ip1 = 0.5 * (mdij - d) * (1.0 - wi);

				double del_r_j = 0.5 * (d - mdij) * wj;
				double del_r_jp1 = 0.5 * (d - mdij) * (1.0 - wj);

				double mi = rod->massArray(4 * j);
				double mip1 = rod->massArray(4 * (j + 1));

				double mj = rod->massArray(4 * l);
				double mjp1 = rod->massArray(4 * (l + 1));

				Vector3d f1 = n * del_r_i * mi / (deltaTime * deltaTime);
				Vector3d f2 = n * del_r_ip1 * mip1 / (deltaTime * deltaTime);

				Vector3d f3 = n * del_r_j * mj / (deltaTime * deltaTime);
				Vector3d f4 = n * del_r_jp1 * mjp1 / (deltaTime * deltaTime);

				m_externalContactForce->getContactForce(j + 0, f1);
				m_externalContactForce->getContactForce(j + 1, f2);
				m_externalContactForce->getContactForce(l + 0, f3);
				m_externalContactForce->getContactForce(l + 1, f4);

				contactNum = contactNum + 1;
			}
		}
	}
}

double world::ClosestPtSegmentSegment(const Vector3d &p1, const Vector3d &q1, const Vector3d &p2, const Vector3d &q2, double &s, double &t, Vector3d &c1, Vector3d &c2)
{
	Vector3d d1 = q1 - p1; // Direction vector of segment S1
	Vector3d d2 = q2 - p2; // Direction vector of segment S2
	Vector3d r = p1 - p2;
	double a = d1.dot(d1); // Squared length of segment S1, always nonnegative
	double e = d2.dot(d2); // Squared length of segment S2, always nonnegative
	double f = d2.dot(r);

	// Check if either or both segments degenerate into points
	if (a <= smallDist2 && e <= smallDist2)
	{
		// Both segments degenerate into points
		s = t = 0.0;
		c1 = p1;
		c2 = p2;
		return (c1 - c2).dot(c1 - c2);
	}

	if (a <= smallDist2)
	{
		// First segment degenerates into a point
		s = 0.0;
		t = f / e; // s = 0 => t = (b*s + f) / e = f / e
		double zero = 0.0, one = 1.0;
		t = clamp(t, zero, one);
	}
	else
	{
		double c = d1.dot(r);
		if (e <= smallDist2)
		{
			// Second segment degenerates into a point
			t = 0.0;
			s = clamp(-c / a, 0.0, 1.0); // t = 0 => s = (b*t - c) / a = -c / a
		}
		else
		{
			// The general nondegenerate case starts here
			double b = d1.dot(d2);
			double denom = a * e - b * b; // Always nonnegative

			// If segments not parallel, compute closest point on L1 to L2, and
			// clamp to segment S1. Else pick arbitrary s (here 0)
			if (denom != 0.0)
			{
				s = clamp((b * f - c * e) / denom, 0.0, 1.0);
			}
			else
				s = 0.0;

			// Compute point on L2 closest to S1(s) using
			// t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
			t = (b * s + f) / e;

			// If t in [0,1] done. Else clamp t, recompute s for the new value
			// of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
			// and clamp s to [0, 1]
			if (t < 0.0)
			{
				t = 0.0;
				s = clamp(-c / a, 0.0, 1.0);
			}
			else if (t > 1.0)
			{
				t = 1.0;
				s = clamp((b - c) / a, 0.0, 1.0);
			}
		}
	}

	c1 = p1 + d1 * s;
	c2 = p2 + d2 * t;
	return (c1 - c2).dot(c1 - c2);
}

void world::computeReactionForce()
{

	reactionForce = VectorXd::Zero(rod->ndof);

	if (useRSS == true && useRFT == false) // RSS
	{
		reactionForce = m_RegularizedStokeslet->ForceVec;
	}

	else if (useRSS == false && useRFT == true) // RFT
	{
		reactionForce = m_resistiveForce->ForceVec;
	}

	else if (useRSS == false && useRFT == false)
	{
		reactionForce = m_dampingForce->ForceVec;
	}

	else if (useRSS == true && useRFT == true)
	{
		std::cout << "In correct option, please check useRSS useRFT" << endl;
	}
}
