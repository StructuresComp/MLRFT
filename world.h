#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "RegularizedStokeslet.h" // RSS-based viscous force
#include "dampingForce.h"		  // MLRFT-based viscous force
#include "ResistiveForce.h"		  // RFT-based viscous force
#include "externalContactForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

// include cppflow
#include "cppflow/cppflow.h"


// clamp function is copied from BASim code (Columbia Univ.)
/** Clamps scalar to the range [min,max]. */
template <typename T>
T clamp(T &scalar, T &min, T &max) // template <typename T> inline T clamp( const T& scalar, const T& min, const T& max)
{
	if (scalar < min)
		return min;
	if (scalar > max)
		return max;
	return scalar;
}

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int j);
	double getCurrentTime();
	double getTotalTime();
	std::vector<double> Ct, Cp, Cz;
	bool isRender();

	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

private:
	// Physical parameters
	double RodLength;
	double helixradius, helixpitchRatio, helixpitch;
	double rodRadiusRatio, rodRadius;
	double deltaLengthInput;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;
	double epsilon;
	double distf;
	double nTurn;
	double flagellaLength, contourRatio;
	string modelname;
	// Damping force parameters
	double eta_per;
	double eta_par;
	double headSize, omega;
	double C_translation, C_rotation; // numerical prefactor to multiply the force against head rotation and translation
	bool translation;				  // option for MLRFT-translation model load.
	double tol, stol;
	int maxIter; // maximum number of iterations
	int maxIterContact;
	double characteristicForce;
	double forceTol;

	// Geometry
	MatrixXd vertices;
	double currentTime;

	// Rod
	elasticRod *rod;

	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;

	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;

	bool useRSS, useRFT;		  // if true, RSS will be used. RFT will be used, if both false, MLRFT used.
	dampingForce *m_dampingForce; // Resistive force theory based model
	resistiveForce *m_resistiveForce;
	RegularizedStokeslet *m_RegularizedStokeslet; // Regularized Stokeslet Segment (RSS) based model

	bool includeContact;						  // if true, contact will be handled. Otherwise, it will be ignored.
	externalContactForce *m_externalContactForce; // contact force

	inertialForce *m_inertialForce;		  // inertial force
	externalGravityForce *m_gravityForce; // gravity

	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();

	// Variables about angular velocity
	double deltaTwist; // actuation
	string inputName;  // name of file containing omega
	string inputName2; // name of file containing coefficient based RSS
	std::vector<double> timeSeries, omegaSeries;
	int numOmegaPoints, currentOmegaIndex;
	void ReadOmegaData();
	void ReadCoefficientData(elasticRod &m_rod);
	bool render;   // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	void updateEachRod();

	double axisLength, axisLengthInput;

	double ClosestPtSegmentSegment(const Vector3d &p1, const Vector3d &q1, const Vector3d &p2, const Vector3d &q2, double &s, double &t, Vector3d &c1, Vector3d &c2);
	double smallDist2;

	void prepareForContact();
	int ne, nv;

	int contactNum;

	double minDistace;

	VectorXd reactionForce;

	void computeReactionForce();

	// input to NN
	std::vector<float> inputNN;
	int numVert;
	VectorXd s1_NN, s2_NN, Vnorm_NN;
	VectorXd VelocityT_NN, VelocityP_NN;
	VectorXd ForceT_RSS, ForceN_RSS, ForceZ_RSS;
	Vector3d vP_NN, t_NN, u_N, u_NN, p1_NN, p2_NN, p_NN, Vel_NN;

	// denormalizing the output coefs
	VectorXd Ct_e, Cp_e, Cz_e, Ct_norm, Cp_norm, Cz_norm;
	vector<float> output_coef;

	Vector3d xCurrent1;
	Vector3d xCurrent2;
	Vector3d rcurrent;
	Vector3d Torque;
	double Torquenorm, Forcehead, normTorque, Forcenorm;

	Vector3d f_0_rod_0;
	Vector3d f_1_rod_0;
	Vector3d f_2_rod_0;
	Vector3d f_rod_0;

	Vector3d f_0_rod_1;
	Vector3d f_1_rod_1;
	Vector3d f_2_rod_1;
	Vector3d f_rod_1;

	Vector3d nextpos;
};

#endif
