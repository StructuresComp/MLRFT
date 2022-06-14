#ifndef SETINPUT_H
#define SETINPUT_H

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "Option.h"
#include "eigenIncludes.h"

class setInput
{
public:
	typedef std::map<std::string, Option> OptionMap;
	OptionMap m_options;

	setInput();
	~setInput();

	template <typename T>
	int AddOption(const std::string &name, const std::string &desc, const T &def);

	Option *GetOption(const std::string &name);
	bool &GetBoolOpt(const std::string &name);
	int &GetIntOpt(const std::string &name);
	double &GetScalarOpt(const std::string &name);
	Vector3d &GetVecOpt(const std::string &name);
	string &GetStringOpt(const std::string &name);

	int LoadOptions(const char *filename);
	int LoadOptions(const std::string &filename)
	{
		return LoadOptions(filename.c_str());
	}
	int LoadOptions(int argc, char **argv);

private:
	double helixradius;
	double helixpitchRatio;
	double rodRadiusRatio;
	double contourRatio;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double tol, stol;

	double distance;

	// Damping force Input
	double headSize;

	int maxIter; // maximum number of iterations
	int maxIterContact;
	double density;
	double epsilon;
	Vector3d gVector;
	double viscosity;
	bool render;
	bool saveData;
	std::string inputName;			  // name of file containing omega
	std::string inputName2;			  // name of file contaiining the coefficient for ML RSS
	bool useRSS, useRFT;			  // should RSS/RFT be used?
	bool translation;				  // if 1 translation MLRFT model being used. if 0 rotation.
	bool includeContact;			  // should contact be handled?
	double C_translation, C_rotation; // numerical prefactor for translation and rotation
};

#include "setInput.tcc"

#endif // PROBLEMBASE_H
