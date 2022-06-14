/**
 * simDER
 * simDER stands for "[sim]plified [D]iscrete [E]lastic [R]ods"
 * Dec 2017
 * This code is based on previous iterations.
 * */

// This line is for mac
//#include <GLUT/glut.h>

// This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "setInput.h"
#include "cppflow/cppflow.h"

world myWorld;
int NPTS;
ofstream outfile;

static void Key(unsigned char key, int x, int y)
{
	switch (key) // ESCAPE to quit
	{
	case 27:
		exit(0);
	}
}

/* Initialize OpenGL Graphics */
void initGL()
{
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
	glClearDepth(2.0f);					  // Set background depth to farthest
	// glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	// glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);						   // Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); // Nice perspective corrections

	glLoadIdentity();
	gluLookAt(0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 10);
	// void gluLookAt(	GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
	//		GLdouble centerX, GLdouble centerY, GLdouble centerZ, GLdouble upX,
	//		GLdouble upY, GLdouble upZ);

	glPushMatrix();

	// glMatrixMode(GL_MODELVIEW);
}

void display(void)
{
	while (myWorld.simulationRunning() > 0)
	{
		//  Clear screen and Z-buffer
		glClear(GL_COLOR_BUFFER_BIT);

		// draw axis
		double axisLen = 2.0;
		glLineWidth(0.5);

		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(axisLen, 0.0, 0.0);

		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, axisLen, 0.0);

		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, axisLen);
		glEnd();

		// draw a line
		glColor3f(0.1, 0.1, 0.1);
		glLineWidth(3.0);

		glBegin(GL_LINES);
		for (int j = 0; j < NPTS - 1; j++)
		{
			glVertex3f(myWorld.getScaledCoordinate(4 * j), myWorld.getScaledCoordinate(4 * j + 1), myWorld.getScaledCoordinate(4 * j + 2));
			glVertex3f(myWorld.getScaledCoordinate(4 * (j + 1)), myWorld.getScaledCoordinate(4 * (j + 1) + 1), myWorld.getScaledCoordinate(4 * (j + 1) + 2));
		}
		glEnd();

		glFlush();

		// Update step
		myWorld.updateTimeStep();
		myWorld.CoutData(outfile);
	}
	exit(1);
}

int main(int argc, char *argv[])
{

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	setInput inputData;
	inputData = setInput();
	inputData.LoadOptions(argv[1]);
	inputData.LoadOptions(argc, argv);
	// read input parameters from txt file and cmd

	myWorld = world(inputData);
	myWorld.setRodStepper();

	myWorld.OpenFile(outfile);

	bool render = myWorld.isRender();
	if (render) // if OpenGL visualization is on
	{
		NPTS = myWorld.numPoints();

		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500);
		glutInitWindowPosition(150, 100);
		glutCreateWindow("simDER");
		initGL();
		glutKeyboardFunc(Key);
		glutDisplayFunc(display);
		glutMainLoop();
	}
	else
	{
		while (myWorld.simulationRunning() > 0)
		{
			myWorld.updateTimeStep();  // update time step
			myWorld.CoutData(outfile); // write data to file
		}
	}

	// Close (if necessary) the data file
	myWorld.CloseFile(outfile);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
	return 0;
}
