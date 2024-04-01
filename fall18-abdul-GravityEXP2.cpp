// This file is part of CNCSVision, a computer vision related library
// This software is developed under the grant of Italian Institute of Technology
//
// Copyright (C) 2011 Carlo Nicolini <carlo.nicolini@iit.it>
//
// CNCSVision is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// CNCSVision is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// CNCSVision. If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <queue>

/********* BOOST MULTITHREADED LIBRARY ****************/
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>	//include asio in order to avoid the "winsock already declared problem"

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#endif

#ifdef __linux__
#include <GL/glut.h>
#include <SOIL/SOIL.h>
#endif

#ifdef _WIN32
#include <windows.h>
#include <gl\gl.h>            // Header File For The OpenGL32 Library
#include <gl\glu.h>            // Header File For The GLu32 Library
#include "glut.h"            // Header File For The GLu32 Library
#include <MMSystem.h>
#include "SOIL.h"
#endif
/********* INCLUDE CNCSVISION LIBRARY HEADERS **********/
//#include "Optotrak.h"
#include "Optotrak2.h"
#include "Marker.h"
#include "Mathcommon.h"
#include "GLUtils.h"
#include "VRCamera.h"
#include "CoordinatesExtractor.h"
#include "CylinderPointsStimulus.h"
#include "EllipsoidPointsStimulus.h"
#include "StimulusDrawer.h"
#include "GLText.h"
#include "BalanceFactor.h"
#include "ParStaircase.h"
#include "Staircase.h"
#include "ParametersLoader.h"
#include "TrialGenerator.h"
#include "Util.h"
#include "BrownMotorFunctions.h"
#include "BrownPhidgets.h"

/***** CALIBRATION FILE *****/
#include "LatestCalibration.h"

/***** DEFINE SIMULATION *****/
//#define SIMULATION
#ifndef SIMULATION
#include <direct.h> // mkdir
#endif

/********* NAMESPACE DIRECTIVES ************************/
using namespace std;
using namespace mathcommon;
using namespace Eigen;
using namespace util;
using namespace BrownMotorFunctions;
using namespace BrownPhidgets;

/********* #DEFINE DIRECTIVES **************************/
#define TIMER_MS 11                               // 85 hz
#define SCREEN_WIDTH  1024                  // 1024 pixels
#define SCREEN_HEIGHT 768                   // 768 pixels

double interoculardistance = 0;
const float DEG2RAD = 3.14159/180;
Screen screen;
double mirrorAlignment = 0.0;
double screenAlignmentY = 0.0;
double screenAlignmentZ = 0.0;

/********* VARIABLES OBJECTS  **************************/
VRCamera cam;
Optotrak2 optotrak;
CoordinatesExtractor headEyeCoords, thumbCoords, indexCoords, thumbJointCoords, indexJointCoords;
Timer timer;
Timer globalTimer;
clock_t t;

/********* VISUALIZATION AND STIMULI *******************/
static const bool gameMode=true;
static const bool stereo=true;

/********* MARKERS AND 3D VECTORS ****************************/
// fingers markers numbers
int ind0 = 3;
int ind1 = 13, ind2 = 14, ind3 = 16;
int thu1 = 15, thu2 = 17, thu3 = 18;
int calibration1 = 1, calibration2 = 2;
int screen1 = 19, screen2 = 20, screen3 = 21;
int mirror1 = 6, mirror2 = 22;
int centercalMarker = 4;

Vector3d eyeLeft, eyeRight;
Vector3d ind, thm;
Vector3d indexCalibrationPoint(0,0,0), thumbCalibrationPoint(0,0,0);

vector <Marker> markers;

/********* VISIBILITY BOOLEANS ************************/
bool markers_status = true;
bool allVisibleIndex=markers_status;
bool allVisibleThumb=markers_status;
bool allVisibleFingers=markers_status;

bool visibleInfo=true;

/********* FILE STREAMS *************************************/
ofstream trialFile;

/*************************************************************************************/
/*** Everything above this point stays more or less the same between experiments.  ***/
/*************************************************************************************/

/********* VARIABLES THAT CHANGE IN EACH EXPERIMENT *************************/
// Experiment variables
ParametersLoader parameters; //high level variables from parameters file
TrialGenerator<double> trial; //chnaged from balancefactor 

// FOR NEW CALIBRATION:
double markerXOffset = 10;
double markerYOffset = 10;

// Variables for counting trials, frames, and lost frames
int trialNumber = 0;
int frameN = 0;
int trialsPerBlock = 36;

// Flags for important states in the experiment
bool fingersCalibrated = true;
bool ProbePhase = false;
bool ProbeBallEdge = false;
bool CueBallEdge = false;


bool cueVelSet = false;
bool finished = false;
float ballPos_y;
float ballPos_z;

// Display
double displayDepth = -400;
//double distanceCue2TargetZ= 100;
// finger pos
float indXNow,indYNow,indZNow;

// Virtual target objects

//Ball
float cueCenter_x;
float cueCenter_y;
float cueCenter_z;
float cueRadius = 8;
float cueVel_x = 0;
float cueVel_y = 0;
float cueVel_z = 0;
float ballStartPos_z = displayDepth - 375;

//Probe ball 
float probeCenter_x;
float probeCenter_y;
float probeCenter_z;
float probeDistance;

//Table
float Tablex1 = -100;
float Tablex2 = 100; 
float Tabley1 = -70;
float TableZ1 = displayDepth-200;
float TableZ2 = displayDepth-400;

//leg
float Legx1 = -100;
float Legx2 = 100;
float Legy1 = Tabley1;
float Legy2 = Tabley1 - 67.83; //XX
float LegZ = TableZ1;

//floor
float Floorx1 = Tablex1;
float Floorx2 = Tablex2;
float Floory1 = Legy2;
float Floorz1 = LegZ + 170;
float Floorz2 = LegZ;

//Physics & misc
double speed;
double Velocity_y;
double Gravity;
int Phase;
int Order;
double probeSpeed; //subject's response changes this
bool response;
double timeStart;
double dist;
double lastFrameCue;
double lastTimeProbe;
double lastFrameProbe;
double ProbeFrame; //for gravity timing in probe 

///////////////////////////////////////


float distanceBetween; 
double elapsed;
double frameOfFall;
double lastFrame;
double timeOfImpact;
double Probe2CueDelay = 1000; 
float distanceToStart;




// Incremented when stepping thru calibration procedure
// Make sure that drawInfo() and handleKeypress() are in agreement about this variable!
int fingerCalibrationDone = 0;

GLfloat LightAmbient[] = {0.5f, 0.5f, 0.5f, 1.0f}; 
GLfloat LightDiffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
GLfloat LightSpecular[] = {1.0f, 1.0f, 1.0f, 1.0f};
GLfloat LightPosition[] = {-50.0f, 150.0f, -350.0f, 1.0f};
GLfloat specularMaterial[] = {1.0, 0.0, 0.0, 1.0};
GLfloat tableMaterial[] = {1.0, 0.0, 0.0, 1.0};
GLfloat ballMaterial[] = {1.0, 0.0, 0.0, 1.0};
GLfloat noiseMaterial[] = {1.0, 0.0, 0.0, 1.0};
GLfloat lineMaterial[] = {0.0, 0.0, 0.0, 1.0};
GLfloat shininessMaterial = 32.0f;
// To use the motors:
// 0. homeEverything(arm speed, screen speed) brings all motors to their start positions, nearest to the mirror.
// 1. Before trying to move, you must pick a point that is rigidly attached to the workspace robotic arm.
// 2. Call this location 'centercal' and simply hold onto this initial location of your 'point of interest'.
// 3. When moving the arm you will decide where to put this point in particular.
// 4. Provide the desired final location and the initial location of the point of interest: moveObjectAbsolute(moveTo, centercal)
// Repeat: Set centercal to be the starting position (after homeEverything) of the point you want to control, then select locations for this point
Vector3d centercal(0.0, 0.0, 0.0);
//double centercalToRotationAxisAtStart = 23; // centercal marker has a y-offset from the desired point of interest
// Settable Target Locations
Vector3d moveTo(0.0,-5000,-500);

/********** FUNCTION PROTOTYPES *****/
void advanceTrial();
void beepOk(int tone);
void calibration_fingers(int phase);
void cleanup();
void drawGLScene();
void drawInfo();
void drawFingers();
void drawStimulus();
void drawResponseGrid();
void handleKeypress(unsigned char k, int x, int y);
void handleResize(int w, int h);
void idle();
void initMotors();
void initOptotrak();
void initProjectionScreen(double _focalDist, const Affine3d &_transformation=Affine3d::Identity(),bool synchronous=true);
void initRendering();
void initStreams();
void initTrial();
void initVariables();
void update(int value);
void updateTheMarkers();

// online operations
bool cue();
bool probe();
bool sleep();
void online_apparatus_alignment();
void online_fingers();
void online_trial();


/*************************** EXPERIMENT SPECS ****************************/
// experiment directory
string experiment_directory = "R:/CLPS_Domini_Lab/abdul/fall18-GravityEXP2/";
// parameters file directory and name
string parametersFile_directory = experiment_directory + "fall18-GravityEXP2Parameters.txt";


/*************************** FUNCTIONS ***********************************/
// First, make sure the filenames in here are correct and that the folders exist.
// If you mess this up, data may not be recorded!
void initStreams()
{
	ifstream parametersFile;
	parametersFile.open(parametersFile_directory.c_str());
	parameters.loadParameterFile(parametersFile);

	string subjectName = parameters.find("SubjectName");
	string Phase_index = parameters.find("Phase");
	Phase = str2num<int>(Phase_index);

	// trialFile directory
	string dirName  = experiment_directory + subjectName;
	mkdir(dirName.c_str()); // windows syntax

	if (util::fileExists(dirName+"/"+subjectName + ".txt"))
	{
		string error_on_file_io = dirName+"/"+subjectName+".txt" + string(" already exists");
		cerr << error_on_file_io << endl;
		MessageBox(NULL, (LPCSTR)"FILE ALREADY EXISTS\n Please check the parameters file.",NULL, NULL);
		exit(0);
	}

	globalTimer.start();
	string trialFile_headers;
	//Will fix this later. Needs to change headers depending on testing phase. 
	// trial file headers
	if (Phase == 1){
		trialFile_headers = "subjName\ttrialN\tPhase\tspeed\telapsed\tframeN\tProbePhase\tProbeBallEdge\tprobeSpeed\tresponse\tballPos_z\tballPos_y";
	}
	else if (Phase == 2){
		trialFile_headers = "subjName\ttrialN\tPhase\tGravity\telapsed\tframeN\tProbePhase\tProbeBallEdge\tprobeSpeed\tresponse\tballPos_z\tballPos_y";
	}
	else{
		trialFile_headers = "subjName\ttrialN\tPhase\tspeed\tGravity\telapsed\tframeN\tProbePhase\tProbeBallEdge\tprobeSpeed\tresponse\tballPos_z\tballPos_y";
	}
	string trialFileName = dirName + "/" + subjectName + ".txt";
	trialFile.open(trialFileName.c_str());
	trialFile << fixed << trialFile_headers << endl;
}

// Edit case 'f' to establish calibration procedure
// Also contains other helpful button presses (ESC for quit, i for info)
void handleKeypress(unsigned char key, int x, int y){
	switch (key){

		case 'o':
		case 'O':
			{
				visibleInfo=!visibleInfo;
			}
			break;

		case 'm':
		case 'M':
			{
				interoculardistance += 0.5;
				headEyeCoords.setInterOcularDistance(interoculardistance);
			}
			break;

		case 'n':
		case 'N':
			{
				interoculardistance -= 0.5;
				headEyeCoords.setInterOcularDistance(interoculardistance);
			}
			break;

		case 27:	// ESC
			{   
				if(trialFile.is_open()){
					trialFile.close();
				}
				homeEverything(5000,4500);
				cleanup();
				exit(0);
			}
			break;

		case 'f':
		case 'F':
			{
				//if ( fingerCalibrationDone==0 && allVisibleFingers  )
				//{
					/* calibration on the X
					indexCalibrationPoint=markers.at(1).p;
					indexCalibrationPoint[0] = indexCalibrationPoint[0] - 7;
					indexCalibrationPoint[1] = indexCalibrationPoint[0] + 10;
					indexCoords.init(indexCalibrationPoint, markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );

					fingerCalibrationDone=3;
					fingersCalibrated=true;*/
					visibleInfo=false;
					beepOk(0);
					trialNumber++;
					initTrial();
					//break;
				//}
			}
			break;


		case '2': //lower speed 
			{
				if (elapsed > lastTimeProbe + 80){
					if (Order == 1){
					response = true;
					}
					else if (Order == 2){
					response = false;
					}
					beepOk(19);	
					advanceTrial(); 


				}
			}
			break;

		case '1': // raise speed
			{
				if (elapsed > lastTimeProbe + 80){
					if (Order == 1){
					response = false;
					}
					else if (Order == 2){
					response = true;
					}
					beepOk(19);	
					advanceTrial(); 
				}

			}
			break;
	}
}

/*** GRASP ***/
void calibration_fingers(int phase)
{

}


// Provide text instructions for calibration, as well as information about status of experiment
void drawInfo()
{
	if (finished)
		visibleInfo = true;

	if ( visibleInfo )
	{
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_BLEND);
		glDisable(GL_LIGHTING);
		GLText text;
		if ( gameMode )
			text.init(SCREEN_WIDTH,SCREEN_HEIGHT,glWhite,GLUT_BITMAP_HELVETICA_12);
		else
			text.init(640,480,glWhite,GLUT_BITMAP_HELVETICA_12);
		text.enterTextInputMode();

		if (finished) {
			glColor3fv(glWhite);
			text.draw("The experiment is over. Thank you! :)");
		}else{
			if(!fingersCalibrated){
				switch (fingerCalibrationDone)
				{
				case 0:
					text.draw("Press F when index finger is on the X.");
					break;
				} // end switch(fingerCalibrationDone)
			}

			/////// Header ////////
			text.draw("####### ####### #######");
			text.draw("#");
			text.draw("# Name: " + parameters.find("SubjectName"));
			text.draw("# IOD: " + stringify<double>(interoculardistance));
			text.draw("# trial: " + stringify<float>(trialNumber));
			text.draw("# Phase:" +stringify<int>(Phase));
			text.draw("# Order:" +stringify<int>(Order));
			text.draw("# Displayed Velocity:" + stringify<double>(speed));
			text.draw("# Displayed Acceleration:" + stringify<double>(Gravity));
			text.draw("# probeDistance:" + stringify<double>(probeDistance));
			text.draw("# ProbePhase? " + stringify<bool>(ProbePhase));
			text.draw("# CueBallEdge? " + stringify<bool>(CueBallEdge));
			text.draw("# Probe ball Edge? " + stringify<double>(ProbeBallEdge));
			text.draw("# ballPos_z " + stringify<float>(cueCenter_z));
			text.draw("# ballPos_y " + stringify<float>(cueCenter_y));
			text.draw("# ballPos_x " + stringify<float>(cueCenter_x));
			text.draw("# Response Velocity/Acceleration :" +stringify<float>(probeSpeed));
			text.draw("# response:" +stringify<bool>(response));



			/////// Mirror and Screen Alignment ////////
			if ( abs(mirrorAlignment - 45.0) < 0.2 )
				glColor3fv(glGreen);
			else
				glColor3fv(glRed);
			text.draw("# Mirror Alignment = " +stringify<double>(mirrorAlignment));

			if ( isVisible(markers[mirror1].p) )
				glColor3fv(glGreen);
			else
				glColor3fv(glRed);
			text.draw("Mirror 1 " +stringify< Eigen::Matrix<double,1,3> > (markers[mirror1].p.transpose()));

			if ( isVisible(markers[mirror2].p) )
				glColor3fv(glGreen);
			else
				glColor3fv(glRed);
			text.draw("Mirror 2 " +stringify< Eigen::Matrix<double,1,3> > (markers[mirror2].p.transpose()));

			if ( screenAlignmentY < 89.0 )
				glColor3fv(glRed);
			else
				glColor3fv(glGreen);
			text.draw("# Screen Alignment Y = " +stringify<double>(screenAlignmentY));
			if ( abs(screenAlignmentZ) < 89.0 )
				glColor3fv(glRed);
			else
				glColor3fv(glGreen);
			text.draw("# Screen Alignment Z = " +stringify<double>(screenAlignmentZ));


			glColor3fv(glWhite);
			if (fingerCalibrationDone==0){

				if ( isVisible(markers[ind0].p) )
					glColor3fv(glGreen);
				else
					glColor3fv(glRed);
				text.draw("Index Calibration Point " +stringify< Eigen::Matrix<double,1,3> > (markers[ind0].p.transpose()));

				///// INDEX FINGER ///////
				glColor3fv(glWhite);
				text.draw(" " );
				text.draw("Index" );
				if ( isVisible(markers[13].p) && isVisible(markers[14].p) && isVisible(markers[16].p) )
					glColor3fv(glGreen);
				else
					glColor3fv(glRed);
				text.draw("Marker " + stringify<int>(13)+stringify< Eigen::Matrix<double,1,3> > (markers[13].p.transpose())+ " [mm]" );
				text.draw("Marker " + stringify<int>(14)+stringify< Eigen::Matrix<double,1,3> > (markers[14].p.transpose())+ " [mm]" );
				text.draw("Marker " + stringify<int>(16)+stringify< Eigen::Matrix<double,1,3> > (markers[16].p.transpose())+ " [mm]" );

				/////// THUMB //////
				glColor3fv(glWhite);
				text.draw(" " );
				text.draw("Thumb" );
				if ( isVisible(markers[15].p) && isVisible(markers[17].p) && isVisible(markers[18].p) )
					glColor3fv(glGreen);
				else
					glColor3fv(glRed);
				text.draw("Marker " + stringify<int>(15)+stringify< Eigen::Matrix<double,1,3> > (markers[15].p.transpose())+ " [mm]" );
				text.draw("Marker " + stringify<int>(17)+stringify< Eigen::Matrix<double,1,3> > (markers[17].p.transpose())+ " [mm]" );
				text.draw("Marker " + stringify<int>(18)+stringify< Eigen::Matrix<double,1,3> > (markers[18].p.transpose())+ " [mm]" );
			}

			/////// Index and Thumb Positions ////////
			if (fingersCalibrated){
				glColor3fv(glWhite);
				text.draw("--------------------");
				if (allVisibleIndex)
					glColor3fv(glGreen);
				else
					glColor3fv(glRed);
				text.draw("Index= " +stringify< Eigen::Matrix<double,1,3> >(ind.transpose()));
				/*if (allVisibleThumb)
				glColor3fv(glGreen);
				else
				glColor3fv(glRed);
				text.draw("Thumb= " +stringify< Eigen::Matrix<double,1,3> >(thm.transpose()));*/
				glColor3fv(glWhite);
				text.draw("--------------------");
			}

			//////// OTHER INFO /////
			glColor3fv(glGreen);
			text.draw("Timer= " + stringify<int>(timer.getElapsedTimeInMilliSec()) );
			text.draw("Frame= " + stringify<int>(frameN));
			glColor3fv(glWhite);
			text.draw("--------------------");
		}
		text.leaveTextInputMode();
		glEnable(GL_LIGHTING);
		glEnable(GL_BLEND);
	}
}

// This will be called at 85hz in the main loop
// Not too much to change here usually, sub-functions do the work.
void drawGLScene() 
{
	online_apparatus_alignment();
	online_fingers();
	online_trial();

	if (stereo)
    {   glDrawBuffer(GL_BACK);
		// Draw left eye view
        glDrawBuffer(GL_BACK_LEFT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0,0.0,0.0,1.0);
        cam.setEye(eyeLeft);
        drawStimulus();
		drawInfo();

        // Draw right eye view
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0,0.0,0.0,0.0);
        cam.setEye(eyeRight);
        drawStimulus();
		drawInfo();

        glutSwapBuffers();
    }
    /*else
    {   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0,0.0,0.0,1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        cam.setEye(eyeRight);
        drawStimulus();
		drawInfo();
        glutSwapBuffers();
    }*/
}

// Can check for various conditions that might affect how the graphics in here
void drawStimulus()
{

	//1. Draw Table Surface 
	glBegin(GL_QUADS);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tableMaterial);
	glVertex3f(Tablex1, Tabley1,TableZ1); // vertex 1
	glVertex3f(Tablex2, Tabley1,TableZ1); // vertex 2
	glVertex3f(Tablex2, Tabley1,TableZ2); // vertex 3
	glVertex3f(Tablex1, Tabley1,TableZ2); // vertex 4
	glEnd();

	//2. Draw Table Leg	
	glBegin(GL_QUADS);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tableMaterial);
	glVertex3f(Tablex1, Legy1,LegZ ); // vertex 1
	glVertex3f(Tablex2, Legy1,LegZ); // vertex 2
	glVertex3f(Tablex2, Legy2,LegZ); // vertex 3
	glVertex3f(Tablex1, Legy2,LegZ); // vertex 4
	glEnd();

	//3. Draw Floor
	glBegin(GL_QUADS);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tableMaterial);
	glVertex3f(Floorx1, Floory1,Floorz1); // vertex 1
	glVertex3f(Floorx2, Floory1,Floorz1); // vertex 2
	glVertex3f(Floorx2, Floory1,Floorz2); // vertex 3
	glVertex3f(Floorx1, Floory1,Floorz2); // vertex 4
	glEnd();
	
	// 4. Draw cue ball
	if((Order == 1 && !ProbePhase && !cue())||(Order == 2 && !ProbePhase && !cue() &&(elapsed > lastTimeProbe + Probe2CueDelay))) {
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPushMatrix();
		glLoadIdentity();
		GLUquadricObj* cueBall = gluNewQuadric();
		gluQuadricDrawStyle(cueBall, GLU_FILL);
		glTranslated(cueCenter_x,cueCenter_y,cueCenter_z);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ballMaterial);
		gluSphere(cueBall, cueRadius, 64, 64);
		gluDeleteQuadric(cueBall);
		glPopMatrix();
	}
			
	// 5. Draw response ball
	if((Order == 2 && ProbePhase && !ProbeBallEdge)|| (Order == 1 && ProbePhase && !ProbeBallEdge && (elapsed > lastFrameCue + Probe2CueDelay))) {//  after display period NEEDS TO BE FIXED FOR OTHER ORDERES || (ProbePhase && Order == 2 )
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPushMatrix();
		glLoadIdentity();
		GLUquadricObj* cueBall = gluNewQuadric();
		gluQuadricDrawStyle(cueBall, GLU_FILL);
		glTranslated(probeCenter_x,probeCenter_y,probeCenter_z);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ballMaterial);
		gluSphere(cueBall, cueRadius, 64, 64);
		gluDeleteQuadric(cueBall);
		glPopMatrix();

		
	}
}

//Probe2CueDelay



// called at the beginning of every trial
void initTrial()
{
	// initializing all variables

	/*string subjectName = parameters.find("SubjectName");
	string trialFileName = experiment_directory_desktop + subjectName + "/" + subjectName + "_trial" + stringify<double>(trialNumber) + ".txt";
	trialFile.open(trialFileName.c_str());
	trialFile << fixed << trialFile_headers << endl;*/

	frameN=0;
	Order = trial.getCurrent().first["Order"];
	cueVelSet = false;
	ProbePhase = (Order == 1) ? false : true;
	CueBallEdge = false; 
	ProbeBallEdge = false;
	response = -1;

	//1. Horizontal Test for Vz
	if (Phase ==1){
		double speed_index = trial.getCurrent().first["Speed"];
		speed = speed_index;
		cueCenter_x = 0;
		cueCenter_y = -62;
		cueCenter_z = ballStartPos_z;

	} 
	//2.Horizontal Test for Vy
	else if (Phase == 2){
		Gravity = trial.getCurrent().first["Gravity"];
		cueCenter_x = 0;
		cueCenter_y = -62;
		cueCenter_z = TableZ1 +12;

	}

	//3. Horizontal Test, full trajectory 
	else {
		speed = trial.getCurrent().first["Speed"];
		Gravity = trial.getCurrent().first["Gravity"];
		cueCenter_x = 0;
		cueCenter_y = -62;
		cueCenter_z = ballStartPos_z;

		}
	probeDistance = rand() % 175 + 68 ;
	probeCenter_x = -1*(probeDistance/2);
	probeCenter_y = -62;
	probeCenter_z = TableZ1 -20;
	probeSpeed = trial.getCurrent().second->getCurrentStaircase()->getState();

	initProjectionScreen(displayDepth);

	// roll on
	drawGLScene();
	timer.start();
}

// This function handles the transition from the end of one trial to the beginning of the next.
void advanceTrial() {
	if (Phase == 1){
		trialFile << fixed <<
			parameters.find("SubjectName") << "\t" <<		//subjName
			trialNumber << "\t" <<							//trialN
			Phase << "\t" <<
			speed << "\t" <<
			elapsed << "\t" <<
			frameN << "\t" <<
			ProbePhase << "\t" <<
			ProbeBallEdge << "\t" <<
			probeSpeed << "\t" <<
			response << "\t"<<
			ballPos_z << "\t" <<
			ballPos_y << endl;
	}
	else if (Phase == 2){
		trialFile << fixed <<
			parameters.find("SubjectName") << "\t" <<		//subjName
			trialNumber << "\t" <<							//trialN
			Phase << "\t" <<
			Gravity << "\t" <<
			elapsed << "\t" <<
			frameN << "\t" <<
			ProbePhase << "\t" <<
			ProbeBallEdge << "\t" <<
			probeSpeed << "\t" <<
			response << "\t"<<
			ballPos_z << "\t" <<
			ballPos_y << endl;
	}
	else{
		trialFile << fixed <<
			parameters.find("SubjectName") << "\t" <<		//subjName
			trialNumber << "\t" <<							//trialN
			Phase << "\t" <<
			speed << "\t" <<
			Gravity << "\t" <<
			elapsed << "\t" <<
			frameN << "\t" <<
			ProbePhase << "\t" <<
			ProbeBallEdge << "\t" <<
			probeSpeed << "\t" <<
			response << "\t"<<
			ballPos_z << "\t" <<
			ballPos_y << endl;
	}

	if(!trial.isEmpty()){
		trial.next(response);
		trialNumber++;
		initTrial();
	}
	else{
		finished=true;
	}
	//if(trialFile.is_open())
	//	trialFile.close();
}

void idle() {

	elapsed = timer.getElapsedTimeInMilliSec();

	// get new marker positions from optotrak
	updateTheMarkers();

	// eye coordinates
	eyeRight = Vector3d(interoculardistance/2,0,0);//0
	eyeLeft = Vector3d(-interoculardistance/2,0,0);//0
}



/*** Online operations ***/
void online_apparatus_alignment()
{
	if(visibleInfo){
		// mirror alignment check
		mirrorAlignment = asin(
			abs((markers.at(mirror1).p.z()-markers.at(mirror2).p.z()))/
			sqrt(
			pow(markers.at(mirror1).p.x()-markers.at(mirror2).p.x(), 2) +
			pow(markers.at(mirror1).p.z()-markers.at(mirror2).p.z(), 2)
			)
			)*180/M_PI;

		// screen Y alignment check
		screenAlignmentY = asin(
			abs((markers.at(screen1).p.y()-markers.at(screen3).p.y()))/
			sqrt(
			pow(markers.at(screen1).p.x()-markers.at(screen3).p.x(), 2) +
			pow(markers.at(screen1).p.y()-markers.at(screen3).p.y(), 2)
			)
			)*180/M_PI;

		// screen Z alignment check
		screenAlignmentZ = asin(
			abs(markers.at(screen1).p.z()-markers.at(screen2).p.z())/
			sqrt(
			pow(markers.at(screen1).p.x()-markers.at(screen2).p.x(), 2) +
			pow(markers.at(screen1).p.z()-markers.at(screen2).p.z(), 2)
			)
			)*180/M_PI*
			abs(markers.at(screen1).p.x()-markers.at(screen2).p.x())/
			(markers.at(screen1).p.x()-markers.at(screen2).p.x());
	}
}

void online_fingers()
{
	// Visibility check
	allVisibleIndex = isVisible(markers.at(ind1).p) && isVisible(markers.at(ind2).p) && isVisible(markers.at(ind3).p);
	//allVisibleThumb = isVisible(markers.at(thu1).p) && isVisible(markers.at(thu2).p) && isVisible(markers.at(thu3).p);
	allVisibleFingers = allVisibleIndex;

	// fingers coordinates, fingersOccluded and framesOccluded
	if ( allVisibleFingers )
	{
		indexCoords.update(markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );
		//thumbCoords.update(markers.at(thu1).p, markers.at(thu2).p, markers.at(thu3).p );
		//indexJointCoords.update(markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );
		//thumbJointCoords.update(markers.at(thu1).p, markers.at(thu2).p, markers.at(thu3).p );
	}

	if(fingersCalibrated)
	{
		// index coordinates
		if(allVisibleIndex){
			ind = indexCoords.getP1();
			//if(fingerCalibrationDone>=3)
			//	indJoint = indexJointCoords.getP1();
		}
		// thumb coordinates
		//if(allVisibleThumb){
		//thm = thumbCoords.getP1();
		//if(fingerCalibrationDone>=3)
		//	thmJoint = thumbJointCoords.getP1();
		//}
	}
}

bool cue(){ // Updates the ball position and returns true only if the cue phase is done
	
	if(!ProbePhase){
		
		//Check for distance with table edge Cue
		if (Phase == 1){

			

			float distanceBetween_z = TableZ1 - cueCenter_z;
			if (distanceBetween_z <= 0){ 

				// when the ball gets to table edge
				lastFrameCue = elapsed;
				ballPos_z = cueCenter_z;
				ballPos_y = cueCenter_y;

			}else{
				//update ball positions		
				cueCenter_z += speed;
			}

			return distanceBetween_z <= 0; // Cue phase is over
		
		} else if (Phase == 2){
			//Check for distance with floor Cue

			float distanceBetween_y = Floory1 - cueCenter_y;
			
			if (cueCenter_y <= -129.83){ 
				// when the ball touches the ground 
				lastFrameCue = elapsed;
				cueCenter_y = -129.83;
				
			}else {
				// update cueball positions 
				if(Order == 1){
				cueCenter_y = (Tabley1+cueRadius) - 0.5*(Gravity/1000)*pow(11.76*(frameN -1), 2);
				}else if(Order ==2){
				cueCenter_y = (Tabley1+cueRadius) - 0.5*(Gravity/1000)*pow(11.76*(frameN -lastFrameProbe), 2);
				}
			}
			return distanceBetween_y >= -1*cueRadius;
		
		} else { // Phase 3
			
			//Check for distance with edge and Cue
			if(!CueBallEdge){
				float distanceBetween_z = TableZ1 - cueCenter_z;
				if ((distanceBetween_z <=(-.3*cueRadius))){ 
					// when the ball touches the ground 
					CueBallEdge = true;
					frameOfFall = frameN+1;
					std::cout << frameOfFall << "  " << frameN << std::endl;
				}

				// Update position for constant velocity
				cueCenter_z += speed;

			} else { // falling

				if(cueCenter_y <= -129.83){ // stopped falling

					ballPos_z = cueCenter_z;
					ballPos_y = cueCenter_y;
					lastFrame = frameN;
					lastFrameCue = elapsed; 

					return true; // Start probe phase

				} else { // still falling

					cueCenter_y = (Tabley1+cueRadius) - 0.5*(Gravity/1000)*pow(11.76*(frameN-frameOfFall), 2);
						if (cueCenter_y < Floory1 + cueRadius){
						cueCenter_y = Floory1 + cueRadius;
						cueCenter_z = (speed/11.76)* sqrt(75*(2000/Gravity)) -600;
					}

					return cueCenter_y <= -129.83;
				}		
			}
		}
	}

	return false;
}


bool probe(){ // Updates the ball position and returns true only if the probe phase is done

	if(ProbePhase){

		//Check for distance with table edge Probe
		if(!ProbeBallEdge){ // Is true until the end of the Probe phase
			float distanceBetween_x = (probeDistance/2) - probeCenter_x; 
			std:: cout<< distanceBetween_x<< std::endl;
			if (distanceBetween_x <= 0){  
				// when the ball hits edge
				ProbeBallEdge = true;
				lastTimeProbe = elapsed;
				lastFrameProbe = frameN + (Probe2CueDelay/11.76); //includes delay frames
				}
			if(!ProbeBallEdge){//update probe movement 
				probeCenter_x += probeSpeed;
			}
		}

		return ProbeBallEdge; // True only if the end of the probe 
	}

	return false;
}

bool sleep(){
	for(int i = Probe2CueDelay; i>0; --i){
				std:: cout<< i<< std::endl;
			}
	return true;
}
void online_trial()
{

	// while the experiment is running
	if (fingersCalibrated && !finished)
	{	
		
		
		// set the cue velocity after striking, this only happens once per trial.
		if(!cueVelSet){
			timeStart = elapsed;
			cueVelSet = true;
		}

		if(Order == 1){
			
			if(!ProbePhase){
				ProbePhase = cue();
			}else if(ProbePhase  && elapsed > lastFrameCue + Probe2CueDelay){
			probe();
			}
			
		} else if(Order == 2){
			probe();
			if (ProbeBallEdge && elapsed > lastTimeProbe + Probe2CueDelay){
				ProbePhase = false;
				cue();
			}
	
		}
		// Advance frame number
		frameN++;
	}
}



///////////////////////////////////////////////////////////
/////// USUALLY DON'T NEED TO EDIT THESE FUNCTIONS ////////
///////////////////////////////////////////////////////////

void updateTheMarkers()
{
	optotrak.updateMarkers();
	markers = optotrak.getAllMarkers();
}

void initVariables() 
{
	trial.init(parameters);
	interoculardistance = str2num<double>(parameters.find("IOD"));
}

void update(int value)
{
	glutPostRedisplay();
	glutTimerFunc(TIMER_MS, update, 0);
}

void handleResize(int w, int h)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0,0,SCREEN_WIDTH, SCREEN_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}

void initProjectionScreen(double _focalDist, const Affine3d &_transformation, bool synchronous)
{
	focalDistance = _focalDist;	
	screen.setWidthHeight(SCREEN_HIGH_SIZE*SCREEN_WIDTH/SCREEN_HEIGHT, SCREEN_HIGH_SIZE);//(SCREEN_WIDE_SIZE, SCREEN_WIDE_SIZE*SCREEN_HEIGHT/SCREEN_WIDTH);
	screen.setOffset(alignmentX,alignmentY);
	screen.setFocalDistance(_focalDist);
	screen.transform(_transformation);
	cam.init(screen);
	if ( synchronous )
		moveScreenAbsolute(_focalDist,homeFocalDistance,4500);
	else
		moveScreenAbsoluteAsynchronous(_focalDist,homeFocalDistance,4500);
}

void initRendering()
{   
	// Clear buffers
	glClearColor(0.0,0.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* Set depth buffer clear value */
	glClearDepth(1.0);

	/* Enable depth test */
	glEnable(GL_DEPTH_TEST);

	// Not sure...for meshes which faces are away from camera
	//glEnable(GL_CULL_FACE);

	/* Set depth function */
	glDepthFunc(GL_LEQUAL);

	// Nice perspective calculations
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	// Set up the lighting
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	//glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	//glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);

	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, LightSpecular);
	glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
	//glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.3f);
	glEnable(GL_LIGHT1);

	// Clean modelview matrix to start
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void initMotors()
{
	homeEverything(6000,4000);
}

void initOptotrak()
{
	optotrak.setTranslation(calibration);

	if ( optotrak.init(LastAlignedFile, OPTO_NUM_MARKERS, OPTO_FRAMERATE, OPTO_MARKER_FREQ, OPTO_DUTY_CYCLE,OPTO_VOLTAGE) != 0)
	{   cerr << "Something during Optotrak initialization failed, press ENTER to continue. A error log has been generated, look \"opto.err\" in this folder" << endl;
	cin.ignore(1E6,'\n');
	exit(0);
	}

	// Read 10 frames of coordinates and fill the markers vector
	for (int i=0; i<10; i++)
	{
		updateTheMarkers();
	}
}

void cleanup()
{
	// Stop the optotrak
	optotrak.stopCollection();
}

void beepOk(int tone)
{
#ifndef SIMULATION
	switch(tone)
	{
	case 0:
		// Remember to put double slash \\ to specify directories!!!
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-1.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 1:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\calibrate.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 2:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-8.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 3:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-reject.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 4:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-twoBlips.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 7:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-left.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 8:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-right.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 9:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-home.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 10:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-grasp.wav", 
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 11:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-marker.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 12:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-estimate.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 13:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-8_lowpass.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 14:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-8_double.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 15:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-rising.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 16:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-falling.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 17:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-highBubblePop.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 18:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-lowBubblePop.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 19:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-success.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	case 20:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\spoken-watch.wav",
			NULL, SND_FILENAME | SND_ASYNC);
		break;
	}
#endif
	return;
}

///////////////////////////////////////////////////////////
////////////////////// MAIN FUNCTION //////////////////////
///////////////////////////////////////////////////////////

int main(int argc, char*argv[])
{
	mathcommon::randomizeStart();

	// Initializes the optotrak and starts the collection of points in background
	initMotors();
	initOptotrak();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
	glutGameModeString(GAME_MODE_STRING);
	glutEnterGameMode();
	//glutFullScreen();

	initRendering();
	initStreams(); // parameters file is loaded

	initVariables(); // staircases are built
	glutDisplayFunc(drawGLScene);
	glutKeyboardFunc(handleKeypress);
	glutReshapeFunc(handleResize);
	glutIdleFunc(idle);
	glutTimerFunc(TIMER_MS, update, 0);
	glutSetCursor(GLUT_CURSOR_NONE);

	//boost::thread initVariablesThread(&initVariables);

	glutMainLoop();

	homeEverything(6000,4000);
	cleanup();
	return 0;
}
