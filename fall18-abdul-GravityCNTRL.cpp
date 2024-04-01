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
BalanceFactor<double> trial;

// FOR NEW CALIBRATION:
double markerXOffset = 10;
double markerYOffset = 10;

// Variables for counting trials, frames, and lost frames
int trialNumber = 0;
int frameN = 0;
int trialsPerBlock = 36;

// Flags for important states in the experiment
bool fingersCalibrated = true;
bool cueBallFalls = false;
bool floorTouch = false;
bool cueBallStruck = false;
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

//Table
float Tablex1 = -100;
float Tablex2 = 100; 
float Tabley1 = -70;
float TableZ1 = displayDepth- ;
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
//Noise
//float NoiseX1 = Tablex2;
//float NoiseX2 = 150;
//float NoiseY1 = Legy1;
//float NoiseY2 = Legy2;

//Shot Glass
float probePos;
float shotRadius = 14;
float shotHeight = 20;


//Physics & misc
float Gravity;
double speed;
double timeOfImpact;
double timeOfFall;
double impact_z;
double fallDuration;
double frameOfFall;
double timeStart;
double dist;
double lastFrame;
double responseGravity; //Subject response to compare to 9.81
double responseVal; 

//NOISE

const int numTilesX = 50;
const int numTilesY = 50;
const int numTiles = numTilesX * numTilesY;
float tileSize = 2;

float mask1_x1[numTiles];
float mask1_x2[numTiles];
float mask1_y1[numTiles];
float mask1_y2[numTiles];

float mask2_x1[numTiles];
float mask2_x2[numTiles];
float mask2_y1[numTiles];
float mask2_y2[numTiles];

float mask1_colors[numTiles];
float mask2_colors[numTiles];
float maskX_offset = 50, maskY_offset = -75;
///////////////////////////////////////
float responseMag;
float responseDif_x;
float responseDif_z;
float cPointX,cPointY,cPointZ;
float movementX = 0;
float movementY = 0;
float movementZ = 0;
double lineVelMagnitude;
double responseDelay = 1500; 
double timeOfContact;
float distanceBetween; 
float pointRotation;
double elapsed;

int cueBallPosition;
double aimingAngle;
double deflectionAngle;
double cueZX_ratio;
double targetZX_ratio;

double vecLength_unnorm;


double distanceToCueBall = 999;

bool handAtStart;
float distanceToStart;

float beta = DEG2RAD*15; // this determines where the target points are
float gamma = DEG2RAD*8; // this determines the deflectionAngle of bad physics
float point_x;
float point_z;
float delta_x2;
float delta_z2;
float alpha2;
float point_x_cue;
float point_z_cue;
float point_x_cue2;
float point_z_cue2;

int response;


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
void online_apparatus_alignment();
void online_fingers();
void online_trial();

/*************************** EXPERIMENT SPECS ****************************/
// experiment directory
string experiment_directory = "R:/CLPS_Domini_Lab/abdul/fall18-GravityCNTRL/";
// parameters file directory and name
string parametersFile_directory = experiment_directory + "fall18-GravityEXP1Parameters.txt";
// trial file headers
string trialFile_headers = "subjName\ttrialN\tGravity\tspeed\telapsed\tframeN\tcueBallFalls\ttimeOfFall\tframeOfFall\ttimeOfImpact\tlastFrame\tprobePos\tballPos_y\tballPos_z\timpact_z"; 
/*************************** FUNCTIONS ***********************************/
// First, make sure the filenames in here are correct and that the folders exist.
// If you mess this up, data may not be recorded!
void initStreams()
{
	ifstream parametersFile;
	parametersFile.open(parametersFile_directory.c_str());
	parameters.loadParameterFile(parametersFile);

	string subjectName = parameters.find("SubjectName");
	
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
				// calibration on the X
				indexCalibrationPoint=markers.at(ind0).p;
				indexCalibrationPoint[0] = indexCalibrationPoint[0] - 25;
				indexCoords.init(indexCalibrationPoint, markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );

				fingerCalibrationDone=3;
				fingersCalibrated=true;
				visibleInfo=false;
				beepOk(0);
				trialNumber++;
				trial.next();
				initTrial();
				break;
			}
		break;
		

		case '8': // higher gravity
		{
			if (probePos >= TableZ1 + shotRadius*1.5 ){
		  probePos = probePos - 3;
			
			}
			else {
				beepOk(3);	
			}
			
		}
		break;

		case '2': // lower gravity
		{
			if (probePos <= TableZ1 + 170 ){
		  probePos = probePos + 3;
			
			}
			else {
				beepOk(3);	
			}
			
		}
		break;

		case '+':
		{
			if ((elapsed > timeOfImpact + responseDelay) && cueBallFalls && !(response == 2)){
				response++;
				if (response == 2){
					advanceTrial();
					beepOk(19);
				}
				else { 
					beepOk(17);
				}
			}
		}
		
		break;
		
	}
}

/*** GRASP ***/
void calibration_fingers(int phase)
{
	/*
	switch (phase)
	{
		case 1:
		{
			if(isVisible(markers[calibration1].p) && isVisible(markers[calibration2].p))
			{
				indexCalibrationPoint=markers.at(calibration1).p;
				indexCalibrationPoint[0] = indexCalibrationPoint[0] - markerXOffset;
				indexCalibrationPoint[1] = indexCalibrationPoint[1] + markerYOffset;
				thumbCalibrationPoint=markers.at(calibration2).p;
				thumbCalibrationPoint[0] = thumbCalibrationPoint[0] - markerXOffset;
				thumbCalibrationPoint[1] = thumbCalibrationPoint[1] - markerYOffset;
			}
		} break;
		case 2:
		{
			indexCoords.init(indexCalibrationPoint, markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );
			thumbCoords.init(thumbCalibrationPoint, markers.at(thu1).p, markers.at(thu2).p, markers.at(thu3).p );
		} break; 
		case 3:
		{
			indexJointCoords.init(indexCalibrationPoint, markers.at(ind1).p, markers.at(ind2).p, markers.at(ind3).p );
			thumbJointCoords.init(thumbCalibrationPoint, markers.at(thu1).p, markers.at(thu2).p, markers.at(thu3).p );
		} break;
		
	}
	*/
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
			text.draw("# Name: " +parameters.find("SubjectName"));
			text.draw("# IOD: " +stringify<double>(interoculardistance));
			text.draw("# Finger to Cue Distance: " +stringify<double>(distanceToCueBall));
			text.draw("# Gravity: " +stringify<float>(Gravity));
			text.draw("# Horizontal Velocity:" +stringify<double>(speed));
			text.draw("# floorTouch:" +stringify<double>(floorTouch));
			text.draw("# cueBallFalls? " +stringify<double>(cueBallFalls));
			text.draw("# ballPos_z " +stringify<float>(ballPos_z));
			text.draw("# frameOfFall " +stringify<double>(frameOfFall));
			text.draw("# lastFrame " +stringify<double>(lastFrame));
			text.draw("# Z-position: " +stringify<float>(cueCenter_z));
			text.draw("# Y-position: " +stringify<float>(cueCenter_y));
			text.draw("# probePos:" +stringify<float>(probePos));
			
			
			
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
    else
    {   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0,0.0,0.0,1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        cam.setEye(eyeRight);
        drawStimulus();
		drawInfo();
        glutSwapBuffers();
    }
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
	
	
	// 4. Draw target ball
	if(!floorTouch){
	
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

	/*5. Draw noise
	 glBegin(GL_QUADS);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, noiseMaterial);
	glVertex3f(NoiseX1, NoiseY2,displayDepth-50); // vertex 1
	glVertex3f(NoiseX1, NoiseY1,displayDepth-50); // vertex 2
	glVertex3f(NoiseX2, NoiseY1,displayDepth-50); // vertex 3
	glVertex3f(NoiseX2, NoiseY2,displayDepth-50); // vertex 4
	glEnd(); */
	/*if((elapsed < timeOfImpact) && cueBallFalls){
		for(int i=0; i<numTiles; i++){
			glColor3f(mask2_colors[i], 0.0f, 0.0f);
			glBegin(GL_POLYGON); 
				glVertex3f(mask2_x1[i],mask2_y1[i],displayDepth+50); 
				glVertex3f(mask2_x2[i],mask2_y1[i],displayDepth+50); 
				glVertex3f(mask2_x2[i],mask2_y2[i],displayDepth+50); 
				glVertex3f(mask2_x1[i],mask2_y2[i],displayDepth+50); 
			glEnd();
		}
	}*/

	
	if((elapsed > lastFrame + responseDelay) && floorTouch){//  after display period
	
	// 5. Draw response point
		
		glPushMatrix();
		glLoadIdentity();
		GLUquadricObj* qobj = gluNewQuadric();
		gluQuadricDrawStyle(qobj, GLU_FILL);
		glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE, lineMaterial);
		glTranslated(cueCenter_x,-137.83,probePos); 
		gluSphere(qobj, 2, 6, 6);
		gluDeleteQuadric(qobj);
		glPopMatrix();
	}
	/*glPushMatrix();
	glLoadIdentity();
	GLUquadricObj* shotGlass = gluNewQuadric();
	gluQuadricDrawStyle(shotGlass, GLU_FILL);
	glTranslated(cueCenter_x,-145,probePos);
	glRotatef(-90, 1.0f, 0.0f, 0.0f); 
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tableMaterial);
	gluCylinder(shotGlass, cueRadius,shotRadius,shotHeight, 64, 64);
	gluDeleteQuadric(shotGlass);
	glPopMatrix();
	glEnd();
	
	
	 draw curve
	glPushMatrix();
	glLoadIdentity();
	glTranslated(.3*cueRadius+speed,Tabley1+cueRadius ,displayDepth); 
	glLineWidth(3.0f);
	double step_size = .1; 
	glBegin(GL_LINE_STRIP);
	glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,lineMaterial);
	for (double y = 0; y >= -100; y = y - step_size){ 
		glVertex3d(sqrt(-1*y*responseVal),y,0); 
	}
	glEnd();
	glPopMatrix();*/
	}


void build_masks()
{
	/*int tileNum = 0;

	for (int xpos=0-(numTilesX/2); xpos<(numTilesX/2); xpos++){
		for (int ypos=0-(numTilesY/2); ypos<(numTilesY/2); ypos++){

			double randNum = unifRand(0,2);
			if (randNum>2)
				mask1_colors[tileNum] = 1.0;
			if (randNum>1.0 && randNum<=2.0)
				mask1_colors[tileNum] = 0.75;
			if (randNum<=1.0)
				mask1_colors[tileNum] = 0.0;

			mask1_x1[tileNum] = xpos*tileSize+maskX_offset;
			mask1_x2[tileNum] = (xpos+1)*tileSize+maskX_offset;
			mask1_y1[tileNum] = ypos*tileSize+maskY_offset;
			mask1_y2[tileNum] = (ypos+1)*tileSize+maskY_offset;
			
			tileNum++;
		}
	}*/

	int tileNum = 0;

	for (int xpos=0-(numTilesX/2); xpos<(numTilesX/2); xpos++){
		for (int ypos=0-(numTilesY/2); ypos<(numTilesY/2); ypos++){

			double randNum = unifRand(0,2);
			if (randNum>2)
				mask2_colors[tileNum] = 1.0;
			if (randNum>1.0 && randNum<=2.0)
				mask2_colors[tileNum] = 0.75;
			if (randNum<=1.0)
				mask2_colors[tileNum] = 0.0;

			mask2_x1[tileNum] = xpos*tileSize+maskX_offset;
			mask2_x2[tileNum] = (xpos+1)*tileSize+maskX_offset;
			mask2_y1[tileNum] = ypos*tileSize+maskY_offset;
			mask2_y2[tileNum] = (ypos+1)*tileSize+maskY_offset;
			
			tileNum++;
		}
	}
}

// called at the beginning of every trial
void initTrial()
{
	// initializing all variables
	frameN=0;
	cueVelSet = false;
	cueBallFalls = false;
	floorTouch = false;
	timeOfImpact = 0;
	timeOfFall = 0;
	ballPos_y = 0;
	ballPos_z = 0;
	probePos = probePos = -1*(rand()% 150+450); 
	response = 0;
	ballMaterial[3] = 1;
	dist = 0;
	build_masks();

	Gravity = trial.getCurrent()["Gravity"];
	speed = trial.getCurrent()["Speed"];
	cueCenter_x = 0;
	cueCenter_y = -62;
	cueCenter_z = ballStartPos_z; 
	
	
	initProjectionScreen(displayDepth);
	
	// roll on
	drawGLScene();
	timer.start();
}

// This function handles the transition from the end of one trial to the beginning of the next.
void advanceTrial() {
	
	trialFile << fixed <<
	parameters.find("SubjectName") << "\t" <<		//subjName
	trialNumber << "\t" <<							//trialN
	Gravity << "\t" <<
	speed << "\t" <<
	elapsed << "\t" <<
	frameN << "\t" <<
	cueBallFalls << "\t" <<
	timeOfFall << "\t" <<
	frameOfFall << "\t" <<
	timeOfImpact << "\t" <<
	lastFrame << "\t" <<
	probePos << "\t" <<
	ballPos_y << "\t" <<
	ballPos_z << "\t" <<
	impact_z << endl;

	//if(trialFile.is_open())
	//	trialFile.close();

	if(!trial.isEmpty()){
		trial.next();
		trialNumber++;
		initTrial();
	}else{
		trialFile.close();
		finished=true;
	}
}

void idle() {

	elapsed = timer.getElapsedTimeInMilliSec();

	// get new marker positions from optotrak
	updateTheMarkers();

	// eye coordinates
	eyeRight = Vector3d(interoculardistance/2,0,0);//0
	eyeLeft = Vector3d(-interoculardistance/2,0,0);//0

	/* Write to trialFile once calibration is over
	if (fingersCalibrated) // write every frame if grasping
	{
	trialFile << fixed <<
	parameters.find("SubjectName") << "\t" <<		//subjName
	trialNumber << "\t" <<							//trialN
	Gravity << "\t" <<
	speed << "\t" <<
	elapsed << "\t" <<
	frameN << "\t" <<
	cueVelSet << "\t" <<
	cueBallFalls << "\t" <<
	cueCenter_x << "\t" <<
	cueCenter_y << "\t" <<
	cueCenter_z << "\t" <<
	fallDuration << "\t" <<
	frameOfFall << "\t" <<
	timeOfImpact << "\t" <<
	ind.transpose() << "\t" <<
	responseGravity << endl;
	}*/
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

void online_trial()
{

	// while the experiment is running
	if (fingersCalibrated && !finished && !floorTouch)
	{
		

		/* in the active phase, get the marker-based index finger position
		if(activePhase){
			indXNow = ind.x();
			indYNow = ind.y();
			indZNow = ind.z();
		// in the passive phase, get the remembered trajectories XX
		}else{
			if(!cueBallStruck ){ // (only when the cueBall isn't yet struck)
					indXNow = indXNow + movementX;
					indYNow = indYNow + movementY ;
					indZNow = indZNow + movementZ ;
			}
		}

		// test if behind virtual wall at 280 mm
		distanceToStart = -280 - indZNow;
		handAtStart = distanceToStart<0;
		*/

		// set the cue velocity after striking, this only happens once per trial.
		if(!cueVelSet){
			timeStart = elapsed;
			cueVelSet = true;
		}
		//Check for contact with table surface
		if(!cueBallFalls){
			float distanceBetween_x = Tablex2 - cueCenter_x; 
			float distanceBetween_y = Tabley1 - cueCenter_y ;
			float distanceBetween_z = TableZ1 - cueCenter_z;
			if (!cueBallFalls && (distanceBetween_z <=(-.3*cueRadius))){ 
			// when the ball falls
				frameOfFall = frameN+1;
				timeOfFall = elapsed;
				cueBallFalls = true;
				floorTouch = false;
				
			}
		}
	
		//Check for contact with floor surface
		if(!floorTouch){
			std::cout << frameN << "  " << cueCenter_y << std::endl;
			//std::cout << cueCenter_z << std::endl;
			  if (cueBallFalls && !floorTouch && (cueCenter_y <= -129.83)){ //XX
				// when the ball touches the floor
					cueCenter_z += 0;
					cueCenter_y += 0;
					ballPos_z = cueCenter_z;
					ballPos_y = cueCenter_y;
					lastFrame = frameN;
					timeOfImpact = elapsed; 
					floorTouch = true;
			    
		  }
		}
		
		// update ball positions 
		if(cueBallFalls && !floorTouch){
			double oldY = cueCenter_y;
			cueCenter_z += speed;
			cueCenter_y = (Tabley1+cueRadius) - 0.5*(Gravity/1000)*pow(11.76*(frameN-frameOfFall+1), 2); //starts at 
			if (cueCenter_y <= -129.83)
			{
				cueCenter_y = -129.83;
				cueCenter_z = (speed/11.76)* sqrt(75*(2000/Gravity)) -600;
			}
			/*if (frameN > frameOfFall + 9){ //ONLY SHOWS LAST 3 FRAMES

				ballMaterial[3] =  alpha;
			}*/
			//dist += sqrt(pow(speed,2) +pow(cueCenter_y- oldY,2));
			//ballMaterial[3] =  alphaVal - alphaVal*(dist/occlusion);
			
		}
		if (!cueBallFalls){
			cueCenter_z += speed;
			cueCenter_y += 0;
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
	/*for(int d=0; d<360; d++){
		// Frontoparallel circle at display depth
		goalX[d] = cos(DEG2RAD*d)*targetRadius;
		goalY[d] = sin(DEG2RAD*d)*targetRadius;
		goalZ[d] = displayDepth-goalDepth;
	}*/
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
