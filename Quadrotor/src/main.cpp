#define _CRT_RAND_S


#define X_DIM 13  //
#define U_DIM 4  //
#define Z_DIM 13  //  
#define INFTY 9e9

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <ppl.h>
#include "matrix.h"
#include "utils.h"
#include "callisto.h"
#include "Dynamics.h"
#include "RRT.h"

typedef Matrix<X_DIM> State;
typedef Matrix<U_DIM> Input;
typedef Matrix<Z_DIM> Observation;
typedef Matrix<3,3> Rotation;

//world constants
double dt;

//derivative constants;
Matrix<3,3> invInertia;
double nominalInput;

double QuadrotorR;

//callisto visualization
int cal_quadrotor;
int cal_scene;
int cal_floor;
int cal_camera1;
int cal_camera2;
int cal_viewframe1;
int cal_viewframe2;
int cal_ellipse;
int cal_path;
int cal_sphere;
int cal_goal;
int cal_obstacles;
int cal_environment;
int cal_point;
int cal_box;
int cal_rrt;
Matrix<3,1> goal;
Matrix<X_DIM, 1> start;
Matrix<3,3> Rotinitial;
double goal_radius;

void initEnvironment()
{
	double tmpdt = 0.1;
	Dynamics dy(tmpdt);


	CAL_SetViewParams(0, 2.5, -18, 5, 2.5, 0, 5, 0, 0, 1); 
	goal[0] = 4; goal[1] = 8; goal[2] = 4;
	goal_radius = 0.5;
	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 1, 0, 0);
	CAL_CreateGroup(&cal_box, 0, true, "cal_box");
	CAL_SetGroupColor(cal_box, 0,0,0, 0.6);
	CAL_CreateGroup(&cal_goal, 0, true, "cal_goal");
	CAL_SetGroupColor(cal_goal, 0,0,0, 0.6);
	CAL_CreateGroup(&cal_quadrotor, 0, true, "cal_quadrotor");
	CAL_SetGroupColor(cal_quadrotor, 1.0, 0.0, 0.0, 1);
	CAL_CreateGroup(&cal_rrt, 0, true, "cal_rrt");
	CAL_SetGroupColor(cal_rrt, 1,0,0, 0.6);
	QuadrotorR = 0.35;

	start = zeros<X_DIM,1>();
	start[0] = 0; start[1] = 0; start[2] = 0; start[3] = 0.0; start[4] = 0.0; start[5] = 0.0; start[6] = 0; start[7] = 0; start[8] = 0; 
	start[9] = 0.0; start[10] = 0.0; start[11] = 0.0; start[12] = dy.nominalInput;
	
	Rotinitial = identity<3>();

	//CAL_CreateSphere(cal_quadrotor, QuadrotorR, start[0], start[1], start[2]);


	int np[1] = {2};
	Matrix<3,1> bbmin, bbmax;
	bbmin[0] = 0; bbmin[1] = 0; bbmin[2] = 0;
	bbmax[0] = 5; bbmax[1] = 10; bbmax[2] = 5;
	
	float side1[6] = {(float)bbmin[0], (float)bbmin[1], (float)bbmin[2], (float)bbmax[0], (float)bbmin[1], (float)bbmin[2]};
	float side2[6] = {(float)bbmin[0], (float)bbmin[1], (float)bbmin[2], (float)bbmin[0], (float)bbmax[1], (float)bbmin[2]};
	float side3[6] = {(float)bbmin[0], (float)bbmin[1], (float)bbmin[2], (float)bbmin[0], (float)bbmin[1], (float)bbmax[2]};

	float side4[6] = {(float)bbmax[0], (float)bbmax[1], (float)bbmax[2], (float)bbmin[0], (float)bbmax[1], (float)bbmax[2]};
	float side5[6] = {(float)bbmax[0], (float)bbmax[1], (float)bbmax[2], (float)bbmax[0], (float)bbmin[1], (float)bbmax[2]};
	float side6[6] = {(float)bbmax[0], (float)bbmax[1], (float)bbmax[2], (float)bbmax[0], (float)bbmax[1], (float)bbmin[2]};
	
	float side7[6] = {(float)bbmin[0], (float)bbmax[1], (float)bbmax[2], (float)bbmin[0], (float)bbmin[1], (float)bbmax[2]};
	float side8[6] = {(float)bbmin[0], (float)bbmax[1], (float)bbmax[2], (float)bbmin[0], (float)bbmax[1], (float)bbmin[2]};
	float side9[6] = {(float)bbmax[0], (float)bbmin[1], (float)bbmax[2], (float)bbmin[0], (float)bbmin[1], (float)bbmax[2]};
	
	float side10[6] = {(float)bbmax[0], (float)bbmin[1], (float)bbmax[2], (float)bbmax[0], (float)bbmin[1], (float)bbmin[2]};
	float side11[6] = {(float)bbmax[0], (float)bbmax[1], (float)bbmin[2], (float)bbmin[0], (float)bbmax[1], (float)bbmin[2]};
	float side12[6] = {(float)bbmax[0], (float)bbmax[1], (float)bbmin[2], (float)bbmax[0], (float)bbmin[1], (float)bbmin[2]};

	CAL_CreateGroup(&cal_box, 0, true, "cal_box");
	CAL_SetGroupColor(cal_box, 0,0,0, 0.6);
	CAL_CreatePolyline(cal_box, 1, np, side1);
	CAL_CreatePolyline(cal_box, 1, np, side2);
	CAL_CreatePolyline(cal_box, 1, np, side3);
	CAL_CreatePolyline(cal_box, 1, np, side4);
	CAL_CreatePolyline(cal_box, 1, np, side5);
	CAL_CreatePolyline(cal_box, 1, np, side6);
	CAL_CreatePolyline(cal_box, 1, np, side7);
	CAL_CreatePolyline(cal_box, 1, np, side8);
	CAL_CreatePolyline(cal_box, 1, np, side9);
	CAL_CreatePolyline(cal_box, 1, np, side10);
	CAL_CreatePolyline(cal_box, 1, np, side11);
	CAL_CreatePolyline(cal_box, 1, np, side12);

	int obj_id;
	CAL_CreateSphere(cal_goal, goal_radius, (float)goal[0], (float)goal[1], (float)goal[2], &obj_id);


	int obj;
	// Visualization parameters
	double beamWidth     = 0.015; // m
	double beamHeight    = 0.0065; // m
	double beamRadius    = 0.02; // m
	double motorRadius   = 0.015; // m
	double motorHeight   = 0.02; // m
	double rotorRadius   = 0.10; // m
	double rotorHeight   = 0.005; // m
	double centerSide    = 0.0889; // m
	double centerHeight  = 0.0365; // m
	double centerTopSide = 0.03; // m
	double flagLength    = 0.0508; // m
	double tileSize      = 1;  // m
	double length      = 0.3429/2;  
	// Quadrotor
	CAL_CreateGroup(&cal_quadrotor, 0, false, "QuadRotor");
	CAL_SetGroupColor(cal_quadrotor, 0.05, 0.05, 0.05);
	CAL_CreateBox(cal_quadrotor, 2*length, beamWidth, beamHeight, 0, 0, 0);
	CAL_CreateBox(cal_quadrotor, beamWidth, 2*length, beamHeight, 0, 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, length, 0, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, -length, 0, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, length, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, -length, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, length, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, -length, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, length, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, -length, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, -length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, -length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateBox(cal_quadrotor, centerSide, centerSide, beamHeight, 0, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, 0, 0, (float) (M_PI*0.25));
	CAL_CreateBox(cal_quadrotor, flagLength, beamWidth + 0.001, beamHeight + 0.001, length / 1.65, 0, 0, &obj);
	CAL_SetObjectColor(obj, 1, 0.15, 0);

	float flagTriangle[18] = {length / 1.65 - flagLength / 2, 0, -beamHeight / 2,
		length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
		length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
		length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
		length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
		length / 1.65 - flagLength / 2, 0, -beamHeight / 2};
	CAL_CreateTriangles(cal_quadrotor, 2, flagTriangle, &obj);
	CAL_SetObjectColor(obj, 1, 0.15, 0);

	float polygon1[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight,
		sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
		sqrt(2.0)*centerSide/2, 0, 0,
		sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon1, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon2[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight,
		sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
		sqrt(2.0)*centerSide/2, 0, 0,
		sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon2, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon3[18] = {0, -sqrt(2.0)*centerSide/2, 0,
		0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight,
		0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
		0, sqrt(2.0)*centerSide/2, 0,
		0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
		0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon3, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon4[18] = {0, -sqrt(2.0)*centerSide/2, 0,
		0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight,
		0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
		0, sqrt(2.0)*centerSide/2, 0,
		0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
		0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon4, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	double ballRadius  = 0.034528125;  
	CAL_CreateSphere(cal_quadrotor, ballRadius, 0, 0, centerHeight + ballRadius, &obj);
	CAL_SetObjectColor(obj, 1, 1, 0);

	moveQuadrotor(cal_quadrotor, start, Rotinitial);

}



int main(int argc, char *argv[])
{
	srand((unsigned int) time(NULL));
	CAL_Initialisation (true, true, true);
	initEnvironment();

	dt = 0.25;
	double plantime = 2.0;


	Matrix<X_DIM,X_DIM> P0 = 0.0001 * 0.0001 * identity<X_DIM>(); // initial state variance
	P0(0,0) = P0(1,1) = 0.1*0.1; P0(2,2) = 0.0001*0.0001; // position (z very certain, because on floor)
	P0(3,3) = P0(4,4) = P0(5,5) = 0.0001*0.0001; // velocity (very certain, because on floor)
	P0(6,6) = P0(7,7) = P0(8,8) = 0.0001*0.0001; // orientation (zero by definition)
	P0(9,9) = P0(10,10) = P0(11,11) = 0.0001*0.0001; // angular speed (very certain, because on floor)
	P0(12,12) = 0.1*0.1; // force

	//RRT rrt(start, Rotinitial, goal, dt, plantime, goal_radius, cal_obstacles);
	RRT rrt(start, Rotinitial, goal, dt, plantime, goal_radius, cal_obstacles);
	rrt.setPlannerDefaults();
	rrt.Plan_K_Seconds();
	rrt.showPath(cal_rrt);


	int num;
	std::cin>>num;

	CAL_End();
}
