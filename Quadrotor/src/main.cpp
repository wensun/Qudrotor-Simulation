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

typedef Matrix<X_DIM> State;
typedef Matrix<U_DIM> Input;
typedef Matrix<Z_DIM> Observation;
typedef Matrix<3,3> Rotation;

//world constants
double dt;

//derivative constants;
Matrix<3,3> invInertia;
double nominalInput;

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

Matrix<3,1> goal;
Matrix<X_DIM, 1> start;

void initEnvironment()
{

	CAL_SetViewParams(0, 0, 0, 14, 0, 0, 0, 0, 1, 0); 
	goal[0] = 0; goal[1] = 1; goal[2] = -7;
	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 1, 0, 0);
	CAL_CreateGroup(&cal_box, 0, true, "cal_box");
	CAL_SetGroupColor(cal_box, 0,0,0, 0.6);
	CAL_CreateGroup(&cal_goal, 0, true, "cal_goal");
	CAL_SetGroupColor(cal_goal, 0,0,0, 0.6);


	int np[1] = {2};
	Matrix<3,1> bbmin, bbmax;
	bbmin[0] = -5; bbmin[1] = -5; bbmin[2] = -8;
	bbmax[0] = 5; bbmax[1] = 5; bbmax[2] = 8;
	
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
	CAL_CreateSphere(cal_goal, 0.2, (float)goal[0], (float)goal[1], (float)goal[2], &obj_id);

}



int main(int argc, char *argv[])
{
	srand((unsigned int) time(NULL));
	CAL_Initialisation (true, true, true);
	initEnvironment();


	int num;
	std::cin>>num;

	CAL_End();
}
