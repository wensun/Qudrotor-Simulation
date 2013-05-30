#define _CRT_RAND_S

#define DIM 2
#define X_DIM 4  //X = (x,y,theta, v)  theta: orientation, v speed.
#define U_DIM 2  //a, phi, a: acceleration, phi: steering wheel angel
#define Z_DIM 3  //it can observe position x and position y.  
#define INFTY 9e9

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <ppl.h>

#include "RRTNode.h"
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Controller.h"
#include "LQGMP.h"

/***********************Environment settings***************************/
int cal_environment;
int cal_obstacles;
int cal_goal;
int cal_rrt;
int cal_paths;
int cal_ellipse;
int cal_point;
int cal_cienvironment, cal_ciobstacles, cal_cipoint;

Matrix<U_DIM> u_min, u_max;  //control limits. acceleration and wheel angel
Matrix<X_DIM> x_min, x_max;  //the bounds of the environment
Matrix<X_DIM> start;         //start position 4-D
Matrix<X_DIM,X_DIM> P0; //start covariance of distribution.
Matrix<DIM> goal;            //goal position  2-D, just the x y position
double goalRadius;             
double dt;
double car_l;

/******************End of Environment settings************************/


/*******************RRT stuff setting**************************/
std::vector<std::vector<RRTNode>> PathSet;
std::vector<RRTNode> rrttree;
std::vector<int> paths;

int maxchildren = 20;
int maxattempts = 10;
double plantime = 0.5;
/*****************End of RRT stuff setting**********************/

double Random()
{
	double a;
	a = rand()%1000;
	return a / 1000;
}

/**********************init environment***************************/
void initEnvironment() 
{
	goalRadius = 0.5;
	dt = 0.5;

	// car length
	car_l = 1;

	//range of x: -5 ~ 5 
	x_min[0] = -5; 
	x_max[0] = 5;
	//range of y: -5 ~ 5
	x_min[1] = -5; 
	x_max[1] = 5;	
	//range of theta: the orientation of the car 0:360.
	x_min[2] = 0;
	x_max[2] = 2*M_PI;
	//range of speed: -1:1
	x_min[3] = 0;
	x_max[3] = 1.0;
	//range of acceleration: -0.5~0.5 
	u_min[0] = -0.4; 
	u_max[0] = 0.8;
	//range of wheel angel : -45 ~ 45 degree.
	u_min[1] = -M_PI*0.35;
	u_max[1] = M_PI*0.35;


	// callisto params
	CAL_SetViewParams(0, 0, 0, 14, 0, 0, 0, 0, 1, 0); 

	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 1, 0, 1);

	CAL_CreateGroup(&cal_paths, 0, false, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.7f);


	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 1, 0, 0);

	CAL_CreateBox(cal_obstacles, 10.2f, 0.25f, 0.2f, 0, -5.f, -0.1f);
	CAL_CreateBox(cal_obstacles, 0.25f, 10.2f, 0.2f, -5.f, 0, -0.1f);
	CAL_CreateBox(cal_obstacles, 0.25f, 10.2f, 0.2f, 5.f, 0, -0.1f);
	CAL_CreateBox(cal_obstacles, 10.2f, 0.25f, 0.2f, 0, 5.f, -0.1f);
	
	//CAL_CreateBox(cal_obstacles, 2.f, 2.f, 0.05f, -2.75f, -2.75f, 0);
	CAL_CreateBox(cal_obstacles, 3.f, 0.6f, 0.05f, -2.f, 2.75f, 0);
	CAL_CreateBox(cal_obstacles, 2.f, 2.f, 0.05f, 3.f, -2.5f, 0);
	//CAL_CreateBox(cal_obstacles, 2.f, 4.4f, 0.05f, 1.5f, 1.8f, 0);

	CAL_CreateBox(cal_obstacles, 1.5f, 0.6f, 0.05f, 3.f, 2.5f, 0);

	//start as (-4.5, 0.3, 0, 0). orientation:0, velocity:0
	start[0] = -4.5; start[1] = 0; start[2] = 0; start[3] = 0.1;
	P0 = identity<X_DIM>() * 0.001;
	//goal center: (4.25, 4.25)
	goal[0] = 4.25; goal[1] = -4.25;

	
	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	CAL_SetGroupColor(cal_goal, 0, 1, 1, 0.5);
	CAL_CreateCylinder(cal_goal, (float) goalRadius, 0.05f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.025f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	
}

/*********************End of environment init**************************/

/************All RRT functions*****************/

////////////////////Distance between two RRT nodes/////////////////
double dist(const Matrix<X_DIM>& x1, const Matrix<X_DIM>& x2)
{
	double distance = 0;

	//calculate the Eucliead distance x^2+y^2
	double d = (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) + (x1[3]-x2[3])*(x1[3]-x2[3]);
	distance = sqrt(d);
	return distance;
}
///////////////////////////////////////////////////////////////////

//////////////////////Find nearest neighbor/////////////////////////
int nearestNeighbor(const Matrix<X_DIM>& pt)
{
	int closest = -1;
	double mindist = INFTY;
	double d;

	int nn = (int)rrttree.size();
	for(int i = 0; i < nn; i++){
		RRTNode & tn = rrttree[i];
		d = dist(tn.x, pt);
		if(d < mindist && (int)tn.children.size() < maxchildren && !tn.isgoal){
			closest = i;
			mindist = d;
		}
	}
	return closest;
}
////////////////////////////////////////////////////////////////////

///////////////////RRT initionlisation/////////////////////////
void initTree(const Matrix<X_DIM>& xinit)
{
	rrttree.clear();
	RRTNode n;
	n.x = xinit;
	rrttree.push_back(n);
}
///////////////////////////////////////////////////////////////

///////////////////////RRTstep/////////////////////////////////
bool rrtStep(bool& found)
{
	int size = (int) rrttree.size();
	int n_node = -1;
	int tries = 0;
	int maxtries = 500;
	Matrix<X_DIM> pt;
	//Random sample a node
	do{
		tries++;
		//goal bias
		if(Random() < 0.1){
			
			pt[0] = goal[0] + 0.1*(2.0*Random() - 1.0);
			pt[1] = goal[1] + 0.1*(2.0*Random() - 1.0);
			pt[2] = x_min[2] + Random() * (x_max[2] - x_min[2]);
			pt[3] = x_min[3] + Random() * (x_max[3] - x_min[3]);
		}
		else{
			pt[0] = x_min[0] + Random() * (x_max[0] - x_min[0]);
			pt[1] = x_min[1] + Random() * (x_max[1] - x_min[1]);
			pt[2] = x_min[2] + Random() * (x_max[2] - x_min[2]);
			pt[3] = x_min[3] + Random() * (x_max[3] - x_min[3]);
		}

		int col1 = -1;
		CAL_CheckPointCollision(cal_environment, (float)pt[0], (float)pt[1], 0, false, &col1);
		int col2 = -1;
		CAL_CheckPointCollision(cal_environment, (float)(pt[0]+dt*pt[3]*cos(pt[2])), (float)(pt[1]+dt*pt[3]*sin(pt[2])), 0, false, &col2);
		if(col1 == 0 && col2 == 0)
			n_node = nearestNeighbor(pt);
		else
			continue;

	}while(n_node == -1 && tries < maxtries);
	
	if(tries == maxtries){
		return false;
	}
	
	//found a node in rrttree: n_node


	//Random select  control for node p and find the next node x_new.
	Matrix<X_DIM> x_new;
	Matrix<X_DIM> x_old;
	Matrix<U_DIM> rancontrol;

	rancontrol[0] = u_min[0] + (u_max[0] - u_min[0])*Random(); //acceleration.
	rancontrol[1] = u_min[1] + (u_max[1] - u_min[1])*Random(); //wheel angel.
	
	
	bool valid = true;

	Dynamics dyn(dt/5, car_l);

	x_old = rrttree[n_node].x;

	for(int i = 1; i <= 5; i++){

		x_new = dyn.f(x_old, rancontrol, zeros<U_DIM,1>());

		int col1 = -1;
		CAL_CheckPointCollision(cal_environment, (float) x_new[0], (float) x_new[1], 0, false, &col1);
		int col2 = -1;
		CAL_CheckPointCollision(cal_environment, (float)x_new[0]+dt/5*x_new[3]*cos(x_new[2]), (float)x_new[1]+dt/5*x_new[3]*sin(x_new[2]), 0, false, &col2);
		int colf = -1;
		CAL_CheckLineCollision(cal_environment,(float)x_old[0],(float)x_old[1],0,(float)x_new[0],(float)x_new[1],0,false,&colf);
		if(col1 != 0 || col2 != 0 || colf != 0){
			if(col2 != 1 && n_node == 0){
				rrttree[n_node].attempts++;
				if(rrttree[n_node].attempts == 100){
					std::cout<<"The initial positon is in virtual obstacles"<<std::endl;
					return false;
				}
				std::cout<<"aaa"<<std::endl;
			}
			
			valid = false;
			break;
		}

		//check if the node satifiys bounds.
		for(int i = 0; i < X_DIM; i++){
			if(x_new[i] < x_min[i] || x_new[i] > x_max[i]){
				std::cout<<"aaa"<<std::endl;
				valid = false;
				break;
			}
		}
		if(valid == false)
			break;

		x_old = x_new;
	}

	if(valid){
		
		RRTNode newnode;
		newnode.x = x_new;
		//newnode.u = rancontrol[best_con];  
		newnode.u = rancontrol;  
		newnode.bp = n_node;

		if((newnode.x[0] - goal[0])*(newnode.x[0] - goal[0]) + (newnode.x[1] - goal[1])*(newnode.x[1] - goal[1]) < goalRadius*goalRadius /2 )
			newnode.isgoal = true;

		rrttree[n_node].children.push_back((int)rrttree.size());
		rrttree.push_back(newnode);

		if(newnode.isgoal == true){
			paths.push_back((int)rrttree.size() - 1);
			found = true;
		}
	
		
		/*Matrix<X_DIM> curr = rrttree[n_node].x;
		Matrix<U_DIM> u = rancontrol;
		for(int k = 1; k <= 5; k++){
			int np[1] = {2};
			Matrix<X_DIM> tmpx;
			tmpx = dyn.f(curr, u, zeros<U_DIM,1>());
			float p[6] = {curr[0], curr[1], 0.0, tmpx[0], tmpx[1], 0.0};
			CAL_CreatePolyline(cal_rrt, 1, np, p);
			curr = tmpx;
		}*/
	}

	return true;
}
///////////////////////////////////////////////////////////////

//////////////////////Build RRT by iterative calling rrtstep//////////////////////////
bool rrtBuild(const Matrix<X_DIM>& xinit)
{
	initTree(xinit);
	paths.clear();
	//build tree
	clock_t startClk = clock();
	bool stepRet = true;
	bool found = false;
	
	for(int i = 0; stepRet && !found; i++){
		stepRet = rrtStep(found);
		
	}

	int num_paths = (int)paths.size();
	if(num_paths == 0)
		return false;
	return true;
}
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////get K paths////////////////////////////////////////////////////////
bool get_K_Paths(int K, const Matrix<X_DIM, 1>& Start)
{
	clock_t startClk = clock();
	int colf = 0;
	CAL_CheckPointCollision(cal_environment, (float)Start[0], (float)Start[1], 0, false, &colf);
	//the starting point for RRT is in collision, exit
	if(colf != 0){
		return false;
		
	}
	int c = 0;

	//30 seconds at most for replanning.
	while(((double)(clock() - startClk) / CLOCKS_PER_SEC) < 30){
		int k = 0;
		for(k = 0; k < K; k++){
			std::cout<<k<<std::endl;
			//std::cout<<k<<std::endl;
			rrttree.clear();
		
			bool s = rrtBuild(Start);
			if(s == false){
				k--;
				c++;
				if(c == 20)
					return false;
				continue;
			}
			else if(s == true){
				std::vector<RRTNode> onePath;
				onePath.clear();
				int tmp = paths[0];
				onePath.push_back(rrttree[tmp]);

				int parent = -1;
				while(parent != 0){
					parent = rrttree[tmp].bp;
					onePath.push_back(rrttree[parent]);
					tmp = parent;
				}
				std::reverse(onePath.begin(), onePath.end());
	
				//calculate the control inputs on this path.
				for(int i = 0; i < (int)onePath.size() - 1; i++){
					onePath[i].u = onePath[i+1].u;
				}
				onePath[onePath.size() - 1].u[0] = 0;
				onePath[onePath.size() - 1].u[1] = 0;
		
				//put the current path into the path set.
				PathSet.push_back(onePath);
			}
		}
		if(k>=K)
			break;
		
	}
	if(PathSet.size() == 0){
		std::cout<<"No path found by RRT"<<std::endl;
		return false;
	}
	
	rrttree.clear();	
	return true;
}
////////////////////////Output Paths////////////////////////////////
void showPath(int& num, const std::vector<RRTNode>& path)
{
		
		CAL_CreateGroup(&num, 0, false, "Paths");
		CAL_SetGroupColor(num, 0, 1, 0);
		
		Dynamics dy(dt/5, car_l);
		for(int j = 0; j < (int)path.size()-1; j++){
			Matrix<X_DIM> curr = path[j].x;
			for(int k = 1; k <= 5; k++){
				int np[1] = {2};
				Matrix<X_DIM> tmpx;
				tmpx = dy.f(curr, path[j].u, zeros<U_DIM,1>());
				float p[6] = {curr[0], curr[1], 0.0, tmpx[0], tmpx[1], 0.0};
				CAL_CreatePolyline(num, 1, np, p);
				curr = tmpx;
			}
		}
		
}

void show_line(int& num, const Matrix<X_DIM>& x1, const Matrix<X_DIM>& x2)
{
	CAL_CreateGroup(&num, 0, false, "Paths");
	CAL_SetGroupColor(num, 0, 1, 1);
	int np[1] = {2};
	float p[6] = {x1[0], x1[1], 0.0, x2[0], x2[1], 0.0};
	CAL_CreatePolyline(num, 1, np, p);
}
////////////////////////////////////////////////////////////////////

//////////////////////Replanning/////////////////////////////////
void replanning(){
	
	//num_show and num_est are used for visualization.
	int num_show = -1;
	int num_est = -1;

	Matrix<X_DIM> x_true;
	Matrix<X_DIM> x_est = start;
	Matrix<X_DIM, X_DIM> cov = P0;
	bool success = true;
	bool replan = true;
	
	//simulte the true position of the robot by sampling from a gaussian distribution.
	x_true = sampleGaussian(x_est, cov);
	show_line(num_est, x_est, x_true);

	replan = get_K_Paths(100, x_est); //do motion planning from the current estimate position of the robot.
	if(replan == false){
		std::cout<<"Replan fail, cound not find a path"<<std::endl;
	}

	//choose the best path in the path set.
	CAL_SuspendVisualisation();

	int best = -1;
	double quality = -1;
	for(int i = 0; i < PathSet.size(); i++){
		LQGMP lqgmp(PathSet[i], dt, car_l, goalRadius);
		lqgmp.initial(x_est,cov,goal);
		double q = lqgmp.computeQuality(cal_obstacles,cal_point, cal_environment);
		if(q > quality){
			best = i;
			quality = q;
		}
	}

	CAL_ResumeVisualisation();

	//executepaht is the path that gonna be executed.
	std::vector<RRTNode> executepath;
	executepath = PathSet[best];
	showPath(num_show, executepath); //show the path before taking the first step

	//compute one step.
	Controller control(executepath, dt, car_l, goalRadius);
	success = control.simulateOneStep(x_true, x_est, cov, cal_rrt, cal_environment, goal);
	
	if(success == false){
		std::cout<<"Fail"<<std::endl;
		return;
	}
	
	while(((x_true[0]-goal[0])*(x_true[0]-goal[0])+(x_true[1]-goal[1])*(x_true[1]-goal[1])) > goalRadius* goalRadius){
			
			//remove the first node of the path that was executed one step before.
			std::vector<RRTNode> currentpath = executepath;
			std::reverse(currentpath.begin(), currentpath.end());
			currentpath.pop_back();
			std::reverse(currentpath.begin(), currentpath.end());
			
			PathSet.clear();

			std::cout<<x_est[0]<<" "<<x_est[1]<<" "<<x_est[2]<<" "<<x_est[3]<<std::endl;
			
			//replanning multiple paths, stored in pathset.
			replan = get_K_Paths(100, x_est); //do replanning from the current estimate position of the robot.
			if(replan == false){
				std::cout<<"replan fail, could not find a path"<<std::endl;
				return;
			}
			//add the hidden path to the pathSet
			if(currentpath.size() > 1){
				LQGMP lqgmp(currentpath, dt, car_l, goalRadius);
				lqgmp.initial(x_est, cov, goal);
				lqgmp.addpath(PathSet, cal_rrt, cal_environment);
			}
			std::cout<<"Planning finished"<<std::endl;

			CAL_SuspendVisualisation();
			//choose the best path in the path set.
			int best = -1;
			double quality = 0;
			for(int i = 0; i < PathSet.size(); i++){
				LQGMP lqgmp(PathSet[i], dt, car_l, goalRadius);
				lqgmp.initial(x_est,cov,goal);
				double q = lqgmp.computeQuality(cal_obstacles,cal_point, cal_environment);
				if(q > quality){
					best = i;
					quality = q;
				}
			}
			CAL_ResumeVisualisation();

			executepath = PathSet[best];			
			
			CAL_DestroyGroup(num_est);
			show_line(num_est, x_est, x_true);

			CAL_DestroyGroup(num_show);
			showPath(num_show, executepath);
			Controller control(executepath, dt, car_l, goalRadius);
			success = control.simulateOneStep(x_true, x_est, cov, cal_rrt, cal_environment, goal);
			if(success == false){
				std::cout<<"Fail"<<std::endl;
				return;
			}

	}
}
////////////////////End of replanning////////////////////////////////
int main()
{

	//srand((unsigned int) time(NULL));
	srand(412);
	// Create Obstacles in Callisto
	CAL_Initialisation (true, true, true);
	initEnvironment();
	int a = 0;
	get_K_Paths(1, start);
	showPath(a, PathSet[0]);
	
	Matrix<X_DIM> new_start = start;
	new_start[0] += 0.1;
	new_start[1] += -0.9;
	LQGMP lqgmp(PathSet[0],dt, car_l, goalRadius);
	lqgmp.initial(new_start, P0, goal);
	lqgmp.addpath(PathSet, cal_rrt, cal_environment);
	lqgmp.draw_prior_distribution(cal_ellipse);
//	LQGMP lqgmp1(PathSet[1],dt, car_l, goalRadius);
//	lqgmp1.initial(new_start, P0, goal);
//	lqgmp1.draw_prior_distribution(cal_ellipse);
//	showPath(a, PathSet[1]);
	

	for(int i = 0; i < 500; i++){
		Controller control(PathSet[0], dt, car_l, goalRadius);
		control.simulateLQG(new_start,P0, cal_rrt, cal_environment, goal);
	}

//	lqgmp.draw_prior_distribution(cal_ellipse);
	
	//get_K_Paths(200, start);
	//replanning();
//	Parallel_for_each(parallel_rrt.begin(), parallel_rrt.end(), [&](RRT& i) { 
//		i.get_K_Paths(100);
//	});
	
	int num;
	std::cin>>num;

	CAL_End();
	return 0;
}