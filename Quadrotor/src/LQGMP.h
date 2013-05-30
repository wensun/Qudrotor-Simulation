#ifndef _LQGMP_
#define _LQGMP_

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>

#include "RRTNode.h"
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Controller.h"

class LQGMP
{
public:
	int path_size;
	double car_l;
	Matrix<X_DIM> start; //the mean of the start distribution
	Matrix<X_DIM, X_DIM> P0; //the covariance of the start distribution
	Matrix<DIM> goal;
	double goal_radius;
	double dt;
	std::vector<RRTNode> path;
	std::vector<Matrix<2*X_DIM, 2*X_DIM>> F;
	std::vector<Matrix<2*X_DIM, U_DIM+Z_DIM>> G;
	std::vector<Matrix<X_DIM+U_DIM, X_DIM+U_DIM>> Y;
	std::vector<Matrix<2*X_DIM>> y_bar; 
	std::vector<Matrix<X_DIM>> x_bar;
	std::vector<Matrix<X_DIM + U_DIM>> x_u;
	std::vector<Matrix<2*X_DIM, 2*X_DIM>> R;
	Matrix<U_DIM+Z_DIM, U_DIM+Z_DIM> Q;
	

	LQGMP(const std::vector<RRTNode>& Path, const double& deltat, const double& CARL, const double& g_r){
		path = Path;
		path_size = path.size();

		F.resize(path_size);
		G.resize(path_size);
		Y.resize(path_size);
		R.resize(path_size);
		y_bar.resize(path_size);
		x_bar.resize(path_size);
		x_u.resize(path_size);
		dt = deltat;
		car_l = CARL;
		goal_radius = g_r;

	}

	void initial(const Matrix<X_DIM>& Start, const Matrix<X_DIM, X_DIM>& P, const Matrix<DIM>& Goal){
		start = Start;
		P0 = P;
		goal = Goal;
		Controller contr(path, dt, car_l, goal_radius);
		Q.reset();
		Q.insert(0,0, contr.M);
		Q.insert(U_DIM, U_DIM, contr.N);
	}

	void compute_F_G_R_Y(){
		Matrix<X_DIM+U_DIM, 2*X_DIM> z;
		z.reset();
		for(int i = 0; i < path_size; i++){
			F[i].reset();
			G[i].reset();
			R[i].reset();
			Y[i].reset();
			y_bar[i].reset();
			x_bar[i].reset();
			x_u[i].reset();
		}

		
		Controller control(path, dt, car_l, goal_radius);
		control.compute_A_B_V();
		control.compute_K(P0);
		control.compute_L();

		y_bar[0].insert(0,0, start - path[0].x);
		y_bar[0].insert(X_DIM, 0, start - path[0].x);
		x_bar[0] = y_bar[0].subMatrix<X_DIM,1>(0,0) + path[0].x;
		R[0].insert(0,0, P0);
		F[0] = zeros<2*X_DIM, 2*X_DIM>();
		G[0] = zeros<2*X_DIM, U_DIM+Z_DIM>();

		z.insert(0,0, identity<X_DIM>());
		z.insert(X_DIM, X_DIM, control.L[0]);
		Matrix<X_DIM+U_DIM> tmp;
		tmp.insert(0,0, path[0].x);
		tmp.insert(X_DIM,0, path[0].u);
		x_u[0] = z*y_bar[0] + tmp;
		Y[0] = z*R[0]*~z;

		for(int t = 1; t < path_size; t++){
			//compute F_t
			F[t].insert(0,0, control.A[t]);
			F[t].insert(0, X_DIM, control.B[t]*control.L[t-1]);
			F[t].insert(X_DIM, 0, control.K[t]*control.H*control.A[t]);
			F[t].insert(X_DIM, X_DIM, control.A[t]+control.B[t]*control.L[t-1]-control.K[t]*control.H*control.A[t]);

			//compute G_t
			G[t].insert(0,0, control.V[t]);
			G[t].insert(0, U_DIM, zeros<X_DIM, Z_DIM>());
			G[t].insert(X_DIM, 0, control.K[t]*control.H*control.V[t]);
			G[t].insert(X_DIM, U_DIM, control.K[t]*control.W);

			y_bar[t] = F[t]*y_bar[t-1];
			x_bar[t] = y_bar[t].subMatrix<X_DIM,1>(0,0) + path[t].x;
			R[t] = F[t]*R[t-1]*~F[t] + G[t]*Q*~G[t];
			
			z.insert(0,0, identity<X_DIM>());
			z.insert(X_DIM, X_DIM, control.L[t]);
			tmp.insert(0,0, path[t].x);
			tmp.insert(X_DIM,0, path[t].u);
			x_u[t] = z*y_bar[t] + tmp;
			Y[t] = z*R[t]*~z;
		}
	}

	void draw_prior_distribution(const int& cal_ellipse){

		compute_F_G_R_Y();
		
		for(int i = 0; i < path_size; i++){
			Matrix<DIM, DIM> Evec, Eval;
			jacobi(Y[i].subMatrix<DIM, DIM>(0,0), Evec, Eval);

			float or[3];
			float sc[3];

			or[0] = 0.0;
			or[1] = 0.0;
			or[2] = (float)mod2pi(atan2(Evec(1,0), Evec(0,0)));

			sc[0] = 3.0*(float)(sqrt(Eval(0,0)));
			sc[1] = 3.0*(float)(sqrt(Eval(1,1)));
			sc[2] = 1.0;

			int id;
			CAL_CreateUserDrawn(cal_ellipse, drawUnitCircle, NULL, 0.0f, 0.0f, 0.0f, &id);
			CAL_SetObjectPosition(id, x_bar[i][0], x_bar[i][1], 0.0);
			CAL_SetObjectOrientation(id, or[0], or[1], or[2]);
			CAL_SetObjectScaling(id, sc[0], sc[1], sc[2]);
		}
	}

	//compute the number of the standard deviation 2D: x,y
	double computeObsConfidence(const Matrix<DIM>& x_pos, const Matrix<DIM,DIM>& x_cov, int & cal_obstacles,int& cal_point, int& cal_environment){
		Matrix<DIM, DIM> EVec, EVal;
		jacobi(x_cov, EVec, EVal);
		Matrix<3,3> Temp = identity<3>();
		Temp.insert(0,0, ~EVec);
		Matrix<4,1> q = quatFromRot(Temp);
		CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
		CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

		Matrix<DIM,DIM> invScale = zeros<DIM,DIM>();
		invScale(0,0) = 1/(float)sqrt(EVal(0,0));
		invScale(1,1) = 1/(float)sqrt(EVal(1,1));
		Matrix<DIM> transPos =  invScale * ~EVec * x_pos;

		CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

		int num_pairs;
		CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);
		SCALResult* results = new SCALResult[num_pairs];
		CAL_GetResults(results);
		double distance = results[0].distance;
		delete[] results;


		CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
		CAL_SetGroupScaling(cal_environment, 1,1,1);
		return distance;

	}

	double computeQuality(int & cal_obstacles,int& cal_point, int& cal_environment){
		compute_F_G_R_Y();
		double quality = 1;
		for(int t = 0; t < path.size(); t++){
			double dis = computeObsConfidence(x_bar[t].subMatrix<DIM,1>(0,0),Y[t].subMatrix<DIM,DIM>(0,0),cal_obstacles, cal_point, cal_environment);
			double prob = incompletegamma(0.5*X_DIM, 0.5*dis*dis);
			quality = quality * prob;
		}
		//std::cout<<"The probability of collision avoidance:"<<quality * 100 <<std::endl;
		return quality;
	}

	//use the lqg simulate without noise to get the hidden path, that path will be added to the path set.
	void addpath(std::vector<std::vector<RRTNode>>& PathSet,  const int& cal_rrt, const int& cal_environment)
	{

		std::vector<RRTNode> tmppath;
		tmppath.resize(path_size);
		//compute the control around this path.
		Controller control(path, dt, car_l, goal_radius);
		control.compute_A_B_V();
		control.compute_K(P0);
		control.compute_L();
		Matrix<U_DIM> u, u_d; //u_d is deviation u - u*.
		Matrix<X_DIM, X_DIM> P;
		Matrix<X_DIM> x_true, x_est, x_true_old, x_d_est;
		Matrix<Z_DIM> z_d; //observation z - h(x*,0)
		P = P0;
		//at the starting position, x_est is the mean, P0 is the covairance,x_true is the simulated true position
		x_est = start;
		x_d_est = x_est - path[0].x; //x_0 - x_0*
		x_true = start; //no noise
		tmppath[0].x = x_true; 

		for(int k = 1; k < path_size; k++){
			u_d = control.L[k-1] * x_d_est;
			u = path[k-1].u + u_d;
			//check if the control violates the constraints, 
		//	if(u[0] < -0.4)
		//		u[0] = -0.4;
		//	if(u[0] > 0.8)
		//		u[0] = 0.8;
		//	if(u[1] < -M_PI*0.35)
		//		u[1] = -M_PI*0.35;
		//	if(u[1] > M_PI*0.35)
		//		u[1] = M_PI*0.35;

			tmppath[k-1].u = u;

			//use the dynamic function f to get the next true state of the robot. without any noise.
			Dynamics dy(dt/5, car_l);		
			for(int iter = 1; iter <=5; iter++){
				x_true_old = x_true;
				x_true = dy.f(x_true_old, u, zeros<U_DIM,1>());
			}
			tmppath[k].x = x_true;

			//kalman Filter
			x_d_est = control.A[k]*x_d_est + control.B[k]*u_d;
			z_d = dy.h(x_true, zeros<Z_DIM,1>()) - dy.h(path[k].x, zeros<Z_DIM,1>());
			x_d_est = x_d_est + control.K[k]*(z_d - control.H * x_d_est);
		}
		PathSet.push_back(tmppath);
	}

	
};


#endif