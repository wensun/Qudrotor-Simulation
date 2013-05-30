#ifndef _DYNAMICS_
#define _DYNAMICS_


#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>

#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "RRTNode.h"


class Dynamics
{
public:
	double dt;
	double car_l;

	Dynamics(const double& d_t, const double& carl){
		dt = d_t;
		car_l = carl;
	}

	Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& m){

		Matrix<X_DIM> x_new;

		// RK4 integration
		Matrix<X_DIM> x1, x2, x3, x4, xtmp;

		xtmp = x;
		//k1
		x1[0] = xtmp[3]*cos(xtmp[2]); 
		x1[1] = xtmp[3]*sin(xtmp[2]); 
		x1[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
		x1[3] = (u[0] + m[0]);
		
		xtmp = x + 0.5*dt*x1;

		//k2
		x2[0] = xtmp[3]*cos(xtmp[2]); 
		x2[1] = xtmp[3]*sin(xtmp[2]); 
		x2[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
		x2[3] = (u[0] + m[0]);

		xtmp = x + 0.5*dt*x2;

		//k3
		x3[0] = xtmp[3]*cos(xtmp[2]); 
		x3[1] = xtmp[3]*sin(xtmp[2]); 
		x3[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
		x3[3] = (u[0] + m[0]);

		xtmp = x + dt*x3;
		//k4
		x4[0] = xtmp[3]*cos(xtmp[2]); 
		x4[1] = xtmp[3]*sin(xtmp[2]); 
		x4[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
		x4[3] = (u[0] + m[0]);

		x_new = x + (dt/6.0)*(x1 + 2.0*(x2 + x3) + x4);

		if(x_new[2] < 0)
			x_new[2] = x_new[2] + 2 * M_PI;

		if(x_new[2] > 2 * M_PI)
			x_new[2] = x_new[2] - 2 * M_PI;

		if(x_new[3] > 1)
			x_new[3] = 1;
		if(x_new[3] < 0)
			x_new[3] = 0;

		return x_new;
	}

	

	Matrix<Z_DIM,1> h(const Matrix<X_DIM>& x, const Matrix<Z_DIM>&n)
	{
		Matrix<Z_DIM,1> obs = zeros<Z_DIM,1>();
		Matrix<Z_DIM, X_DIM> H = zeros<Z_DIM, X_DIM>();
		H(0,0) = 1; H(0,1) = 0; H(0,2) = 0; H(0,3) = 0;
		H(1,0) = 0; H(1,1) = 1; H(1,2) = 0; H(1,3) = 0;
		H(2,0) = 0; H(2,1) = 0; H(2,2) = 0; H(2,3) = 1;

		obs = H*x + n ;
		return obs;
	}

	Matrix<X_DIM> f_d(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<U_DIM>& m){
		
		double deltat = dt/5;
		Matrix<X_DIM> x_start = x;
		Matrix<X_DIM> x_new;
		for(int i = 1; i <= 5; i++){
		
			Matrix<X_DIM> x1, x2, x3, x4, xtmp;

			xtmp = x_start;
			//k1
			x1[0] = xtmp[3]*cos(xtmp[2]); 
			x1[1] = xtmp[3]*sin(xtmp[2]); 
			x1[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
			x1[3] = (u[0] + m[0]);
		
			xtmp = x_start + 0.5*deltat*x1;

			//k2
			x2[0] = xtmp[3]*cos(xtmp[2]); 
			x2[1] = xtmp[3]*sin(xtmp[2]); 
			x2[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
			x2[3] = (u[0] + m[0]);

			xtmp = x_start + 0.5*deltat*x2;

			//k3
			x3[0] = xtmp[3]*cos(xtmp[2]); 
			x3[1] = xtmp[3]*sin(xtmp[2]); 
			x3[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
			x3[3] = (u[0] + m[0]);
				
			xtmp = x_start + deltat*x3;
			//k4
			x4[0] = xtmp[3]*cos(xtmp[2]); 
			x4[1] = xtmp[3]*sin(xtmp[2]); 
			x4[2] = xtmp[3]*tan(u[1] + m[1])/car_l; 
			x4[3] = (u[0] + m[0]);

			x_new = x_start + (deltat/6.0)*(x1 + 2.0*(x2 + x3) + x4);

			if(x_new[2] < 0)
				x_new[2] = x_new[2] + 2 * M_PI;

			if(x_new[2] > 2 * M_PI)
				x_new[2] = x_new[2] - 2 * M_PI;

			if(x_new[3] > 1)
				x_new[3] = 1;
			if(x_new[3] < 0)
				x_new[3] = 0;

			x_start = x_new;
		}
		return x_new;
	}

};

#endif