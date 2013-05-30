#ifndef _DYNAMICS_
#define _DYNAMICS_

#define _CRT_RAND_S

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

class Dynamics{

public:
	double gravity;
	Matrix<3,1> north;

	double dt;
	//quadrotor constants
	double mass;
	Matrix<3,3> inertia;
	double forceConst;
	double momentConst;
	double dragConst;
	double roll_latency;
	double pitch_latency;
	double yaw_latency;
	double thrust_latency;
	double length;
	double minForce;
	double maxForce;
	double ballRadius;

	// camera constants
	double tanOpeningAngleX;
	double tanOpeningAngleY;
	Matrix<3,3> CamRot1;
	Matrix<3,3> CamRot2;
	Matrix<3,1> CamPos1;
	Matrix<3,1> CamPos2;

	//derivative constants;
	Matrix<3,3> invInertia;
	double nominalInput;

	Dynamics(double tau){
		dt = tau;

		gravity =	9.80665;
		north[0] =	0.70710678118654752440084436210485; north[1] = 0; north[2] = -0.70710678118654752440084436210485;
		// Quadrotor parameters
		mass        = 0.500;             // mass, kg  (source: paper)
		inertia     = 0.1*identity<3>(); // moment of inertia matrix 
		forceConst  = 1; //6.11e-8;      // N/(r/min^2)
		momentConst = 1.5e-9 / 6.11e-8;  // Nm/(r/min^2)
		dragConst   = 0.25;

		roll_latency = 20.0;
		pitch_latency = 20.0;
		yaw_latency = 20.0;
		thrust_latency = 20.0;
		length      = 0.3429/2;          // distance between center and rotor, m
		minForce = 0.09*4;

		maxForce = 3.71*4;
		ballRadius  = 0.034528125;       // radius of yellow ball, m

		// Derivative constants
		invInertia   = !inertia;
		nominalInput = gravity*mass/forceConst;

		// Camera parameters
		tanOpeningAngleX = 2.0/3.0;
		tanOpeningAngleY = 0.5;
		CamPos1[0] = 0; CamPos1[1] = -1.25; CamPos1[2] = 3;
		CamPos2[0] = 0; CamPos2[1] = 1.25;  CamPos2[2] = 3;

		CamRot1(0,0) = 1; CamRot1(1,0) = 0;  CamRot1(2,0) = 0;
		CamRot1(1,0) = 0; CamRot1(1,1) = -1; CamRot1(2,1) = 0;
		CamRot1(2,0) = 0; CamRot1(2,1) = 0;  CamRot1(2,2) = -1;
	
		CamRot2 = CamRot1;
	}

	//observation model.
	Matrix<X_DIM> h(const Matrix<X_DIM>& x, const Matrix<3,3>& R){
		Matrix<X_DIM> z;
		// accelerometer
		z[0] = -(dragConst/mass)*((R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*x[3] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*x[4] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*x[5]);
		z[1] = -(dragConst/mass)*((R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*x[3] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*x[4] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*x[5]);
		z[2] = (forceConst/mass)*x[12];

		// magnetic compass
		z[3] = (R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*north[0] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*north[1] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*north[2];
		z[4] = (R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*north[0] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*north[1] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*north[2];
		z[5] = (R(0,2) + R(0,0)*x[7] - R(0,1)*x[6])*north[0] + (R(1,2) + R(1,0)*x[7] - R(1,1)*x[6])*north[1] + (R(2,2) + R(2,0)*x[7] - R(2,1)*x[6])*north[2];

		// rate-gyros
		z[6] = x[9];
		z[7] = x[10];
		z[8] = x[11];

		// height
		z[9] = x[2];

		// pixel (in meters -- as if focal length is 1m)
		z[10] = (CamRot1(0,0)*(x[0] - CamPos1[0]) + CamRot1(1,0)*(x[1] - CamPos1[1]) + CamRot1(2,0)*(x[2] - CamPos1[2])) / (CamRot1(0,2)*(x[0] - CamPos1[0]) + CamRot1(1,2)*(x[1] - CamPos1[1]) + CamRot1(2,2)*(x[2] - CamPos1[2]));
		z[11] = (CamRot1(0,1)*(x[0] - CamPos1[0]) + CamRot1(1,1)*(x[1] - CamPos1[1]) + CamRot1(2,1)*(x[2] - CamPos1[2])) / (CamRot1(0,2)*(x[0] - CamPos1[0]) + CamRot1(1,2)*(x[1] - CamPos1[1]) + CamRot1(2,2)*(x[2] - CamPos1[2]));

		// radius (in radians)
		z[12] = asin(ballRadius / sqrt((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2])));

		return z;
	}

	//H: jacobian of the observation model.
	Matrix<Z_DIM, X_DIM> Jacobian_hx(const Matrix<X_DIM>& x, const Matrix<3,3>& R){
		Matrix<Z_DIM, X_DIM> result = zeros<Z_DIM, X_DIM>();

		result(0,3) = -(dragConst/mass)*(R(0,0) + R(0,1)*x[8] - R(0,2)*x[7]);
		result(0,4) = -(dragConst/mass)*(R(1,0) + R(1,1)*x[8] - R(1,2)*x[7]);
		result(0,5) = -(dragConst/mass)*(R(2,0) + R(2,1)*x[8] - R(2,2)*x[7]);
		result(0,7) =  (dragConst/mass)*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5]);
		result(0,8) = -(dragConst/mass)*(R(0,1)*x[3] + R(1,1)*x[4] + R(2,1)*x[5]);

		result(1,3) = -(dragConst/mass)*(R(0,1) - R(0,0)*x[8] + R(0,2)*x[6]);
		result(1,4) = -(dragConst/mass)*(R(1,1) - R(1,0)*x[8] + R(1,2)*x[6]);
		result(1,5) = -(dragConst/mass)*(R(2,1) - R(2,0)*x[8] + R(2,2)*x[6]);
		result(1,6) = -(dragConst/mass)*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5]);
		result(1,8) =  (dragConst/mass)*(R(0,0)*x[3] + R(1,0)*x[4] + R(2,0)*x[5]);

		result(2,12) = forceConst/mass;

		result(3,7) = -R(0,2)*north[0] - R(1,2)*north[1] - R(2,2)*north[2];
		result(3,8) =  R(0,1)*north[0] + R(1,1)*north[1] + R(2,1)*north[2];
	
		result(4,6) =  R(0,2)*north[0] + R(1,2)*north[1] + R(2,2)*north[2];
		result(4,8) = -R(0,0)*north[0] - R(1,0)*north[1] - R(2,0)*north[2];
			
		result(5,6) = -R(0,1)*north[0] - R(1,1)*north[1] - R(2,1)*north[2];
		result(5,7) =  R(0,0)*north[0] + R(1,0)*north[1] + R(2,0)*north[2];

		result(6,9) = 1;
		result(7,10) = 1;
		result(8,11) = 1;

		result(9,2) = 1;

		result(10,0) = ( CamPos1[1]*CamRot1(0,2)*CamRot1(1,0) - CamPos1[1]*CamRot1(0,0)*CamRot1(1,2) + CamPos1[2]*CamRot1(0,2)*CamRot1(2,0) - CamPos1[2]*CamRot1(0,0)*CamRot1(2,2) - CamRot1(0,2)*CamRot1(1,0)*x[1] + CamRot1(0,0)*CamRot1(1,2)*x[1] - CamRot1(0,2)*CamRot1(2,0)*x[2] + CamRot1(0,0)*CamRot1(2,2)*x[2])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));
		result(10,1) = (-CamPos1[0]*CamRot1(0,2)*CamRot1(1,0) + CamPos1[0]*CamRot1(0,0)*CamRot1(1,2) + CamPos1[2]*CamRot1(1,2)*CamRot1(2,0) - CamPos1[2]*CamRot1(1,0)*CamRot1(2,2) + CamRot1(0,2)*CamRot1(1,0)*x[0] - CamRot1(0,0)*CamRot1(1,2)*x[0] - CamRot1(1,2)*CamRot1(2,0)*x[2] + CamRot1(1,0)*CamRot1(2,2)*x[2])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));
		result(10,2) = (-CamPos1[0]*CamRot1(0,2)*CamRot1(2,0) - CamPos1[1]*CamRot1(1,2)*CamRot1(2,0) + CamPos1[0]*CamRot1(0,0)*CamRot1(2,2) + CamPos1[1]*CamRot1(1,0)*CamRot1(2,2) + CamRot1(0,2)*CamRot1(2,0)*x[0] - CamRot1(0,0)*CamRot1(2,2)*x[0] + CamRot1(1,2)*CamRot1(2,0)*x[1] - CamRot1(1,0)*CamRot1(2,2)*x[1])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));

		result(11,0) = ( CamPos1[1]*CamRot1(0,2)*CamRot1(1,1) - CamPos1[1]*CamRot1(0,1)*CamRot1(1,2) + CamPos1[2]*CamRot1(0,2)*CamRot1(2,1) - CamPos1[2]*CamRot1(0,1)*CamRot1(2,2) - CamRot1(0,2)*CamRot1(1,1)*x[1] + CamRot1(0,1)*CamRot1(1,2)*x[1] - CamRot1(0,2)*CamRot1(2,1)*x[2] + CamRot1(0,1)*CamRot1(2,2)*x[2])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));
		result(11,1) = (-CamPos1[0]*CamRot1(0,2)*CamRot1(1,1) + CamPos1[0]*CamRot1(0,1)*CamRot1(1,2) + CamPos1[2]*CamRot1(1,2)*CamRot1(2,1) - CamPos1[2]*CamRot1(1,1)*CamRot1(2,2) + CamRot1(0,2)*CamRot1(1,1)*x[0] - CamRot1(0,1)*CamRot1(1,2)*x[0] - CamRot1(1,2)*CamRot1(2,1)*x[2] + CamRot1(1,1)*CamRot1(2,2)*x[2])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));
		result(11,2) = (-CamPos1[0]*CamRot1(0,2)*CamRot1(2,1) - CamPos1[1]*CamRot1(1,2)*CamRot1(2,1) + CamPos1[0]*CamRot1(0,1)*CamRot1(2,2) + CamPos1[1]*CamRot1(1,1)*CamRot1(2,2) + CamRot1(0,2)*CamRot1(2,1)*x[0] - CamRot1(0,1)*CamRot1(2,2)*x[0] + CamRot1(1,2)*CamRot1(2,1)*x[1] - CamRot1(1,1)*CamRot1(2,2)*x[1])
			/ ((CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2])*(CamPos1[0]*CamRot1(0,2) + CamPos1[1]*CamRot1(1,2) + CamPos1[2]*CamRot1(2,2) - CamRot1(0,2)*x[0] - CamRot1(1,2)*x[1] - CamRot1(2,2)*x[2]));

		result(12,0) = -ballRadius * (x[0] - CamPos1[0]) / (((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2])) * sqrt((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2]) - ballRadius*ballRadius));
		result(12,1) = -ballRadius * (x[1] - CamPos1[1]) / (((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2])) * sqrt((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2]) - ballRadius*ballRadius));
		result(12,2) = -ballRadius * (x[2] - CamPos1[2]) / (((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2])) * sqrt((x[0] - CamPos1[0])*(x[0] - CamPos1[0]) + (x[1] - CamPos1[1])*(x[1] - CamPos1[1]) + (x[2] - CamPos1[2])*(x[2] - CamPos1[2]) - ballRadius*ballRadius));

		return result;
	}

	//////////////////Dynamics function: f (only one step). the derivative of the dynamics function.
	 Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<3,3>& R, const Matrix<U_DIM>& u) {
		Matrix<X_DIM> xdot;

		xdot[0] = x[3];
		xdot[1] = x[4];
		xdot[2] = x[5];

		xdot[3] =            (R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])* -dragConst/mass*((R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*x[3] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*x[4] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*x[5]) +
			(R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])* -dragConst/mass*((R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*x[3] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*x[4] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*x[5]) +
			(R(0,2) + R(0,0)*x[7] - R(0,1)*x[6])* forceConst/mass*x[12]; 

		xdot[4] =            (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])* -dragConst/mass*((R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*x[3] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*x[4] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*x[5]) +
			(R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])* -dragConst/mass*((R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*x[3] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*x[4] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*x[5]) +
			(R(1,2) + R(1,0)*x[7] - R(1,1)*x[6])* forceConst/mass*x[12]; 

		xdot[5] = -gravity + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])* -dragConst/mass*((R(0,0) + R(0,1)*x[8] - R(0,2)*x[7])*x[3] + (R(1,0) + R(1,1)*x[8] - R(1,2)*x[7])*x[4] + (R(2,0) + R(2,1)*x[8] - R(2,2)*x[7])*x[5]) +
			(R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])* -dragConst/mass*((R(0,1) - R(0,0)*x[8] + R(0,2)*x[6])*x[3] + (R(1,1) - R(1,0)*x[8] + R(1,2)*x[6])*x[4] + (R(2,1) - R(2,0)*x[8] + R(2,2)*x[6])*x[5]) +
			(R(2,2) + R(2,0)*x[7] - R(2,1)*x[6])* forceConst/mass*x[12]; 

		xdot[6] = x[9]  - x[8]*x[10] + x[7]*x[11];
		xdot[7] = x[10] + x[8]*x[9]  - x[6]*x[11];
		xdot[8] = x[11] - x[7]*x[9]  + x[6]*x[10];

		xdot[9] =  roll_latency   * (u[0] - x[9] ); 
		xdot[10] = pitch_latency  * (u[1] - x[10]); 
		xdot[11] = yaw_latency    * (u[2] - x[11]); 

		xdot[12] = thrust_latency * (u[3] - x[12]); 

		return xdot;
	}

	//jacobian of the dynamic function
	Matrix<X_DIM, X_DIM> Jacobian_fx(const Matrix<X_DIM>& x, const Matrix<3,3>& R, const Matrix<U_DIM>& u) {
		Matrix<X_DIM, X_DIM> result = zeros<X_DIM, X_DIM>();

		result(0,3) = 1;
		result(1,4) = 1;
		result(2,5) = 1;

		result(3,3) = -dragConst/mass*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8])*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) 
			-dragConst/mass*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8])*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]);
		result(3,4) = -dragConst/mass*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8])*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) 
			-dragConst/mass*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8])*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]);
		result(3,5) = -dragConst/mass*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			-dragConst/mass*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]);
		result(3,6) = -forceConst/mass*R(0,1)*x[12] - dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) 
			- dragConst/mass*R(0,2)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]));
		result(3,7) =  forceConst/mass*R(0,0)*x[12] + dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) 
			+ dragConst/mass*R(0,2)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(3,8) =  dragConst/mass*(R(0,0)*x[3] + R(1,0)*x[4] + R(2,0)*x[5])*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) 
			-dragConst/mass*(R(0,1)*x[3] + R(1,1)*x[4] + R(2,1)*x[5])*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) 
			+dragConst/mass*R(0,0)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8])) 
			-dragConst/mass*R(0,1)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(3,12) = forceConst/mass*(R(0,2) + R(0,0)*x[7] - R(0,1)*x[6]);

		result(4,3) = -dragConst/mass*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8])*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) 
			-dragConst/mass*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8])*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]);
		result(4,4) = -dragConst/mass*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8])*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) 
			-dragConst/mass*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8])*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]);
		result(4,5) = -dragConst/mass*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8])
			-dragConst/mass*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]);

		result(4,6) = -forceConst/mass*R(1,1)*x[12] - dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) 
			- dragConst/mass*R(1,2)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]));
		result(4,7) =  forceConst/mass*R(1,0)*x[12] + dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) 
			+ dragConst/mass*R(1,2)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(4,8) =  dragConst/mass*(R(0,0)*x[3] + R(1,0)*x[4] + R(2,0)*x[5])*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) 
			-dragConst/mass*(R(0,1)*x[3] + R(1,1)*x[4] + R(2,1)*x[5])*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) 
			+dragConst/mass*R(1,0)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8])) 
			-dragConst/mass*R(1,1)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(4,12) = forceConst/mass*(R(1,2) + R(1,0)*x[7] - R(1,1)*x[6]);

		result(5,3) = -dragConst/mass*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			-dragConst/mass*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]);
		result(5,4) = -dragConst/mass*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			-dragConst/mass*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]);
		result(5,5) = -dragConst/mass*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			-dragConst/mass*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]);
		result(5,6) = -forceConst/mass*R(2,1)*x[12] - dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			- dragConst/mass*R(2,2)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]));
		result(5,7) =  forceConst/mass*R(2,0)*x[12] + dragConst/mass*(R(0,2)*x[3] + R(1,2)*x[4] + R(2,2)*x[5])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]) 
			+ dragConst/mass*R(2,2)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(5,8) =  dragConst/mass*(R(0,0)*x[3] + R(1,0)*x[4] + R(2,0)*x[5])*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8]) 
			-dragConst/mass*(R(0,1)*x[3] + R(1,1)*x[4] + R(2,1)*x[5])*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]) 
			+dragConst/mass*R(2,0)*(x[3]*(R(0,1) + R(0,2)*x[6] - R(0,0)*x[8]) + x[4]*(R(1,1) + R(1,2)*x[6] - R(1,0)*x[8]) + x[5]*(R(2,1) + R(2,2)*x[6] - R(2,0)*x[8])) 
			-dragConst/mass*R(2,1)*(x[3]*(R(0,0) - R(0,2)*x[7] + R(0,1)*x[8]) + x[4]*(R(1,0) - R(1,2)*x[7] + R(1,1)*x[8]) + x[5]*(R(2,0) - R(2,2)*x[7] + R(2,1)*x[8]));
		result(5,12) = forceConst/mass*(R(2,2) + R(2,0)*x[7] - R(2,1)*x[6]);

		result(6,7)  = x[11];
		result(6,8)  = -x[10];
		result(6,9)  = 1;
		result(6,10) = -x[8];
		result(6,11) = x[7];

		result(7,6)  = -x[11];
		result(7,8)  = x[9];
		result(7,9)  = x[8];
		result(7,10) = 1;
		result(7,11) = -x[6];

		result(8,6)  = x[10];
		result(8,7)  = -x[9];
		result(8,9)  = -x[7];
		result(8,10) = x[6];
		result(8,11) = 1;

		result(9,9) = -roll_latency;
		result(10,10) = -pitch_latency;
		result(11,11) = -yaw_latency;
		result(12,12) = -thrust_latency;
		return result;
	}

	//jacobian of the dynamics function related to input.
	Matrix<X_DIM, U_DIM> Jacobian_fu(const Matrix<X_DIM>& x, const Matrix<3,3>& R, const Matrix<U_DIM>& u) {
		Matrix<X_DIM, U_DIM> result = zeros<X_DIM, U_DIM>();

		result(9,0)  = roll_latency;
		result(10,1) = pitch_latency;
		result(11,2) = yaw_latency;
		result(12,3) = thrust_latency;

		return result;
	}

	//compute the linearized version jacobian
	void linearizeDiscretize(const Matrix<X_DIM>& x, const Matrix<3,3>& R, const Matrix<U_DIM>& u, const Matrix<X_DIM,X_DIM>& M, Matrix<X_DIM,X_DIM>& A, Matrix<X_DIM,U_DIM>& B, Matrix<X_DIM,X_DIM>& MM) {
		Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
		Matrix<X_DIM,U_DIM> G = Jacobian_fu(x,R,u);

		A = exp(dt*F);

		Matrix<X_DIM,X_DIM> A2 = exp((dt*0.5)*F);

		// Simpson's integration approximation of Int[0,dt] exp(F*T) dT G
		B = (dt/6) * (G + 4*(A2*G) + A*G);
		MM = (dt/6) * (M + 4*A2*M*~A2 + A*M*~A);
	}

	//without noise: dynamics function.
	void propagate(Matrix<X_DIM>& x, Matrix<3,3>& R, const Matrix<U_DIM>& u) {
		Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
		Matrix<X_DIM> xDot = f(x,R,u);

		Matrix<X_DIM,X_DIM> A = exp(dt*F);

		Matrix<X_DIM,X_DIM> A2 = exp((dt*0.5)*F);
		x = x + (dt/6)*(xDot + 4*(A2*xDot) + A*xDot);
		// reset rotation-error in xNew into RNew
		R = R*exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
		x[6] = 0; x[7] = 0; x[8] = 0;
	}

	void propagate_noise(Matrix<X_DIM>& x, Matrix<3,3>& R, const Matrix<U_DIM>& u, const Matrix<X_DIM,X_DIM>& M) {
		Matrix<X_DIM,X_DIM> F = Jacobian_fx(x,R,u);
		Matrix<X_DIM> xDot = f(x,R,u);
		Matrix<X_DIM,X_DIM> A = exp(dt*F);
		Matrix<X_DIM,X_DIM> A2 = exp((dt*0.5)*F);
		Matrix<X_DIM,X_DIM> MM = (dt/6) * (M + 4*A2*M*~A2 + A*M*~A);
		x = x + (dt/6)*(xDot + 4*(A2*xDot) + A*xDot) + sampleGaussian(zeros<X_DIM,1>(), MM);
		// reset rotation-error in xNew into RNew
		R = R*exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
		x[6] = 0; x[7] = 0; x[8] = 0;
	}

};

#endif