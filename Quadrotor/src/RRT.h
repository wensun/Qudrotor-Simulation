#ifndef _RRT_
#define _RRT_

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
#include "Dynamics.h"


class RRT{

public:

	struct TreeNode{
		Matrix<X_DIM> T; //x, y, theta, v
		Matrix<U_DIM> u; //v, phi
		Matrix<3,3> Rotation; //rotation matrix, computed from the 6,7,8 entries of the vector T.


		int bp;
		bool marked;
		int attempts;

		double depth;
		double clearance;

		std::vector<int> children;

		TreeNode()
		{
			T.reset();
			u.reset();
			bp = -1;
			marked = false;
			depth = 0.0;
			clearance = 9999999;
			attempts = 0;
		}
	};

	struct PathNode{
		Matrix<X_DIM> T;
		Matrix<U_DIM> u;
		Matrix<3,3> Rotation;
	};

	Matrix<X_DIM> start;
	Matrix<3,3> Rotinitial;
	Matrix<2> goal;
	double dt;
	double rthreshold;

	double factor;
	double planbias;
	double maxstep;

	int maxchildren;
	int maxattempts;
	int maxiter;
	int maxNNiter;
	int maxtries;
	double plantime;
	double K;
	std::vector<PathNode> rrtpath;
	std::vector<TreeNode> rrttree;
	std::vector<int> rrtpaths;
	std::vector<std::vector<PathNode>> pathSet;

	int cal_obstacles;

	double AWMAX;
	double AWMIN;
	double ForceMAX;
	double ForceMIN;

	RRT(const Matrix<X_DIM, 1>& s, const Matrix<3,3> RotInitial, const Matrix<2>& gl, const double& timestep, const double& ptime, const double& gr, const int& cal_obs){
		start = s;
		goal = gl;
		Rotinitial = RotInitial;
		dt = timestep;
		plantime = ptime;
		rthreshold = gr;
		pathSet.clear();
		
		Dynamics dy(dt);
		AWMAX = 1/3 * 3.134159; //60 degree/m
		AWMIN = -AWMIN;			
		ForceMAX = dy.nominalInput*2;
		ForceMIN = dy.nominalInput*0.5;
		cal_obstacles = cal_obs;
	}

	void setPlannerDefaults(){
		factor = 2.0;
		K = 0.35 * 3.1415926;
		//w_max = k / ( v + 1).

		planbias = 0.1;
		maxiter = 1000;
		maxNNiter = 1000;
		maxtries = 1000;
		maxstep = 1;
		maxchildren = 15;
		maxattempts = 10;
	}

	void initTree(std::vector<TreeNode>& tree, const Matrix<X_DIM>& pose, const Matrix<3,3>& Rotation){
		TreeNode n;
		n.T = pose;
		n.Rotation = Rotation;
		n.bp = -1;
		n.depth = 0.0;
		n.marked = false;
		n.attempts = 0;
		tree.push_back(n);
	}
	double dist(const Matrix<X_DIM>& p1, const Matrix<X_DIM>& p2){
		double dis = 0;
		Matrix<3> pp1 = p1.subMatrix<3,1>(0,0);
		Matrix<3> pp2 = p2.subMatrix<3,1>(0,0);
		dis = sqrt(tr(~(pp1 - pp2) * (pp1 - pp2)));
		return dis;
	}

	int nearestNeighbor(const std::vector<TreeNode>& tree, const Matrix<X_DIM>& point){
		int closest = -1;
		double mindist = 10000;
		for(int i = 0; i < tree.size(); i++){
			Matrix<X_DIM> tmpnode = tree[i].T;
			double tmpdis = dist(tmpnode, point);
			if(tmpdis < mindist && tree[i].attempts < maxattempts && tree[i].marked == false){
				closest = i;
				mindist = tmpdis;
			}
		}
		return closest;
	}

	bool rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths, bool& found,  const clock_t& sc){
		found = false;
		if((double)((clock() - sc)/CLOCKS_PER_SEC) > plantime)
			return false;

		int node = -1;
		int tries = 0;

		Matrix<X_DIM> point;
		point.reset();
		double rgoal = rthreshold;

		tries = 0;
		do{
			if(tries > 25)
				rgoal = 0.5;

			if(random() < planbias){
				point[0] = goal[0] + 0.1*(2.0*random() - 1.0);
				point[1] = goal[1] + 0.1*(2.0*random() - 1.0);
				point[2] = goal[2] + 0.1*(2.0*random() - 1.0);

			}
			else{
				point[0] = 0 + random() * 5;
				point[1] = 0 + random() * 10;
				point[2] = 0 + random() * 5;
				
			}

			int col = -1;
			CAL_CheckPointCollision(cal_obstacles, point[0], point[1], point[2], false, &col);

			if(col == 0){
				node = nearestNeighbor(tree, point);
				if(node != -1){
					tree[node].attempts ++;
				}
			}

		} while ((node == -1) && (++tries < maxNNiter));
	
		if(tries == maxtries){
			return false;  //return false when cannot sample
		}

		if(node == -1)
			return false;

		Matrix<X_DIM> x_new, x_old;
		Matrix<3,3> Rot_new; Matrix<3,3> Rot_old;
		Rot_new.reset(); Rot_old.reset();
		x_new.reset(); x_old.reset();
		Matrix<U_DIM> rancontrol;
		rancontrol.reset();
		x_old = tree[node].T;
		Rot_old = tree[node].Rotation;

		rancontrol[0] = AWMIN + random() * (AWMAX - AWMIN);
		rancontrol[1] = AWMIN + random() * (AWMAX - AWMIN);
		rancontrol[2] = AWMIN + random() * (AWMAX - AWMIN);
		rancontrol[3] = ForceMIN + random() * (ForceMAX - ForceMIN);

		bool valid = false;
		int col = -1;
		double tau = dt * 1.0 / 5.0;
		for(int seg = 0; seg < 5; seg++){
			Dynamics dyn(tau);
			dyn.propagate(x_old, Rot_old, x_new, Rot_new, rancontrol);
			CAL_CheckLineCollision(cal_obstacles, x_old[0], x_old[1], 0.0, x_new[0], x_new[1], 0.0, false, &col);
			if(col != 0)
				break;
			x_old = x_new;
			Rot_old = Rot_new;
		}
	
		if(col == 0){
			TreeNode newnode;
			newnode.T = x_new;
			newnode.Rotation = Rot_new;
			newnode.bp = node;
			newnode.u = rancontrol;
			newnode.marked = false;
			newnode.attempts = 0;

			int newid = (int) tree.size();
			tree[node].children.push_back(newid);
			tree.push_back(newnode);
			
			double dg = tr(~(x_new.subMatrix<2,1>(0,0) - goal) * (x_new.subMatrix<2,1>(0,0) - goal));
			if(dg < rthreshold * rthreshold)
			{
				TreeNode& tnode = tree[newid];
				tnode.marked = true;
				paths.push_back(newid);
				found = true;
				return true;
			}
		}
			return true;
	}

	bool buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, const clock_t& sc){
		bool stepRet = true;
		// Initialize tree
		initTree(tree, start, Rotinitial);

		bool found = false;
		for (int i = 0; stepRet ; ++i) {
			stepRet = rrtStep(tree, paths, found, sc);
			if(found == true)
				break;
		}
		if (stepRet && !paths.empty()) 
		{
			//std::cout << "\nRRT build: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
			int numpaths = (int)paths.size();
			//std::cout << "Num paths found: " << numpaths << std::endl;

		} 
		else {
			std::cout << "Unable to find solution, reusing previous solution" << std::endl;
		}

		return (stepRet && !paths.empty());
	}

	bool executeRRT(const clock_t& sc)
	{
		rrttree.clear();
		rrtpath.clear();
		rrtpaths.clear();

		return buildRRT(rrttree, rrtpaths, sc);
	}

	bool Plan_K_Seconds(){
		clock_t startCLK = clock();
	
		while(((double)(clock() - startCLK) / CLOCKS_PER_SEC) < plantime){
	
			bool s = executeRRT(startCLK);
			if(s == true){
			
				std::vector<PathNode> onePath;
				onePath.clear();
				int tmp = rrtpaths[0];
				PathNode tmpnode;
				tmpnode.T = rrttree[tmp].T;
				tmpnode.u = rrttree[tmp].u;
				tmpnode.Rotation = rrttree[tmp].Rotation;
				onePath.push_back(tmpnode);

				int parent = -1;
				while(parent != 0){
					parent = rrttree[tmp].bp;
					tmpnode.T = rrttree[parent].T;
					tmpnode.u = rrttree[parent].u;
					tmpnode.Rotation = rrttree[parent].Rotation;
					onePath.push_back(tmpnode);
					tmp = parent;
				}
				std::reverse(onePath.begin(), onePath.end());

				for(int i = 0; i < (int) onePath.size() - 1; i++){
					onePath[i].u = onePath[i+1].u;
				}
				onePath[(int)onePath.size() - 1].u[0] = 0;
				onePath[(int)onePath.size() - 1].u[1] = 0;

				pathSet.push_back(onePath);
				onePath.clear();
			}	
		}

		if(pathSet.size() == 0){
			std::cout<<"No Path find in "<<plantime<<" seconds"<<std::endl;
			return false;
		}
		//std::cout<<pathSet.size()<<std::endl;
		return true;
	}

	void showPath(const int& cal_plan){
		
		for(int i = 0; i < (int)pathSet.size(); i++){
			for(int t = 0; t < (int)pathSet[i].size() - 1; t++){
				double tau = dt*1.0/5;
				Dynamics dyn(tau);
				Matrix<X_DIM> oldtmp = pathSet[i][t].T;
				Matrix<3,3> oldRot = pathSet[i][t].Rotation;
				Matrix<X_DIM> newtmp = oldtmp;
				Matrix<3,3> newRot = oldRot;
				Matrix<U_DIM> u = pathSet[i][t].u;
				for(int k = 0; k < 5; k++){
					dyn.propagate(oldtmp, oldRot, newtmp, newRot, u);
					float line[6] = {oldtmp[0], oldtmp[1], 0, newtmp[0], newtmp[1], 0};
					int np[1] = {2};
					CAL_CreatePolyline(cal_plan, 1, np, line);
					oldtmp = newtmp;
					oldRot = newRot;
				}
			}
		}
	}
};


#endif