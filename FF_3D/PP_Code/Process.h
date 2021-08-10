#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "Eigen\src\Eigenvalues\EigenSolver.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <Eigen\SparseCholesky>
#include "MeshDefinition.h"
#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include "PardisoSolver.h"
#include "Acceleration\Anderson.h"

using namespace Eigen;
using namespace std;

typedef Triplet<double> Tri;

class Process
{
public:
	Process(string filename, string initial_filename);
	~Process();

	bool writeobj(VolumeMesh& out_mesh, const string out_filename);
	bool writeobj(const string out_filename);

	void init();				

	void pck();
	void solve();

	double Energy();
	void local();
	void global();

	void Distortion();
	
protected:
	double	lb;						// lower bound on SVs(-1 = disabled)
	double	ub;						// upper bound on SVs(-1 = disabled)
	double	tol_err;				
	double	tol_energy;	
	int		iter_max;				// maximal number of BD projection iterations

	int		n_constraints;			// the number of constraints
	SparseMatrix<double> A;			// coefficient of constraints
	VectorXd			 b;			// right side value of constraints

	VectorXd Tx;					// lifted variable differentials
	VectorXd pTx;					// projected 
	VectorXd distortions;			// distortions of differentials i.e. condition number
	VectorXd flips;					
	//SimplicialCholesky<SparseMatrix<double>> TCholesky;		// precompile
	SparseMatrix<double> T;			// lifts init_position(d*n) into a higher dimension space of differentials(m*d*d)
	vector<double> T_coef;

	int			F_N;
	int			V_N;
	MatrixXi	F;
	VectorXd	source_position;
	VectorXd	global_position;
	VolumeMesh	source_mesh;
	VolumeMesh	init_mesh;

	VectorXd	cur_position;
	VectorXd	K;					// conformal distortion bound
	double		cur_energy;
	double		pre_energy;
	int			total_iter;			// number of solving equations

	OpenVolumeMesh::IO::FileManager fileManager;

	PardisoSolver* pardiso;
};

