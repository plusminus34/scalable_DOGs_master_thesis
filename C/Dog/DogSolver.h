#pragma once

#include "Dog.h"

#include "igl/serialize.h"

#include "../Optimization/CompositeConstraints.h"
#include "../Optimization/CompositeObjective.h"
#include "../Optimization/EdgesAngleConstraints.h"
#include "../Optimization/EdgePointConstraints.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/PointPairConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"

#include "../Optimization/Solvers/Newton.h"
#include "../Optimization/Solvers/NewtonKKT.h"
#include "../Optimization/Solvers/EqSQP.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/StitchingConstraints.h"
#include "Objectives/MVTangentCreaseAngleConstraints.h"
#include "Objectives/PairedBoundaryVerticesBendingObjective.h"

#include "../Folding/FoldingBinormalBiasConstraints.h"
#include "../Folding/FoldingMVBiasConstraints.h"

#include "Objectives/LinearConstraints.h"
#include "Objectives/ProximalObjective.h"
#include "Objectives/SomeSerialObjective.h"
#include "Objectives/SubEdgesAngleConstraints.h"

using std::vector;
using std::pair;

class DogSolver {
public:

	enum SolverMode {mode_standard, mode_subsolvers, mode_vsadmm, mode_jadmm, mode_proxjadmm, mode_serial, mode_procrustes, mode_experimental};

	struct Params : public igl::Serializable {
		double bending_weight = 1.;
		double paired_boundary_bending_weight = 1.;
		double isometry_weight = 20000;
		double stitching_weight = 10000;
		double soft_pos_weight = 5;
		double dihedral_weight = 1000;
		double pair_weight = 0;
		double fold_bias_weight = 1;
		double mv_bias_weight = 0;
		double merit_p = 10;
		int max_newton_iters = 5;
		double infeasability_epsilon = 1e-3;
		double infeasability_filter = 1e-1;
		double convergence_threshold = 1e-6;
		bool folding_mode = true;
		bool flip_sign = false;
		double admm_rho = 1;
		double admm_gamma = 1;

		void InitSerialization() {
			Add(bending_weight,std::string("bending_weight"));
			Add(paired_boundary_bending_weight,std::string("paired_boundary_bending_weight"));
			Add(isometry_weight,std::string("isometry_weight"));
			Add(stitching_weight,std::string("stitching_weight"));
			Add(soft_pos_weight,std::string("soft_pos_weight"));
			Add(dihedral_weight,std::string("dihedral_weight"));
			Add(pair_weight,std::string("pair_weight"));
			Add(fold_bias_weight,std::string("fold_bias_weight"));
			Add(mv_bias_weight,std::string("mv_bias_weight"));
			Add(merit_p,std::string("merit_p"));
			Add(max_newton_iters,std::string("max_newton_iters"));
			Add(infeasability_epsilon,std::string("infeasability_epsilon"));
			Add(infeasability_filter,std::string("infeasability_filter"));
			Add(convergence_threshold,std::string("convergence_threshold"));
			Add(folding_mode,std::string("folding_mode"));
			Add(flip_sign,std::string("flip_sign"));
		}
	};

	DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
		std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
		std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
		std::vector<std::pair<int,int>>& pairs,
		std::vector<std::pair<int,int>>& bnd_vertices_pairs,
		std::ofstream* time_measurements_log = NULL);

	~DogSolver();

	//DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, DogSolver::Params& p);

	// returns [bending, isometry, sum]
	Eigen::VectorXd get_obj_parts();

	void set_opt_vars(const Eigen::VectorXd& x_i) { x = x_i;}
	Eigen::VectorXd get_opt_vars() { return x;}

	void set_solver_mode(SolverMode mode_new);
	bool is_subsolver(){return (mode != mode_standard && sub_dogsolver.size()<1);}

	void single_iteration(double& constraints_deviation, double& objective);
	void single_iteration_fold(double& constraints_deviation, double& objective);
	void single_iteration_subsolvers(double& constraints_deviation, double& objective);
	void single_iteration_ADMM(double& constraints_deviation, double& objective);
	void single_iteration_normal(double& constraints_deviation, double& objective);
	void single_iteration_serial(double& constraints_deviation, double& objective);
	void single_iteration_procrustes(double& constraints_deviation, double& objective);
	void single_iteration_experimental(double& constraints_deviation, double& objective);

	void update_edge_coords(Eigen::MatrixXd& edgeCoords) {constraints.edgePtConst.update_coords(edgeCoords);}
	void update_point_coords(Eigen::VectorXd& bc);
	void update_edge_angles(const std::vector<double> cos_angles_i);
	void update_mv_cos_angles(const std::vector<double> cos_angles_i) {constraints.mvTangentCreaseAngleConst.set_angles(cos_angles_i);}
	void update_obj_weights(const std::vector<double>& weights_i);

	//for subsolvers
	void update_w_coords(const Eigen::MatrixXd& W);

	Dog& getDog(){return dog;}

	bool is_folded();
	bool is_mountain_valley_correct(const Eigen::VectorXd& x);

	void get_x_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T);
	void get_y_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T);
	void set_x_rotation(Eigen::Matrix3d& R);
	void set_y_rotation(Eigen::Matrix3d& R);

	//VSADMM
	void build_VSADMMObjective(const Eigen::SparseMatrix<double>& A);
	void build_ProximalObjective(const Eigen::SparseMatrix<double>& P);
	void remake_compobj();
	void set_z(const Eigen::VectorXd vec){ admm_z = vec; }
	Eigen::VectorXd get_z(){return admm_z;}
	void set_lambda(const Eigen::VectorXd vec){ admm_lambda = vec; }
	Eigen::VectorXd get_Ax(){return admm_A*x;}
	Eigen::VectorXd get_lambda(){return admm_lambda;}


	struct Constraints {
		Constraints(const Dog& dog, const Eigen::VectorXd& init_x0, Eigen::VectorXi& b, Eigen::VectorXd& bc,
			std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
			std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
			std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
			std::vector<std::pair<int,int>>& pairs);

		DogConstraints dogConst;
		StitchingConstraints stitchingConstraints;
		PositionalConstraints posConst;
		EdgePointConstraints edgePtConst;
		EdgesAngleConstraints edgeAngleConst;
		SubEdgesAngleConstraints subEdgesAngleConst;
		MVTangentCreaseAngleConstraints mvTangentCreaseAngleConst;
		PointPairConstraints ptPairConst;
		CompositeConstraints compConst;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
	  			 Constraints& constraints,
	  			FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
	  			FoldingMVBiasConstraints& foldingMVBiasConstraints,
	  			std::vector<std::pair<int,int>>& bnd_vertices_pairs,
	  			const DogSolver::Params& p);

	  	SimplifiedBendingObjective bending;
	  	IsometryObjective isoObj;
      	QuadraticConstraintsSumObjective pointsPosSoftConstraints;
      	QuadraticConstraintsSumObjective edgePosSoftConstraints;
      	QuadraticConstraintsSumObjective edgeAnglesSoftConstraints;
      	QuadraticConstraintsSumObjective mvTangentCreaseSoftConstraints;
      	QuadraticConstraintsSumObjective ptPairSoftConst;
      	QuadraticConstraintsSumObjective foldingBinormalBiasObj;
      	QuadraticConstraintsSumObjective foldingMVBiasObj;
      	QuadraticConstraintsSumObjective stitchingConstraintsPenalty;
      	PairedBoundaryVerticesBendingObjective pairedBndVertBendingObj;
				//Plus QuadraticConstraintsSumObjective from LinearConstraints
				//Plus proximal objective
				//plus something for the serial mode
      	CompositeObjective compObj;
	};

private:

	Dog& dog;
	Eigen::VectorXd x; // variables
	Eigen::VectorXd x0; // variables at the beginning of time step
	bool is_constrained;
	FoldingBinormalBiasConstraints foldingBinormalBiasConstraints;
	FoldingMVBiasConstraints foldingMVBiasConstraints;

	// The constraints needs to be defined before the objectives, as some of the objective are dependent on constraints
	DogSolver::Constraints constraints;
	DogSolver::Objectives obj;
	DogSolver::Params& p;

	SolverMode mode = mode_experimental;
	int iter_i=0;

	//for submeshes
	vector< Dog* > sub_dog;
	vector< DogSolver* > sub_dogsolver;

	vector< Eigen::VectorXi > sub_b;
	vector< Eigen::VectorXd > sub_bc;
	vector< vector<int> > sub_ij_to_bc;

	vector< vector< EdgePoint > > constrained_edge_points;//uses local indices
	vector< Eigen::MatrixXd > sub_edgeCoords;
	vector< vector< pair<int,int> > > corresponding_edge_points;//stored as <submesh, idx>
	void update_sub_edgeCoords();

	vector< vector<int> > sub_idx_to_angle_idx;
	vector< Eigen::MatrixXi > sub_edge_angle_ww;

	Eigen::VectorXi empty_xi;
	Eigen::VectorXd empty_xd;
	std::vector<EdgePoint> empty_ep;
	Eigen::MatrixXd empty_mat;
	std::vector<std::pair<Edge,Edge>> empty_egg;
	std::vector<MVTangentCreaseFold> empty_thing;
	std::vector<double> empty_d;
	std::vector<pair<int,int>> empty_pair;

	//ADMM
	LinearConstraints vsadmmConst;
	QuadraticConstraintsSumObjective *vsadmm_obj = nullptr;
	ProximalObjective *pjadmm_obj = nullptr;
	SomeSerialObjective *serial_obj = nullptr;
	Eigen::VectorXd admm_lambda;
	Eigen::VectorXd admm_z;
	Eigen::SparseMatrix<double> admm_A, admm_P;

	//Procrustes
	std::vector< Eigen::SparseMatrix<double> > proc_T;

	// Solvers
	//NewtonKKT newtonKKT;
	EqSQP newtonKKT;
	//FeasibleIneqInteriorPoint interiorPt;

	std::ofstream* time_measurements_log;
};
