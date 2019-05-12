#pragma once

#include "Dog.h"

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
#include "../Optimization/Solvers/FeasibleIneqInteriorPoint.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/StitchingConstraints.h"
#include "Objectives/MVTangentCreaseAngleConstraints.h"
//#include "Objectives/PointsRigidAlignmentObjective.h"
#include "Objectives/CurveAffineSymmetryConstraint.h"
#include "Objectives/Curve2AffineSymmetriesConstraint.h"
#include "Objectives/Curve2AffineCommuteConstraint.h"

#include "../Folding/FoldingBinormalBiasConstraints.h"
#include "../Folding/FoldingMVBiasConstraints.h"

using std::vector;


class DogSolver {
public:
	
	struct Params {
		double bending_weight = 1.;
		double isometry_weight = 20000;
		double soft_pos_weight = 1;
		double dihedral_weight = 10;
		double fold_bias_weight = 1;
		double mv_bias_weight = 0;
		double wallpaper_curve_weight = 1;
		int penalty_repetitions = 1;
		double merit_p = 10;
		int max_newton_iters = 5;
		double infeasability_epsilon = 1e-3;
		double infeasability_filter = 1e-1;
		double convergence_threshold = 1e-6;
		bool folding_mode = true;
		bool flip_sign = false;
	};

	DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
		std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
		std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
		std::vector<std::pair<int,int>>& pairs,
		std::pair<vector<int>,vector<int>>& matching_curve_pts_y,
		std::pair<vector<int>,vector<int>>& matching_curve_pts_x,
		std::ofstream* time_measurements_log = NULL);
	
	void single_iteration(double& constraints_deviation, double& objective);
	void single_iteration_fold(double& constraints_deviation, double& objective);
	void single_iteration_normal(double& constraints_deviation, double& objective);
	void update_edge_coords(Eigen::MatrixXd& edgeCoords) {constraints.edgePtConst.update_coords(edgeCoords);}
	void update_point_coords(Eigen::VectorXd& bc) {constraints.posConst.update_coords(bc);}
	void update_edge_angles(const std::vector<double> cos_angles_i) {constraints.edgeAngleConst.set_angles(cos_angles_i);}
	void update_mv_cos_angles(const std::vector<double> cos_angles_i) {constraints.mvTangentCreaseAngleConst.set_angles(cos_angles_i);}

	Dog& getDog(){return dog;}

	bool is_folded();

	void get_x_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T);
	void get_y_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T);
	
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
		MVTangentCreaseAngleConstraints mvTangentCreaseAngleConst;
		PointPairConstraints ptPairConst;
		CompositeConstraints compConst;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
	  			 Constraints& constraints,
	  	/*
	  			PositionalConstraints& posConst,
	  			EdgePointConstraints& edgePtConst,
	  			EdgesAngleConstraints& edgeAngleConst,
	  			PointPairConstraints& ptPairConst,*/
	  			FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
	  			FoldingMVBiasConstraints& foldingMVBiasConstraints,
	  			QuadraticConstraintsSumObjective& affineAlignmentSoft,
	  			QuadraticConstraintsSumObjective& affineCommuteSoft,
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
      	CompositeObjective compObj;
	};

private:
	static Eigen::VectorXd init_variables(const Eigen::VectorXd& init_mesh_vars, 
			std::pair<vector<int>,vector<int>>& matching_curve_pts_x,
			std::pair<vector<int>,vector<int>>& matching_curve_pts_y);
	Dog& dog;
	Eigen::VectorXd x; // variables
	bool is_constrained;
	FoldingBinormalBiasConstraints foldingBinormalBiasConstraints;
	FoldingMVBiasConstraints foldingMVBiasConstraints;
	//CurveAffineSymmetryConstraint affineAlignment;
	Curve2AffineSymmetriesConstraint affineAlignment;
	QuadraticConstraintsSumObjective affineAlignmentSoft;
	Curve2AffineCommuteConstraint affineCommuteConst;
	QuadraticConstraintsSumObjective affineCommuteSoft;

	// The constraints needs to be defined before the objectives, as some of the objective are dependent on constraints
	DogSolver::Constraints constraints;
	DogSolver::Objectives obj;
	DogSolver::Params& p;

	// Solvers
	//NewtonKKT newtonKKT;
	EqSQP newtonKKT;
	//FeasibleIneqInteriorPoint interiorPt;

	std::ofstream* time_measurements_log;
};