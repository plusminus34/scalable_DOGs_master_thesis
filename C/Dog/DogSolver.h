#pragma once

#include "Dog.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGGuess.h"
#include "Solvers/Newton.h"
#include "Solvers/NewtonKKT.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/HEnergy.h"
#include "Objectives/LaplacianSimilarity.h"

#include "Objectives/StitchingConstraints.h"
#include "../Optimization/CompositeConstraints.h"


class DogSolver {
public:
	DogSolver(Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0, const DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	
	void single_iteration(double& constraints_deviation, double& objective);
	void update_edge_coords(Eigen::MatrixXd& edgeCoords& edgeCoords) {constraints.edgePtConst.update_coords(edgeCoords)}
	void update_point_coords(Eigen::VectorXd& bc) {constraints.posConst.update_coords(bc)}
	
	enum SolverType {
		SOLVE_NONE = 0,
		SOLVE_NEWTON_PENALTY = 1,
		SOLVE_NEWTON_FLOW = 2
	};

	struct Params {
		DeformationController::SolverType solverType = SOLVE_NEWTON_FLOW;
		double bending_weight = 1.;
		double isometry_weight = 0.1;
		double laplacian_similarity_weight = 0;
		double soft_pos_weight = 1;
		int penalty_repetitions = 1;
		double merit_p = 1;

		bool align_procrustes = false;
	};

	struct Constraints {
		Constraints(const Dog& dog, const QuadTopology& quadTop,
			Eigen::VectorXi& b, Eigen::VectorXd& bc,
			std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);

		DogConstraints dogConst;
		StitchingConstraints stitchingConstraints;
		PositionalConstraints posConst;
		EdgePointConstraints edgePtConst;
		CompositeConstraints compConst;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0,
	  			EdgePointConstraints& edgePtConst);

	  	SimplifiedBendingObjective bending;
	  	IsometryObjective isoObj;
	  	//LaplacianSimilarity laplacianSimilarity;
      	QuadraticConstraintsSumObjective pointsPosSoftConstraints;
      	QuadraticConstraintsSumObjective edgePosSoftConstraints;
      	CompositeObjective compObj;
	};

private:
	Dog& dog;
	const QuadTopology& quadTop;

	// Optimization parameters
	Eigen::VectorXd init_x0;
	DogSolver::Objectives obj;
	DogSolver::Constraints constraints;
	const DogSolver::Params& p;

	// Solvers
	Newton newton;
	NewtonKKT newtonKKT;
	DOGGuess dogGuess;
};