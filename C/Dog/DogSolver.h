#pragma once

#include "Dog.h"

#include "../Optimization/CompositeConstraints.h"
#include "../Optimization/CompositeObjective.h"
#include "../Optimization/EdgePointConstraints.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"

#include "../Optimization/Solvers/Newton.h"
#include "../Optimization/Solvers/NewtonKKT.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGGuess.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/HEnergy.h"
#include "Objectives/LaplacianSimilarity.h"
#include "Objectives/StitchingConstraints.h"


class DogSolver {
public:
	enum SolverType {
		SOLVE_NONE = 0,
		SOLVE_NEWTON_PENALTY = 1,
		SOLVE_NEWTON_FLOW = 2
	};
	struct Params {
		DogSolver::SolverType solverType = SOLVE_NEWTON_FLOW;
		double bending_weight = 1.;
		double isometry_weight = 0.1;
		double laplacian_similarity_weight = 0;
		double soft_pos_weight = 1;
		int penalty_repetitions = 1;
		double merit_p = 1;
		int max_newton_iters = 10;
		double infeasability_epsilon = 1e-3;
		bool align_procrustes = true;
	};

	DogSolver(Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0, const DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords);
	
	void single_iteration(double& constraints_deviation, double& objective);
	void update_edge_coords(Eigen::MatrixXd& edgeCoords) {constraints.edgePtConst.update_coords(edgeCoords);}
	void update_point_coords(Eigen::VectorXd& bc) {constraints.posConst.update_coords(bc);}

	Dog& getDog(){return dog;}
	const QuadTopology& getQuadTop(){return quadTop;}
	
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
	  			PositionalConstraints& posConst,
	  			EdgePointConstraints& edgePtConst,
	  			const DogSolver::Params& p);

	  	SimplifiedBendingObjective bending;
	  	IsometryObjective isoObj;
	  	//LaplacianSimilarity laplacianSimilarity;
      	QuadraticConstraintsSumObjective pointsPosSoftConstraints;
      	QuadraticConstraintsSumObjective edgePosSoftConstraints;
      	LaplacianSimilarity laplacianSimilarity;
      	CompositeObjective compObj;
	};

private:
	Dog& dog;
	const QuadTopology& quadTop;
	bool is_constrained;

	// Optimization parameters
	Eigen::VectorXd init_x0;
	// The constraints needs to be defined before the objectives, as some of hte objective are dependent on constraints
	DogSolver::Constraints constraints;
	DogSolver::Objectives obj;
	const DogSolver::Params& p;

	// Solvers
	DOGGuess dogGuess;
	Newton newton;
	NewtonKKT newtonKKT;
};