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
#include "Solvers/DOGGuess.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/StitchingConstraints.h"

#include "../Folding/FoldingBinormalBiasConstraints.h"


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
		double soft_pos_weight = 1;
		double fold_bias_weight = 100;
		int penalty_repetitions = 1;
		double merit_p = 10;
		int max_newton_iters = 30;
		double infeasability_epsilon = 1e-3;
		double infeasability_filter = 1e-1;
		bool align_procrustes = false;
	};

	DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, const DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
		std::vector<std::pair<int,int>>& pairs
		/*CurvedFoldingBiasObjective& curvedFoldingBiasObj*/);
	
	void single_iteration(double& constraints_deviation, double& objective);
	void update_edge_coords(Eigen::MatrixXd& edgeCoords) {constraints.edgePtConst.update_coords(edgeCoords);}
	void update_point_coords(Eigen::VectorXd& bc) {constraints.posConst.update_coords(bc);}

	Dog& getDog(){return dog;}

	bool is_folded();
	
	struct Constraints {
		Constraints(const Dog& dog, Eigen::VectorXi& b, Eigen::VectorXd& bc,
			std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
			std::vector<std::pair<int,int>>& pairs);

		DogConstraints dogConst;
		StitchingConstraints stitchingConstraints;
		PositionalConstraints posConst;
		EdgePointConstraints edgePtConst;
		PointPairConstraints ptPairConst;
		CompositeConstraints compConst;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
	  			PositionalConstraints& posConst,
	  			EdgePointConstraints& edgePtConst,
	  			PointPairConstraints& ptPairConst,
	  			FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
	  			const DogSolver::Params& p);

	  	SimplifiedBendingObjective bending;
	  	IsometryObjective isoObj;
      	QuadraticConstraintsSumObjective pointsPosSoftConstraints;
      	QuadraticConstraintsSumObjective edgePosSoftConstraints;
      	QuadraticConstraintsSumObjective ptPairSoftConst;
      	QuadraticConstraintsSumObjective foldingBinormalBiasObj;
      	CompositeObjective compObj;
	};

private:
	Dog& dog;
	bool is_constrained;
	FoldingBinormalBiasConstraints foldingBinormalBiasConstraints;

	// Optimization parameters
	Eigen::VectorXd init_x0;
	// The constraints needs to be defined before the objectives, as some of the objective are dependent on constraints
	DogSolver::Constraints constraints;
	DogSolver::Objectives obj;
	const DogSolver::Params& p;

	// Solvers
	DOGGuess dogGuess;
	Newton newton;
	EqSQP newtonKKT;
	FeasibleIneqInteriorPoint interiorPt;
};