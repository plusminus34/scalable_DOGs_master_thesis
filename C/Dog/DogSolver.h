#pragma once

#include "Dog.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGFlowAndProject.h"
#include "Solvers/DOGGuess.h"
#include "Solvers/Newton.h"
#include "Solvers/NewtonKKT.h"

#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/HEnergy.h"
#include "Objectives/LaplacianSimilarity.h"

#include "../Optimization/Solvers/LBFGS.h"


class DogSolver {
public:
	enum DeformationType {
		DIHEDRAL_FOLDING = 0,
		CURVE_DEFORMATION = 1
	};
	enum SolverType {
		SOLVE_NONE = 0,
		SOLVE_FLOW_PROJECT = 1,
		SOLVE_PENALTY = 2,
		SOLVE_LBFGS = 3,
		SOLVE_NEWTON_PENALTY = 4,
		SOLVE_NEWTON_FLOW = 5
	};

	DogSolver() : state(NULL) {};
	void init_from_new_dog(Dog& dog, const QuadTopology& quadTop);
	
	void single_optimization();
	void update_positional_constraints();
	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};

	struct Params {
		DogSolver::DeformationType deformationType = CURVE_DEFORMATION;
		DogSolver::SolverType solverType = SOLVE_NEWTON_FLOW;
		double bending_weight = 1.;
		double isometry_weight = 0.1;
		double laplacian_similarity_weight = 0;
		double diag_length_weight = 0;
		int max_lbfgs_routines = 400;
		double const_obj_penalty = 1;
		int penalty_repetitions = 1;
		double merit_p = 1;
		bool project_after_flow = true;

		double folding_angle = 0;
		double curve_timestep = 0;

		bool align_procrustes = false;
		bool arap_guess;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0) : 
	  		bending(quadTop), isoObj(quadTop, init_x0), laplacianSimilarity(dog,init_x0) {}
	  // Objectives
	  SimplifiedBendingObjective bending;
	  //HEnergy bending(state->quadTop);
	  IsometryObjective isoObj;
	  LaplacianSimilarity laplacianSimilarity;
	};

	DogSolver::Params p;
	double constraints_deviation;
	double objective;
	
private:
	void init_solver_state(Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0);

	struct State {
		State(Dog& dog, const QuadTopology& quadTop, const DogSolver::Params& p, const Eigen::VectorXd& init_x0);

		Dog& dog;
		const QuadTopology& quadTop;
		const DogSolver::Params& p;
		DOGFlowAndProject flowProject;
		Newton newton;
		NewtonKKT newtonKKT;
		LBFGS lbfsgSolver;
		DOGGuess dogGuess;
		FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;
		CurveInterpolationConstraintsBuilder curveConstraintsBuilder;
		CurveInterpolationConstraintsBuilder geoConstraintsBuilder;

		Eigen::VectorXd init_x0;
		DogSolver::Objectives obj;
	};

	DogSolver::State* state;
	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	// Curve constraints
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};