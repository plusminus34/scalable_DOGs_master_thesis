#pragma once

#include "Dog.h"

#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGFlowAndProject.h"

class DogSolver {
public:
	DogSolver() : state(NULL) {};
	void init_from_new_dog(Dog& dog, const QuadTopology& quadTop);
	
	void single_optimization();
	void update_positional_constraints();
	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};

	struct Params {
		double bending_weight = 1.;
		double isometry_weight = 100.;
		bool fold_mesh = true;
		double folding_angle = 0;
		int max_lbfgs_routines = 400;
		double const_obj_penalty = 100.;
		int penalty_repetitions = 1;
	};

	DogSolver::Params p;
	
private:
	void init_solver_state(Dog& dog, const QuadTopology& quadTop);

	struct State {
		State(Dog& dog, const QuadTopology& quadTop, const DogSolver::Params& p);

		Dog& dog;
		const QuadTopology& quadTop;
		const DogSolver::Params& p;
		DOGFlowAndProject flowProject;// = NULL;
		FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;// = NULL;
	};

	DogSolver::State* state;
	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
};