#pragma once

#include "Dog.h"

#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGFlowAndProject.h"


class DogSolver {
public:
	enum DeformationType {
		DIHEDRAL_FOLDING = 0,
		CURVE_DEFORMATION = 1
	};

	DogSolver() : state(NULL) {};
	void init_from_new_dog(Dog& dog, const QuadTopology& quadTop);
	
	void single_optimization();
	void update_positional_constraints();
	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};

	struct Params {
		DogSolver::DeformationType deformationType;
		double bending_weight = 1.;
		double isometry_weight = 100.;
		int max_lbfgs_routines = 400;
		double const_obj_penalty = 100.;
		int penalty_repetitions = 1;

		double folding_angle = 0;
		double curve_timestep = 0;
	};

	DogSolver::Params p;
	double constraints_deviation;
	double objective;
	
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
	// Curve constraints
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};