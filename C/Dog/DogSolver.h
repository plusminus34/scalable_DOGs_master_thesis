#pragma once

#include "Dog.h"

#include "Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Solvers/DOGFlowAndProject.h"
#include "Solvers/DOGGuess.h"

#include "../Optimization/Solvers/IpOptSolver.h"


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
		SOLVE_IPOPT = 3
	};

	DogSolver() : state(NULL) {};
	void init_from_new_dog(Dog& dog, const QuadTopology& quadTop);
	
	void single_optimization();
	void update_positional_constraints();
	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};

	struct Params {
		DogSolver::DeformationType deformationType = CURVE_DEFORMATION;
		DogSolver::SolverType solverType = SOLVE_PENALTY;
		double bending_weight = 1.;
		double isometry_weight = 100.;
		int max_lbfgs_routines = 400;
		double const_obj_penalty = 100.;
		int penalty_repetitions = 1;

		double folding_angle = 0;
		double curve_timestep = 0;

		bool align_procrustes;
		bool arap_guess;
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
		DOGFlowAndProject flowProject;
		DOGGuess dogGuess;
		FoldingAnglePositionalConstraintsBuilder angleConstraintsBuilder;
		CurveInterpolationConstraintsBuilder curveConstraintsBuilder;
	};

	DogSolver::State* state;
	// Positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	// Curve constraints
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
};