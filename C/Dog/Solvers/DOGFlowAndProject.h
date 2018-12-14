#pragma once

#include "../Dog.h"
#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/LBFGSWithPenalty.h"

class DOGFlowAndProject : public ConstrainedSolver{
  
public:
	DOGFlowAndProject(const Dog& dog,  double flow_t, const int& max_flow_project_iter, const int& max_lbfgs_proj_iter = 400, 
								const int& penalty_repetitions = 1);
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	double solve_single_iter(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x, bool project_after_flow = true);

	void set_max_iter(int max_iter) {max_iter = max_iter;}

	void resetSmoother() {lbfgsWithPenalty.resetSmoother();}
private:
	double flow(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);
	void project(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& obj, double cur_energy = -1);

	const int& max_flow_project_iter; 
	const int& max_lbfgs_proj_iter;
	const int& penalty_repetitions;

	double flow_t;
	// initial dog
	const Dog& dog_init;
	bool first_solve;

	LBFGSWithPenalty lbfgsWithPenalty;

	/*
	// Linear system solving related
	PardisoSolver m_solver;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
	*/
};