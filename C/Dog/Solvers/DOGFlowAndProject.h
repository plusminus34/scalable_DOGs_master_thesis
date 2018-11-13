#pragma once

#include "../Dog.h"
#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/PardisoSolver.h"

class DOGFlowAndProject : public ConstrainedSolver{
  
public:
	DOGFlowAndProject(const Dog& dog,  double flow_t, int max_iter);
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	double solve_single_iter(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	double flow(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);
	void project(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& obj, double cur_energy = -1);

	int max_iter;
	double flow_t;
	// initial dog
	const Dog& dog_init;
	bool first_solve;

	// Linear system solving related
	PardisoSolver m_solver;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
};