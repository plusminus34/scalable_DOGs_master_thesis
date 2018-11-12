#pragma once

#include "../../Optimization/Solver.h"

class DOGFlowAndProject : public ConstrainedSolver{
  
public:
	DOGFlowAndProject(int max_iter): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	double solve_single_iter(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	void flow(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);
	void project(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	int max_iter;
};