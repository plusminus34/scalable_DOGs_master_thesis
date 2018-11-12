#pragma once

#include "../Solver.h"

class LBFGSWithPenalty : public ConstrainedSolver{
  
public:
	LBFGSWithPenalty(int max_iter): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	int max_iter;
	const double OBJ_P_INIT = 0.1;
};