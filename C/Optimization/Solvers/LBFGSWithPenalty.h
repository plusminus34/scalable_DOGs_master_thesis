#pragma once

#include "../Solver.h"

class LBFGSWithPenalty : public ConstrainedSolver {
  
public:
	LBFGSWithPenalty(const int& max_lbfgs_iter, const int& penalty_repetitions): max_lbfgs_iter(max_lbfgs_iter), 
								penalty_repetitions(penalty_repetitions), p (OBJ_P_INIT) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x);
	double solve_single_iter_with_fixed_p(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x);

	//void set_max_iter(int max_lbfgs_iter) {max_lbfgs_iter = max_lbfgs_iter;}
	void resetSmoother() { p = OBJ_P_INIT;}
private:
	const int& max_lbfgs_iter;
	const int& penalty_repetitions;
	const double OBJ_P_INIT = 0.1;
	double p;
};