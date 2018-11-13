#include "LBFGSWithPenalty.h"

#include "../CompositeObjective.h"
#include "../QuadraticConstraintsSumObjective.h"

#include "LBFGS.h"

double LBFGSWithPenalty::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints,
													 Eigen::VectorXd& x) {
	
	x = x0;
	
	// Solve a few times with different penalty methods. Every time set an objective with varying penalty on the constraints
	double f;
	for (int i = 0; i < 10; i++) {
		f = solve_single_iter_with_fixed_p(x0,obj,constraints,x);
		p*=0.5;
	}
	return f;
}

double LBFGSWithPenalty::solve_single_iter_with_fixed_p(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
	QuadraticConstraintsSumObjective constObj(constraints);
	CompositeObjective compObj({&constObj, &obj}, {1,p});
	LBFGS solver(max_iter);
	auto f = solver.solve(x, compObj, x);
	return f;
}