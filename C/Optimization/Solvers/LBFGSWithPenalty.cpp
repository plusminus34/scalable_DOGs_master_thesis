#include "LBFGSWithPenalty.h"

#include "../CompositeObjective.h"
#include "../QuadraticConstraintsSumObjective.h"

#include "LBFGS.h"

double LBFGSWithPenalty::solve_constrained(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints,
													 Eigen::VectorXd& x) {
	auto obj_val = f.obj(x0);
	std::cout << "penalty_repetitions = " << penalty_repetitions << std::endl;
	if (!penalty_repetitions) return obj_val;

	x = x0;
	// Solve a few times with different penalty methods. Every time set an objective with varying penalty on the constraints
	for (int i = 0; i < penalty_repetitions; i++) {
		obj_val = solve_single_iter_with_fixed_p(x0,f,constraints,x);
		std::cout << "solved repetition = " << i << std::endl;
		p*=0.5;
	}
	return obj_val;
}

double LBFGSWithPenalty::solve_single_iter_with_fixed_p(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x) {
	QuadraticConstraintsSumObjective constObj(constraints);
	CompositeObjective compObj({&constObj, &f}, {1,p});
	LBFGS solver(max_lbfgs_iter);
	auto obj_val = solver.solve(x, compObj, x);
	return obj_val;
}