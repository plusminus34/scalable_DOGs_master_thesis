#include "IpOptSolver.h"
/*
#include <ifopt/problem.h>

double IpOptSolver::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, 
            Eigen::VectorXd& x) {
	// 1. define the problem
    Problem nlp;
    nlp.AddVariableSet  (std::make_shared<IpOptVariables(x0)>());
    nlp.AddConstraintSet(std::make_shared<IpOptConstraints(constraints)>());
    nlp.AddCostSet      (std::make_shared<IpOptObjective(obj)>());
    nlp.PrintCurrent();

    // 2. choose solver and options
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    ipopt.SetOption("jacobian_approximation", "exact");

    // 3 . solve
    ipopt.Solve(nlp);
    x = nlp.GetOptVariables()->GetValues();
}
*/