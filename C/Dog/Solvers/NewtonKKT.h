#pragma once

#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/PardisoSolver.h"

// This just used the linearized constraints but the lagrangian but solves [H,J^t;J,0] 
//  meaning H is just the obj Hessian without the second order parts of the constraints
class NewtonKKT : public ConstrainedSolver {
  
public:
	NewtonKKT(const double& merit_p) : merit_p(merit_p), m_solver(ai,aj,K) {m_solver.set_type(-2);}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x);

	const double& merit_p;

	PardisoSolver m_solver;

	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
};