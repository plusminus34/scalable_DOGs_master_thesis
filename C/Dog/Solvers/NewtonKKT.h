#pragma once

#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/PardisoSolver.h"

// This just used the linearized constraints but the lagrangian but solves [H,J^t;J,0] 
//  meaning H is just the obj Hessian without the second order parts of the constraints
class NewtonKKT : public ConstrainedSolver {
  
public:
	NewtonKKT(const double& merit_p) : merit_p(merit_p), m_solver(ai,aj,K) {m_solver.set_type(-2);}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x);

private:
	void build_kkt_system(const Eigen::SparseMatrix<double>& hessian, const Eigen::SparseMatrix<double>& Jacobian,
						Eigen::SparseMatrix<double>& KKT);

	void build_kkt_system_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, 
                                     const std::vector<Eigen::Triplet<double> >& jacobian_IJV,
                                     std::vector<Eigen::Triplet<double> >& kkt_IJV);

	const double& merit_p;

	PardisoSolver m_solver;

	bool first_solve = true;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;

	// cached stuff
	Eigen::SparseMatrix<double> eps_id;
	Eigen::SparseMatrix<double> id_KKT;
};