#pragma once

#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/PardisoSolver.h"

// This just used the linearized constraints but the lagrangian but solves [H,J^t;J,0] 
//  meaning H is just the obj Hessian without the second order parts of the constraints
class FeasibleIneqInteriorPoint : public IneqConstrainedSolver {
  
public:
	FeasibleIneqInteriorPoint(const double& infeasability_epsilon, const double& infeasability_filter, const int& max_newton_iters, const double& merit_p) :
			 infeasability_epsilon(infeasability_epsilon), infeasability_filter(infeasability_filter), 
			 max_newton_iters(max_newton_iters) ,merit_p(merit_p), m_solver(ai,aj,K) {m_solver.set_type(-2);}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints,
			Constraints& ineq_constraints, Eigen::VectorXd& x) = 0;

private:
	double single_homotopy_iter(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
		double mu, Eigen::VectorXd& x, double current_merit);
	double get_feasible_point(const Eigen::VectorXd& x0, double current_merit_p, Objective& obj, Constraints& eq_constraints,
            Constraints& ineq_constraints, Eigen::VectorXd& x);
	void compute_step_direction(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
		double mu, Eigen::VectorXd& d, double current_merit);

	void build_kkt_system(const Eigen::SparseMatrix<double>& hessian, const Eigen::SparseMatrix<double>& Jacobian,
						Eigen::SparseMatrix<double>& KKT);

	void build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, 
								   const std::vector<Eigen::Triplet<double> >& const_lambda_hessian, int var_n, 
								   const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n,
								   const std::vector<Eigen::Triplet<double> >& ineq_lambda_hessian,
								   const std::vector<Eigen::Triplet<double> >& ineq_jacobian_ijv, int ineq_const_n);

	double kkt_mu_error(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
		double mu);

	double get_max_alpha(const Eigen::VectorXd& x, const Eigen::VectorXd& d);
	void update_variables(Eigen::VectorXd& x,const Eigen::VectorXd& d, Constraints& ineq_constraints, double alpha);

	double ineq_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, double mu, Objective& f, 
		Constraints& eq_constraints, Constraints& ineq_constraints, double current_merit);
	double merit_func(Eigen::VectorXd& x, double mu, Objective& f, Constraints& eq_constraints, Constraints& ineq_constraints,
		double current_merit);

	const double& infeasability_epsilon; 
	const double& infeasability_filter;
	const int& max_newton_iters;
	const double& merit_p;

	double tol = 1e-6;

	std::vector<Eigen::Triplet<double> > kkt_IJV; Eigen::VectorXi cached_ijv_data;
	Eigen::SparseMatrix<double> A;

	PardisoSolver m_solver;

	bool is_feasible = false;
	bool first_solve = true;
	// lambda are the equalities multipliers, 's' are the slack variables and 'z' are their corresponding multipliers
	Eigen::VectorXd lambda,s,z;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
};