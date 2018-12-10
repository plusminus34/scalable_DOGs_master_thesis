#pragma once

#include "../Solver.h"
#include "lbfgs/LBFGSSolver.h"

class LBFGS : public Solver {
  
public:
	LBFGS(int max_iter) {
        param.epsilon = 1e-5;
        param.delta = 1e-5;
        param.max_iterations = max_iter;

        solver = new LBFGSpp::LBFGSSolver<double>(param);
    }
    ~LBFGS() {delete solver;}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x);

	void set_max_iter(int max_iter) {
        param.max_iterations = max_iter;
        delete solver;
        solver = new LBFGSpp::LBFGSSolver<double>(param);
    }
private:
    LBFGSpp::LBFGSParam<double> param;
    // Create solver and function object
    LBFGSpp::LBFGSSolver<double>* solver;
};

class LBFGS_obj_grad_interface {
public:
    LBFGS_obj_grad_interface(const Objective& obj) : objective(obj) {}
    double obj(const Eigen::VectorXd& x) {return objective.obj(x);}
    void grad(const Eigen::VectorXd& x, Eigen::VectorXd& grad) { grad = objective.grad(x);}

    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& g) {
        grad(x,g);
        return obj(x);
     }
private:
	const Objective& objective;
};