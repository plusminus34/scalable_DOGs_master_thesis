#pragma once

#include "../Solver.h"

class LBFGS {
  
public:
	LBFGS(int max_iter): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, const Objective& obj, Eigen::VectorXd& x);

	void set_max_iter(int max_iter) {max_iter = max_iter;}
private:
	int max_iter;
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
