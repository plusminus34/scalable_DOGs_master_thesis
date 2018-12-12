#pragma once

#include "../Solver.h"
#include "lbfgs/LBFGSSolver.h"
#include "igl/Timer.h"

class LBFGS : public Solver {
  
public:
	LBFGS(int max_iter) {
        param.epsilon = 1e-10;
        param.delta = 1e-10;
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
    double obj(const Eigen::VectorXd& x) {
        auto t = timer.getElapsedTime();
        auto obj =  objective.obj(x);
        f_time += timer.getElapsedTime()-t; f_cnt++;
        return obj;
    }
    void grad(const Eigen::VectorXd& x, Eigen::VectorXd& grad) { 
        auto t = timer.getElapsedTime();
        grad = objective.grad(x); g_cnt++;
        g_time += timer.getElapsedTime()-t;
    }

    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& g) {
        grad(x,g);
        return obj(x);
     }
     double g_time = 0; int g_cnt = 0;
     double f_time = 0; int f_cnt = 0;
private:
	const Objective& objective;
    igl::Timer timer;
};