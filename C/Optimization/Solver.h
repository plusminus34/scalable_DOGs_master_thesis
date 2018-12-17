#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Objective.h"
#include "Constraints.h"

class Solver {
  
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) = 0;
};

class ConstrainedSolver {
  
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) = 0;
};

double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, double cur_energy = -1);