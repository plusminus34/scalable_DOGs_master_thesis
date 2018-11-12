#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"
#include "Objective.h"

class QuadraticConstraintsSum : public Objective {
  
public:
	QuadraticConstraintsSum(const Constraints& constraints): cnst(constraints) {};
	double obj(const Eigen::VectorXd& x) const {
		return cnst.Vals(x).squaredNorm();
	}
	
	Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// The derivative of the sum of squares is twice times the sume of the gradients times the values
		return 2*cnst.Jacobian(x).transpose()*cnst.Vals(x);
	}

private:
	const Constraints& cnst;
};