#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"
#include "Objective.h"

// No hessian at the moment (see https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm for a matrix notation of the hessian)
class QuadraticConstraintsSumObjective : public Objective {
  
public:
	QuadraticConstraintsSumObjective(const Constraints& constraints): cnst(constraints) {};
	virtual QuadraticConstraintsSumObjective* clone() const {return new QuadraticConstraintsSumObjective(*this);}
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