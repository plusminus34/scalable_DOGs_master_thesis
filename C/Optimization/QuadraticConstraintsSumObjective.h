#pragma once

#include <igl/Timer.h>

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"
#include "Objective.h"

// Not a complete hessian at the moment (see https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm for a matrix notation of the hessian)
// The hessian here only takes into account the first derivative from the constraint itself, and is of the form 2*Jt*J
class QuadraticConstraintsSumObjective : public Objective {
  
public:
	QuadraticConstraintsSumObjective(Constraints& constraints): cnst(&constraints) {};
	virtual QuadraticConstraintsSumObjective* clone() const {return new QuadraticConstraintsSumObjective(*this);}
	double obj(const Eigen::VectorXd& x) const {
		return cnst->Vals(x).squaredNorm();
	}
	
	Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// The derivative of the sum of squares is twice times the sume of the gradients times the values
		return 2*cnst->Jacobian(x).transpose()*cnst->Vals(x);
	}

private:
	// TODO: This completely ignores the hessian of the constraints! (works perfectly for linear constraints such as positions though)
	virtual void updateHessianIJV(const Eigen::VectorXd& x) {
		// Could write it directly maybe, or have an eddificent A*A' at least..
		// Or preallocate the IJV (second time)
		auto J = cnst->Jacobian(x);
		cachedH = 2*J.transpose()*J;
		IJV = to_triplets(cachedH);
	}

	Constraints* cnst;
};