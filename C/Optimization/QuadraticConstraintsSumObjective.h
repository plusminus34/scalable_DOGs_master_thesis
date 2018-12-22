#pragma once

#include <igl/Timer.h>

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
		return cnst->Vals(x).squaredNorm();
	}
	
	Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// The derivative of the sum of squares is twice times the sume of the gradients times the values
		return 2*cnst->Jacobian(x).transpose()*cnst->Vals(x);
	}

//std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M)
	// TODO: This completely ignores the hessian of the constraints! (works perfectly for linear constraints such as positions though)
	virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x) {
		igl::Timer timer; auto init_time = timer.getElapsedTime();
		// The constraints are linear so the hessian is just simple 2*J'*J
		// (Here it's even simpler as just 2*(diagonals with one at constrained indices) but this is simpler to write)
		auto J = cnst->Jacobian(x);
		cachedH = 2*J.transpose()*J;
		double H_time = timer.getElapsedTime()-init_time;
		return cachedH;
	}

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x) {
		// Could write it directly maybe, or have an eddificent A*A' at least..
		// Or preallocate the IJV (second time)
		auto J = cnst->Jacobian(x);
		cachedH = 2*J.transpose()*J;
		IJV = to_triplets(cachedH);
	}

	Constraints* cnst;
};