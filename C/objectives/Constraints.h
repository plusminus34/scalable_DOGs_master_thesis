#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Quad.h"
#include "Energy.h"

class Constraints {
public:
	Constraints(const QuadTopology& quadTop) : Energy(quadTop) {}
	int getConstNum(){return const_n;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) = 0;
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x)  = 0;
	
	//virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) = 0;

	/*
	double obj(const Eigen::VectorXd& x) {
		return Vals(x).squaredNorm();
	}

	Eigen::VectorXd grad(const Eigen::VectorXd& x) {
		// The derivative of the sum of squares is twice times the sume of the gradients times the values
		return 2*Jacobian(x).transpose()*Vals(x);
	}
	*/
protected:
	int const_n;
};