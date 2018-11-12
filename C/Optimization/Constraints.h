#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Constraints {
public:
	int getConstNum(){return const_n;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x) const = 0;
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const = 0;

protected:
	int const_n;
};