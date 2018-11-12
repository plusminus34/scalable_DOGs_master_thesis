#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Constraints {
public:
	int getConstNum(){return const_n;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const = 0;
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const = 0;

	// The user just needs to implement JacobianIJV, so other methods can more efficiently concatenate multiple constraints jacobian
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x) const {
		Eigen::SparseMatrix<double> Jacobian(const_n, x.rows());
		auto IJV = JacobianIJV(x);
		Jacobian.setFromTriplets(IJV.begin(),IJV.end());
		return Jacobian;
	};

protected:
	int const_n;
};