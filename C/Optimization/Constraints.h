#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Constraints {
public:
	int getConstNum(){return const_n;}
	int getApproxNonZeros() {return approx_nnz;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const = 0;
	// By default returns a null matrix
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	// The user just needs to implement JacobianIJV, so other methods could efficiently concatenate multiple constraints jacobian
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x) const {
		Eigen::SparseMatrix<double> Jacobian(const_n, x.rows());
		auto IJV = JacobianIJV(x);
		Jacobian.setFromTriplets(IJV.begin(),IJV.end());
		return Jacobian;
	};

protected:
	int const_n;
	int approx_nnz;
};