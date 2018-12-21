#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>

class Constraints {
public:

	virtual Constraints* clone() const = 0;

	int getConstNum() const {return const_n;}
	int getApproxNonZeros() const {return approx_nnz;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const = 0;
	// By default returns a null matrix
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	double deviation(const Eigen::VectorXd& x) const {return Vals(x).squaredNorm();}

	// The user just needs to implement JacobianIJV, so other methods could efficiently concatenate multiple constraints jacobian
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x) {
		Eigen::SparseMatrix<double> Jacobian(const_n, x.rows());
		auto IJV = JacobianIJV(x);
    	if (cached_ijv_data.rows() == 0) {
      		igl::sparse_cached_precompute(IJV, cached_ijv_data, Jacobian);
    	} else {
      		igl::sparse_cached(IJV, cached_ijv_data, Jacobian);
    	}
    	return Jacobian;
	};

protected:
	int const_n;
	int approx_nnz;
private:
	Eigen::VectorXi cached_ijv_data;
};