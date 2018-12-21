#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>

class Constraints {
public:

	virtual Constraints* clone() const = 0;

	int getConstNum() const {return const_n;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) = 0;
	// By default returns a null matrix
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	int get_IJV_size() {return IJV.size();}

	double deviation(const Eigen::VectorXd& x) const {return Vals(x).squaredNorm();}

	const std::vector<Eigen::Triplet<double>>& JacobianIJV() {return IJV;}

	// The user just needs to implement JacobianIJV, so other methods could efficiently concatenate multiple constraints jacobian
	virtual Eigen::SparseMatrix<double> Jacobian(const Eigen::VectorXd& x) {
		updateJacobianIJV(x);
    	if (cachedJacobian.rows() == 0) {
    		cachedJacobian =  Eigen::SparseMatrix<double>(const_n, x.rows());
      		igl::sparse_cached_precompute(IJV, cached_ijv_data, cachedJacobian);
      		std::cout << "not cached" << std::endl;
    	} else {
      		igl::sparse_cached(IJV, cached_ijv_data, cachedJacobian);
      		std::cout << "yeaaahah" << std::endl;
    	}
    	return cachedJacobian;
	};

protected:
	int const_n;
	std::vector<Eigen::Triplet<double> > IJV;
private:	
	Eigen::SparseMatrix<double> cachedJacobian;
	Eigen::VectorXi cached_ijv_data;
};