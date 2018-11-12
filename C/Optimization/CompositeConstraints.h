#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"

class CompositeConstraints : public Constraints {
public:
	CompositeConstraints(const std::vector<Constraints*>& constraints) : constraints(constraints){
		approx_nnz = 0; const_n = 0; 
		for (auto cnst: constraints) {const_n+=cnst->getConstNum(); approx_nnz+= cnst->getApproxNonZeros();}
	};

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x){
		Eigen::VectorXd vals(const_n);
		int const_cnt = 0; 
		for (auto cnst: constraints) {
			auto cnst_vals = cnst->Vals(x);
			for (int val_const_i = 0; val_const_i < cnst_vals.rows(); val_const_i++) {vals[const_cnt++] = cnst_vals[val_const_i];}
		}
		return vals;
	}
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const {
		std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz);
		int const_cnt = 0;
		for (auto cnst: constraints) {
			auto cnst_IJV = cnst->JacobianIJV(x);
			for (auto val : cnst_IJV) IJV.push_back(Eigen::Triplet<double>(const_cnt++, val.col(), val.value()));
		}
	}
	
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		Eigen::SparseMatrix<double> lambdaHessian(x.rows(),x.rows());
		for (auto cnst: constraints) lambdaHessian+=cnst->LambdaHessian(x,lambda);
		return lambdaHessian;
	};

private:
	std::vector<Constraints*> constraints;
};