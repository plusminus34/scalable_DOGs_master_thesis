#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"

class CompositeConstraints : public Constraints {
public:
	CompositeConstraints() {}
	CompositeConstraints(const std::vector<const Constraints&>& constraints) : constraints(constraints){
		approx_nnz = 0; const_n = 0; 
		for (auto cnst: constraints) {const_n+=cnst.getConstNum(); approx_nnz+= cnst.getApproxNonZeros();}
	};

	void add_constraints(const Constraints& cnst) {constraints.push_back(cnst);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		std::cout << "before vals = " << std::endl;
		Eigen::VectorXd vals(const_n);
		int const_cnt = 0; 
		for (auto cnst: constraints) {
			auto cnst_vals = cnst.Vals(x);
			for (int val_const_i = 0; val_const_i < cnst_vals.rows(); val_const_i++) {vals[const_cnt++] = cnst_vals[val_const_i];}
		}
		std::cout << "after vals = " << std::endl;
		return vals;
	}
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const {
		std::cout << "size here = " << constraints[0].JacobianIJV(x).size() << std::endl;
		std::cout << "before IJV = " << std::endl;
		std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz);
		int const_cnt = 0;
		for (auto cnst: constraints) {
			std::vector<Eigen::Triplet<double> > cnst_IJV = cnst.JacobianIJV(x);
			std::cout <<"cnst_IJV.size() = " << cnst_IJV.size() << std::endl;
			std::cout << "ok?"; int blabla; std::cin >> blabla;
			
			for (auto val : cnst_IJV) {
				int i = const_cnt; int j = val.col(); double value = val.value();
				std::cout << "adding i = " << i << " j = " << j << " value = " << value << std::endl;
				IJV.push_back(Eigen::Triplet<double>(i, j, value));
				const_cnt++;
			}
		}
		std::cout << "after IJV = " << std::endl;
		return IJV;
	}
	
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		Eigen::SparseMatrix<double> lambdaHessian(x.rows(),x.rows());
		for (auto cnst: constraints) lambdaHessian+=cnst.LambdaHessian(x,lambda);
		return lambdaHessian;
	};

private:
	std::vector<const Constraints&> constraints;
};