#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"

class CompositeConstraints : public Constraints {
public:
	CompositeConstraints() {}
	virtual CompositeConstraints* clone() const {return new CompositeConstraints(*this);}

	CompositeConstraints(const std::vector<Constraints*>& constraints_i) {
		constraints.resize(constraints_i.size());
		for (int i = 0; i < constraints.size(); i++) constraints[i] = constraints_i[i]->clone();
		approx_nnz = 0; const_n = 0; 
		for (auto cnst: constraints) {const_n+=cnst->getConstNum(); approx_nnz+= cnst->getApproxNonZeros();}
	};

	void add_constraints(Constraints* cnst) {
		constraints.push_back(cnst->clone());
		const_n+= cnst->getConstNum();
		approx_nnz+= cnst->getApproxNonZeros();
	}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(const_n);
		int const_cnt = 0; 
		for (auto cnst: constraints) {
			auto cnst_vals = cnst->Vals(x);
			for (int val_const_i = 0; val_const_i < cnst_vals.rows(); val_const_i++) {vals[const_cnt++] = cnst_vals[val_const_i];}
		}
		if (const_n!= const_cnt) {
			std::cout << "Error! const_n = " << const_n << " and const_cnt = " << const_cnt << std::endl;
			exit(1);
		}
		return vals;
	}
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const {
		std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz);
		int row_base = 0;
		for (auto cnst: constraints) {
			std::vector<Eigen::Triplet<double> > cnst_IJV = cnst->JacobianIJV(x);
			
			for (auto val : cnst_IJV) {
				IJV.push_back(Eigen::Triplet<double>(val.row() + row_base, val.col(), val.value()));
			}
			row_base += cnst->getConstNum();
		}
		return IJV;
	}
	
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		Eigen::SparseMatrix<double> lambdaHessian(x.rows(),x.rows());
		for (auto cnst: constraints) lambdaHessian+=cnst->LambdaHessian(x,lambda);
		return lambdaHessian;
	};

private:
	std::vector<Constraints*> constraints;
};