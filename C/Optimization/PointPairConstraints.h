#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

class PointPairConstraints : public Constraints {
public:
	PointPairConstraints(const Eigen::VectorXi& b1, const Eigen::VectorXi& b2) : b1(b1),b2(b2) {
		const_n = b1.rows();
		IJV.resize(2*const_n);
	};

	virtual PointPairConstraints* clone() const {return new PointPairConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd constrained_pts_coords1; igl::slice(x,b1,1, constrained_pts_coords1);
		Eigen::VectorXd constrained_pts_coords2; igl::slice(x,b2,1, constrained_pts_coords2);
		return constrained_pts_coords1-constrained_pts_coords2;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_n = 0;
		for (int b_i = 0; b_i < b1.rows(); b_i++ ) {
			int var_const_idx1 = b1(b_i); int var_const_idx2 = b2(b_i);
			// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
			IJV[b_i] = Eigen::Triplet<double>(const_n, var_const_idx1, 1);
			IJV[b_i] = Eigen::Triplet<double>(const_n, var_const_idx2, -1);
			const_n++;
		}
	}
	
	// These are first order constraints and so the hessian is zero and there's no need in overriding the default zero hessian
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	Eigen::VectorXi b1,b2;
};