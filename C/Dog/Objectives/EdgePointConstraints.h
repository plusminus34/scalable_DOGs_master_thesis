#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

class EdgePointConstraints : public Constraints {
public:
	EdgePointConstraints(std::vector<EdgePoint> edgePoints , const Eigen::MatrixXd& curveCoords) {
		const_n = 3*edgePoints.rows(); approx_nnz = 2*const_n;
	};

	virtual PositionalConstraints* clone() const {return new EdgePointConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::MatrixXd curveCoords(EdgePoint::getPositionInMesh(edgePoints, ))
		//Eigen::VectorXd constrained_pts_coords; igl::slice(x,b,1, constrained_pts_coords);
		// A vector of x(b)-bc for the constraints x(b)-bc = 0
		return constrained_pts_coords-bc;
	}
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const {
		std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz);
		int const_n = 0;
		for (int b_i = 0; b_i < b.rows(); b_i++ ) {
			int var_const_idx = b(b_i);
			// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
			IJV.push_back(Eigen::Triplet<double>(const_n++, var_const_idx, 1));
		}
		return IJV;
	}
	
	// These are first order constraints and so the hessian is zero and there's no need in overriding the default zero hessian
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	Eigen::VectorXi getPositionIndices() const {return b;}
	Eigen::VectorXd getPositionVals() const {return bc;}

private:
	std::vector<EdgePoint> edgePoints; Eigen::VectorXd bc;
};