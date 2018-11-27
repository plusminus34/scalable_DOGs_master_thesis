#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"
#include "../QuadMesh/Quad.h"

class EdgePointConstraints : public Constraints {
public:
	EdgePointConstraints(std::vector<EdgePoint> edgePoints , const Eigen::MatrixXd& edgePointCoords) : edgePoints(edgePoints) {
		mat2_to_vec(edgePointCoords, bc); // flatten to a vector
		const_n = bc.rows(); approx_nnz = 2*const_n; // 2 points per edge
	};
	EdgePointConstraints() {const_n = 0; approx_nnz = 0; edgePoints.resize(0); bc.resize(0);} // empty set of constraints c'tor

	virtual EdgePointConstraints* clone() const {return new EdgePointConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd edgeCoords(EdgePoint::getPositionInMesh(edgePoints, x));
		return edgeCoords-bc;
	}

	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const {
		int vn = x.rows()/3;
		std::vector<Eigen::Triplet<double> > IJV; IJV.reserve(approx_nnz);
		int const_n = 0;
		for (int b_i = 0; b_i < edgePoints.size(); b_i++ ) {
			int v1 = edgePoints[b_i].edge.v1, v2 = edgePoints[b_i].edge.v2;
			double t = edgePoints[b_i].t;

			// 1 constraint for every coordinate
			IJV.push_back(Eigen::Triplet<double>(const_n, v1, t));
			IJV.push_back(Eigen::Triplet<double>(const_n, v2, 1-t));

			IJV.push_back(Eigen::Triplet<double>(const_n+1, vn+v1, t));
			IJV.push_back(Eigen::Triplet<double>(const_n+1, vn+v2, 1-t));

			IJV.push_back(Eigen::Triplet<double>(const_n+2, 2*vn+v1, t));
			IJV.push_back(Eigen::Triplet<double>(const_n+2, 2*vn+v2, 1-t));

			const_n+=3;
		}
		return IJV;
	}
	
	// These are first order constraints and so the hessian is zero and there's no need in overriding the default zero hessian
	virtual Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const {
		return Eigen::SparseMatrix<double>(x.rows(),x.rows());
	};

	std::vector<EdgePoint> getEdgePoints() const {return edgePoints;}
	Eigen::VectorXd getEdgePointConstraints() const {return bc;}

private:
	std::vector<EdgePoint> edgePoints; Eigen::VectorXd bc;
};