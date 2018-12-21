#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class DogConstraints: public Constraints {
public:
	DogConstraints(const QuadTopology& quadTop, bool offset_planar = true);

	virtual DogConstraints* clone() const {return new DogConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	//Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const;
private:
	const bool offset_planar;
	const QuadTopology& quadTop;
};