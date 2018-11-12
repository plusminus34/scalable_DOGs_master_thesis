#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class SimplifiedBendingObjective: public Objective {
  
public:
	SimplifiedBendingObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {};
	virtual double obj(const Eigen::VectorXd& x);
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x);
	Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x);

private:
	const QuadTopology& quadTop;
};