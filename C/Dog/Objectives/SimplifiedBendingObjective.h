#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class SimplifiedBendingObjective: public Objective {
  
public:
	SimplifiedBendingObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {};
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;
	Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const;

private:
	const QuadTopology& quadTop;
};