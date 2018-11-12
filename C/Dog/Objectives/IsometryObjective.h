#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

// comparing to squared length
class IsometryObjective: public Objective {
  
public:
	IsometryObjective(const QuadTopology& quadTop);
	virtual double obj(const Eigen::VectorXd& x);
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x);
	virtual void set_ref(const Eigen::VectorXd& x0);
	Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x);

	//Eigen::VectorXd refA,refB;
private:
	const QuadTopology& quadTop;
	Eigen::VectorXd refL;
};