#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class EqualDiagObjective: public Objective {
  
public:
	EqualDiagObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {}
	virtual EqualDiagObjective* clone() const {return new EqualDiagObjective(*this);}
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	const QuadTopology& quadTop;
};