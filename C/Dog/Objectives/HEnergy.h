#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class HEnergy: public Objective {
  
public:
	HEnergy(const QuadTopology& quadTop)  : quadTop(quadTop) {}
	virtual HEnergy* clone() const {return new HEnergy(*this);}
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	const QuadTopology& quadTop;
};