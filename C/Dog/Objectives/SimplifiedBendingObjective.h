#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class SimplifiedBendingObjective: public Objective {
  
public:
	SimplifiedBendingObjective(const QuadTopology& quadTop)  : quadTop(quadTop) {};
	virtual SimplifiedBendingObjective* clone() const {return new SimplifiedBendingObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual std::vector<Eigen::Triplet<double> > hessianIJV(const Eigen::VectorXd& x) const;
	
	const QuadTopology& quadTop;
};