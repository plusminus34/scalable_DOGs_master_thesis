#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class StitchingConstraints: public Constraints {
public:
	StitchingConstraints(const QuadTopology& quadTop,const DogEdgeStitching& edgeStitching);
	virtual StitchingConstraints* clone() const {return new StitchingConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

private:
	const QuadTopology& quadTop;
	const DogEdgeStitching& eS;
};
