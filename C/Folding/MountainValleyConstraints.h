#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class MountainValleyConstraints: public Constraints {
public:
	MountainValleyConstraints(const Eigen::MatrixXd& initV, const DogEdgeStitching& edgeStitching,
								std::vector<bool> is_mountain);
	virtual MountainValleyConstraints* clone() const {return new MountainValleyConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// TODO
	};

private:
	const DogEdgeStitching& eS;
	std::vector<bool> is_mountain;
	std::vector<std::vector<bool>> flip_binormal;
};
