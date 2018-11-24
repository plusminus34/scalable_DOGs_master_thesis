#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../../GeometryPrimitives/Curve.h"

class CurveInterpolationConstraintsBuilder {
public:
	CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			const double& timestep);

	void get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc);
private:

	int const_n;
	// between 0 and 1
	const double& timestep;
};
