#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../../GeometryPrimitives/Curve.h"

class CurveInterpolationConstraintsBuilder {
public:
	CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			const double& timestep);
	~CurveInterpolationConstraintsBuilder();

	void get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc);
private:
	int const_n;
	Curve* srcCurve;
	Curve* dstCurve;
	// Rigid motion
	Eigen::RowVector3d T; Eigen::Matrix3d F;
	// between 0 and 1
	const double& timestep;
};
