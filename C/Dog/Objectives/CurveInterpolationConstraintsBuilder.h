#pragma once

#include "../Dog.h"
#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../../GeometryPrimitives/Curve.h"

class CurveInterpolationConstraintsBuilder {
public:
	CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, int curve_idx,
			const double& timestep, double k_translation, double k_addition, double t_addition, int max_curve_points);
	CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const std::vector<int>& v_indices,
			const double& timestep, double k_addition, double k_mult, double t_addition);
	~CurveInterpolationConstraintsBuilder();
	void get_curve_constraints(SurfaceCurve& surfaceCurve_out, Eigen::MatrixXd& bc);
private:
	void init_from_surface_curve(const Eigen::MatrixXd& V, double k_addition, double k_mult, double t_addition);
	int const_n;
	// The constraint (intrinsic) surface cure given by edge point
	SurfaceCurve surfaceCurve;
	Curve* srcCurve;
	Curve* dstCurve;
	// Rigid motion
	Eigen::RowVector3d T; Eigen::Matrix3d F;
	// between 0 and 1
	const double& timestep;
	int max_curve_points;
};
