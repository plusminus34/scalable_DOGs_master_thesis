#include "CurveInterpolationConstraintsBuilder.h"


CurveInterpolationConstraintsBuilder::CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			const double& timestep) : timestep(timestep) {
	// Create initial curve and dest curve, save the initial frame
	SurfaceCurve surfaceCurve; surfaceCurve.edgePoints = eS.stitched_curves[0];
	srcCurve = new Curve(surfaceCurve.get_curve_coords(V));
	// todo save frame
	// todo create dst curve from the curve parameters
	// todo constructor should deallocate the curves
	// todo interpolate should just interpolate curves and get the coordinates with the frame
	//..std::vector<double> edge_t;
}

void CurveInterpolationConstraintsBuilder::get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc) {
	// Interpolate with timestep
	// Return coordinates
}