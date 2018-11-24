#include "CurveInterpolationConstraintsBuilder.h"


CurveInterpolationConstraintsBuilder::CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			const double& timestep) : timestep(timestep) {
	// Create initial curve and dest curve, save the initial frame
	//SurfaceCurve surfaceCurve(eS.)
	//std::vector<Edge> curve_edges; 
	//..std::vector<double> edge_t;
}

void CurveInterpolationConstraintsBuilder::get_positional_constraints(Eigen::VectorXi& b, Eigen::VectorXd& bc) {
	// Interpolate with timestep
	// Return coordinates
}