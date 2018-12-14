#include "CurveInterpolationConstraintsBuilder.h"


CurveInterpolationConstraintsBuilder::CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			const double& timestep) : timestep(timestep) {
	// Create initial curve and dest curve, save the initial frame
	 surfaceCurve.edgePoints = eS.stitched_curves[0];
	auto initCoords = surfaceCurve.get_curve_coords(V);
	srcCurve = new Curve(initCoords);
	Curve::getTranslationAndFrameFromCoords(initCoords, T, F);
	// todo save frame
	// todo create dst curve from the curve parameters
	std::vector<double> dst_len = srcCurve->len, dst_k = srcCurve->k, dst_t = srcCurve->t;
	for (auto& k: dst_k) k*=4;
	for (auto& t: dst_t) t+=0.1;
	dstCurve = new Curve(dst_len, dst_k, dst_t);
}

CurveInterpolationConstraintsBuilder::~CurveInterpolationConstraintsBuilder() {
	delete srcCurve;
	delete dstCurve;
}

void CurveInterpolationConstraintsBuilder::get_curve_constraints(SurfaceCurve& surfaceCurve_out, Eigen::MatrixXd& bc) {
	// Interpolate with timestep
	Curve intCurve(*srcCurve,*dstCurve,timestep);
	// Return coordinates
	bc = intCurve.getCoords(T,F);
	surfaceCurve_out = surfaceCurve;
	// TODO this gives us linear constraints and we can't use b or bc here since these are edge points
	// So we need to indices at b and at bc

	// We can generalize the previous constraints for these linear constraints (plotting but also optimization)
	// TODO: write EdgePointConstraints? Then viewer should be able to show it and all optimization as well
	// (Besides ARAP now anc procrustes now)

	// The procrustes solver can do multiple rounds of ICP, where one mesh can be forced to not move
	// At every iteration we set positional constraints from the linear constraints
	// Assuming that other mesh is fixed this linear constraints become a standard positional constraint
	// This should be ok
	// How fast can this work?
	// Meaning we just need to write an ICP class (start with assuming 2 connected components)
}