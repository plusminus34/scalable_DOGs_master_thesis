#include "CurveInterpolationConstraintsBuilder.h"

#include "igl/procrustes.h"

CurveInterpolationConstraintsBuilder::CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
			int curve_idx, const double& timestep, double k_addition, double k_mult, double t_addition) : timestep(timestep) {
	// Create initial curve and dest curve, save the initial frame
	 if (eS.edge_const_1.size()) {
	 	surfaceCurve.edgePoints = eS.stitched_curves[curve_idx];
	 	init_from_surface_curve(V, k_addition, k_mult, t_addition);
	 }
}

CurveInterpolationConstraintsBuilder::CurveInterpolationConstraintsBuilder(const Eigen::MatrixXd& V, const std::vector<int>& v_indices,
			const double& timestep) : timestep(timestep) {
	surfaceCurve.edgePoints.resize(v_indices.size());
	int vi_cnt = 0;
	for (auto& edgePt: surfaceCurve.edgePoints) {
		edgePt.edge.v1 = edgePt.edge.v2 = v_indices[vi_cnt++];
		edgePt.t = 1;
	}
	init_from_surface_curve(V,0,2,0);
}

void CurveInterpolationConstraintsBuilder::init_from_surface_curve(const Eigen::MatrixXd& V, 
		double k_addition, double k_mult, double t_addition) {
	auto initCoords = surfaceCurve.get_curve_coords(V);
	srcCurve = new Curve(initCoords);
	Curve::getTranslationAndFrameFromCoords(initCoords, T, F);
	if (isnan(F(1,1))) {
		F.setIdentity(); //Fix for straight x curve..
	}
	// todo save frame
	// todo create dst curve from the curve parameters
	std::vector<double> dst_len = srcCurve->len, dst_k = srcCurve->k, dst_t = srcCurve->t;
	for (auto& k: dst_k) {k = k*k_mult+k_addition;/*k+=0.15;k*=2;*/};
	for (auto& t: dst_t) t+=t_addition;
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
	//Eigen::MatrixXd bc_src = srcCurve->getCoords(T,F);

	//Eigen::MatrixXd R; Eigen::VectorXd t; double scale_dummy;
	//igl::procrustes(bc,bc_src,false,false,scale_dummy,R,t);
	//bc = (bc*R).rowwise() + t.transpose();
	//std::cout << "R = " << R << " t = " << t << std::endl;

	surfaceCurve_out = surfaceCurve;
}