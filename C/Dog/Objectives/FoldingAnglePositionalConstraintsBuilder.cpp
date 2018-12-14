#include "FoldingAnglePositionalConstraintsBuilder.h"

#include "../../GeometryPrimitives/LinearTransformations.h"

using namespace std;
using namespace Eigen;

FoldingAnglePositionalConstraintsBuilder::FoldingAnglePositionalConstraintsBuilder(const Eigen::MatrixXd& V, const DogEdgeStitching& eS,
					const double& alpha) : alpha(alpha) {
	if (eS.edge_const_1.size()) {
		const_n = 9; // Constrain 4 vertices = 12 vars

		//int c_i = 0.5*eS.edge_const_1.size()/4-3; // TODO: This logic should be inside the constraints builder..
		int c_i = eS.edge_const_1.size()/2;

		edge1 = eS.edge_const_1[c_i]; edge2 = eS.edge_const_2[c_i]; edge_t_coordinate = eS.edge_coordinates[c_i];
		edge1_p1 = V.row(edge1.v1);
		edge1_p2 = V.row(edge1.v2);
		edge2_p1 = V.row(edge2.v1);
		edge2_p2 = V.row(edge2.v2);

		// Rotation center (the point)

		center = EdgePoint(edge1,edge_t_coordinate).getPositionInMesh(V);//edge_t_coordinate * edge1_p1 + (1-edge_t_coordinate) * edge1_p2;
		
		//axis = getRotAxis(V, es, center);
		setRotAxis(V,eS,center);
		
		int vnum = V.rows();
		b.resize(const_n);
		b << edge1.v1,edge2.v1,edge2.v2,
			edge1.v1+vnum,edge2.v1+vnum,edge2.v2+vnum,
			edge1.v1+2*vnum,edge2.v1+2*vnum,edge2.v2+2*vnum;
	}
}

void FoldingAnglePositionalConstraintsBuilder::setRotAxis(const Eigen::MatrixXd& V, const DogEdgeStitching& eS, 
				const Eigen::RowVector3d& rotCenter) {

	// Find a point a long a curve whose vertex is rotCenter
	for (auto curve: eS.stitched_curves) {
		for (int i = 1; i < curve.size()-1; i++) {
			if (curve[i].getPositionInMesh(V) == rotCenter) {
				auto pxb = curve[i-1].getPositionInMesh(V), p0 = rotCenter, pxf = curve[i+1].getPositionInMesh(V);
				auto e1 = pxf-p0, e2 = p0-pxb;

				// If the curve is locally straight rotate the edge in the plane
				if (e1.cross(e2).squaredNorm() < 1e-12) {
					// Rotate 90 deg in the plane
					//axis = RowVector3d(e1(1),-e1(0),0.).normalized();
					axis = e1.normalized();
				} else {
					auto tangent = (e1.squaredNorm()*e2 + e2.squaredNorm()*e1).normalized();
					axis = tangent.normalized();
					//auto circle_normal = (e1-e1.dot(tangent)*tangent).normalized();
					//axis = circle_normal;	
				}
				return;
			}
		}
	}
}

void FoldingAnglePositionalConstraintsBuilder::get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) {
	// Rotate edge 1 pt1 around the center (constraint point) in 'axis' with angle alpha
	Eigen::RowVector3d e1_p1_new_loc = rotate_vec(edge1_p1, center, axis, alpha);

	bc_out.resize(const_n);
	bc_out << e1_p1_new_loc(0),edge2_p1(0),edge2_p2(0),
			e1_p1_new_loc(1),edge2_p1(1),edge2_p2(1),
			e1_p1_new_loc(2),edge2_p1(2),edge2_p2(2);
	b_out = b;
}