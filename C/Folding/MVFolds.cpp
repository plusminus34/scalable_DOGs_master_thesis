#include "MVFolds.h"

#include "../GeometryPrimitives/LinearTransformations.h"

using namespace std;

MountainValleyFold::MountainValleyFold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain, bool keep_rigid_motion) :
				 is_mountain(is_mountain), keep_rigid_motion(keep_rigid_motion) {
	const Eigen::MatrixXd& V = dog.getV();
	auto eS = dog.getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_idx];	
	ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_f = foldingCurve[edge_idx+1];
	dog.get_2_inner_vertices_from_edge(ep.edge,v1,v2);

	Eigen::RowVector3d p0 = ep.getPositionInMesh(V);
	len1 = (p0-V.row(v1)).norm();
	len2 = (p0-V.row(v2)).norm();

	Eigen::RowVector3d pf = ep_f.getPositionInMesh(V), pb = ep_b.getPositionInMesh(V);

	Eigen::RowVector3d t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
	Eigen::RowVector3d principal_n = ((pf-p0).normalized()+(pb-p0).normalized()).normalized();

	auto surface_norm_dot_z = t.cross(principal_n).dot(Eigen::RowVector3d(0,0,1));
	if ( surface_norm_dot_z < 0) {
		// Swap the direction of the curve
		std::swap(ep_b,ep_f);
		t = -1*t;
	}
	// The angle the curve folds makes with the folded edge
	auto norm_edge = (p0-V.row(v1))/len1;
	// the angle should be in the base of the tangent and principal_n
	curve_tangents_angle = acos(t.dot(norm_edge));
	auto dot = t.dot(norm_edge);
	auto t_coord = t.dot(norm_edge);
	auto n_coord = principal_n.dot(norm_edge);

	if (n_coord < 0) curve_tangents_angle = 2*M_PI-curve_tangents_angle;

	orig_center = p0;
	orig_t = t;
	orig_edge1 = -len1*norm_edge;
	orig_edge2 = len2*norm_edge;
	//std::cout << "orig_t = " << orig_t << std::endl;
}

void MountainValleyFold::get_constraint_indices(const Dog& dog, Eigen::VectorXi& b, EdgePoint& edgePoint) const {
	// Constraint the 2 vertices  along the curve and fix as the edge point itself
	int vnum = dog.getV().rows();
	b.resize(6); b << v1,v2,v1+vnum,v2+vnum,v1+2*vnum,v2+2*vnum;
	edgePoint = ep;
}

void MountainValleyFold::get_constraint_coords(double folding_angle, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords) const {
	const Eigen::MatrixXd& V = dog.getV();
	bc.resize(6);
	Eigen::RowVector3d p0;
	Eigen::RowVector3d edge1,edge2; Eigen::RowVector3d t;
	if (!keep_rigid_motion) {
		p0 = ep.getPositionInMesh(V);
		Eigen::RowVector3d pf = ep_f.getPositionInMesh(V), pb = ep_b.getPositionInMesh(V);
		//std::cout << "p0 = " << p0 << " but orig_center = " << orig_center << std::endl;
		t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
		Eigen::RowVector3d principal_n = ((pf-p0).normalized()+(pb-p0).normalized()).normalized();

		// First construct the unfolded t1 and t2 based on the frame t,n
		Eigen::RowVector3d t1 = cos(curve_tangents_angle)*t + sin(curve_tangents_angle)*principal_n;
		edge1 = -len1*t1;
		edge2 = len2*t1; // opposite direction
	} else {
		p0 = orig_center;
		t = orig_t;
		edge1 = orig_edge1;
		edge2 = orig_edge2;
	}

	double alpha = dihedral_angle_to_tangent_rotation_angle(*this, folding_angle); if (!is_mountain) alpha = -alpha;
	std::cout << "folding_angle = " << folding_angle << " alpha = " << alpha << std::endl;
	//double alpha = folding_angle; if (!is_mountain) alpha = -alpha;
	Eigen::RowVector3d center(0,0,0);

	// Now rotate them around the center in angle "folding_angle"
	//std::cout << "edge1 before = " << edge1 << std::endl;
	edge1 = rotate_vec(edge1, center, t, alpha);
	//std::cout << "edge1 after = " << edge1 << std::endl;
	edge2 = rotate_vec(edge2, center, t, -alpha);
	Eigen::RowVector3d dest_pos1 = p0+edge1;
	Eigen::RowVector3d dest_pos2 = p0+edge2;

	bc << dest_pos1(0),dest_pos2(0),dest_pos1(1),dest_pos2(1),dest_pos1(2),dest_pos2(2);
	edgeCoords.resize(1,3); edgeCoords.row(0) = p0;
}

double MountainValleyFold::dihedral_angle_to_tangent_rotation_angle(const MountainValleyFold& mvFold, double dihedral_angle) {
	// We rotate both tangents, so we need half
	return 0.5*acos(pow(cos(mvFold.curve_tangents_angle),2) + pow(sin(mvFold.curve_tangents_angle),2)*cos(dihedral_angle));
}

void MVFoldingConstraintsBuilder::add_fold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain, bool keep_center) {
	MountainValleyFold mvF(dog, curve_idx, edge_idx, is_mountain, keep_center);
	folds.push_back(mvF);
}

void MVFoldingConstraintsBuilder::get_folds_constraint_indices(const Dog& dog, Eigen::VectorXi& b, std::vector<EdgePoint>& edgePoints) const {
	int fold_n = folds.size(); b.resize(6*fold_n); edgePoints.resize(fold_n); // 2 vertices with 3 coordinates
	for (int i = 0; i < folds.size(); i++) {
		Eigen::VectorXi bFold; EdgePoint epFold; folds[i].get_constraint_indices(dog, bFold, epFold);
		edgePoints[i] = epFold;
	}
	int cnt = 0; // The order is x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < folds.size(); i++) {
			Eigen::VectorXi bFold; EdgePoint epFold; folds[i].get_constraint_indices(dog, bFold, epFold);
			b(cnt++) = bFold(2*coord_i); b(cnt++) = bFold(2*coord_i+1);
		}	
	}
}

void MVFoldingConstraintsBuilder::get_folds_constraint_coords(std::vector<double>& folding_angles, const Dog& dog, Eigen::VectorXd& bc,
																 Eigen::MatrixXd& edgeCoords) const {
	int fold_n = folds.size(); bc.resize(6*fold_n); edgeCoords.resize(fold_n,3); // 2 vertices with 3 coordinates
	for (int i = 0; i < folds.size(); i++) {
		Eigen::VectorXd bcFold; Eigen::MatrixXd edgeCoord; folds[i].get_constraint_coords(folding_angles[i], dog, bcFold, edgeCoord);
		edgeCoords.row(i) = edgeCoord.row(0);
	}
	int cnt = 0; // The order is x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < folds.size(); i++) {
			Eigen::VectorXd bcFold; Eigen::MatrixXd edgeCoord; folds[i].get_constraint_coords(folding_angles[i], dog, bcFold, edgeCoord);
			bc(cnt++) = bcFold(2*coord_i); bc(cnt++) = bcFold(2*coord_i+1);
		}	
	}
}

void MVFoldingConstraintsBuilder::get_folds_constraint_coords(double folding_angle, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords) const {
	std::vector<double> folding_angles(folds.size(), folding_angle);
	get_folds_constraint_coords(folding_angles, dog, bc, edgeCoords);
}