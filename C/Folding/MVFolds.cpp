#include "MVFolds.h"

#include "../GeometryPrimitives/LinearTransformations.h"

using namespace std;

MountainValleyFold::MountainValleyFold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain) : is_mountain(is_mountain) {
	const Eigen::MatrixXd& V = dog.getV();
	auto eS = dog.getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_idx];	
	ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_f = foldingCurve[edge_idx+1];
	dog.get_2_inner_vertices_from_edge(ep.edge,v1,v2);



	Eigen::RowVector3d p0 = ep.getPositionInMesh(V);
	len1 = (p0-V.row(v1)).norm();
	len2 = (p0-V.row(v2)).norm();

	Eigen::RowVector3d pf = ep_f.getPositionInMesh(V), pb = ep_b.getPositionInMesh(V);
	Eigen::RowVector3d t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
	// The angle the curve folds makes with the folded edge
	curve_tangents_angle = acos(t.dot( (p0-V.row(v1))/len1 ));
}

void MountainValleyFold::get_constraint_indices(const Dog& dog, Eigen::VectorXi& b, EdgePoint& edgePoint) {
	// Constraint the 2 vertices  along the curve and fix as the edge point itself
	int vnum = dog.getV().rows();
	b.resize(6); b << v1,v1+vnum,v1+2*vnum,v2,v2+vnum,v2+2*vnum;
	edgePoint = ep;
}

void MountainValleyFold::get_constraint_coords(double folding_angle, const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords) {
	const Eigen::MatrixXd& V = dog.getV();
	bc.resize(6);

	// Find current tangent and principle normal (for the frame)
	Eigen::RowVector3d p0 = ep.getPositionInMesh(V), pf = ep_f.getPositionInMesh(V), pb = ep_b.getPositionInMesh(V);
	Eigen::RowVector3d t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
	Eigen::RowVector3d principal_n = ((pf-p0).normalized()+(pb-p0).normalized()).normalized();

	// First construct the unfolded t1 and t2 based on the frame t,n
	Eigen::RowVector3d t1 = cos(curve_tangents_angle)*t + sin(curve_tangents_angle)*principal_n;
	Eigen::RowVector3d edge1 = len1*t1;
	Eigen::RowVector3d edge2 = -len2*t1; // opposite direction

	double alpha = curve_tangents_angle; if (!is_mountain) alpha = -alpha;
	Eigen::RowVector3d center(0,0,0);

	// Now rotate them around the center in angle "folding_angle"
	edge1 = rotate_vec(edge1, center, t, alpha);
	edge2 = rotate_vec(edge2, center, t, -alpha);
	Eigen::RowVector3d dest_pos1 = center+edge1;
	Eigen::RowVector3d dest_pos2 = center+edge2;

	bc << dest_pos1(0),dest_pos1(1),dest_pos1(2),dest_pos2(0),dest_pos2(1),dest_pos2(2);


	edgeCoords.resize(1,3); edgeCoords.row(0) = ep.getPositionInMesh(V);
}

void MVFoldingConstraintsBuilder::add_fold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain) {
	MountainValleyFold mvF(dog, curve_idx, edge_idx, is_mountain);
	folds.push_back(mvF);
}
/*
struct MountainValleyFold {
	bool is_mountain;
	EdgePoint ep, ep_b, ep_f;
	int v1,v2; double len1, len2;
	double curve_tangents_angle;
};

class MVFoldingConstraintsBuilder {
public:
	void add_fold(const Dog& dog, int curve_idx, int edge_idx, bool is_mountain);
	int get_folds_num() {return folds.size();}
	void clear_folds() {folds.clear();}

	void get_folds_constraint_indices(const Dog& dog, Eigen::VectorXi& b, std::vector<EdgePoint>& edgePoints);
	void get_folds_constraint_coords(const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords);
private:
	std::vector<MountainValleyFold> folds;
};
*/