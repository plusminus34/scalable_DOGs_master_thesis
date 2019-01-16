#include "ReflectionFolds.h"

using namespace std;

ReflectionFold::ReflectionFold(const Dog& dog, int curve_idx, int edge_idx, std::vector<bool>& submeshes_set) {
	const Eigen::MatrixXd& V = dog.getV();
	auto eS = dog.getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_idx];	
	ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_f = foldingCurve[edge_idx+1];
	dog.get_2_inner_vertices_from_edge(ep.edge,v1,v2);

	if (!submeshes_set[dog.v_to_submesh_idx(v1)]) {
		if (!submeshes_set[dog.v_to_submesh_idx(v2)]) {
			std::cout << "Should not get here! Got to vertices of 2 submeshes we didn't pass on" << std::endl;
			std::cout << "dog.v_to_submesh_idx(v1) = " << dog.v_to_submesh_idx(v1) << std::endl;
			std::cout << "dog.v_to_submesh_idx(v2) = " << dog.v_to_submesh_idx(v2) << std::endl;
			exit(1);
		}
		std::swap(v1,v2);
	}

	Eigen::RowVector3d p0 = ep.getPositionInMesh(V);
	len2 = (p0-V.row(v2)).norm();
}

void ReflectionFold::get_constraint_indices(const Dog& dog, Eigen::VectorXi& b, EdgePoint& edgePoint) const {
	// Constraint the 2 vertices  along the curve and fix as the edge point itself
	int vnum = dog.getV().rows();
	b.resize(3); b << v2,v2+vnum,v2+2*vnum;
}

void ReflectionFold::get_constraint_coords(const Dog& dog, Eigen::VectorXd& bc, Eigen::MatrixXd& edgeCoords) const {
	const Eigen::MatrixXd& V = dog.getV();
	bc.resize(6);
	Eigen::RowVector3d p0;
	Eigen::RowVector3d edge1,edge2; Eigen::RowVector3d t;

	// Get Binormal
	p0 = ep.getPositionInMesh(V);
	Eigen::RowVector3d pf = ep_f.getPositionInMesh(V), pb = ep_b.getPositionInMesh(V);
	//std::cout << "p0 = " << p0 << " but orig_center = " << orig_center << std::endl;
	t = ((pf-p0).normalized()-(pb-p0).normalized()).normalized();
	Eigen::RowVector3d principal_n = ((pf-p0).normalized()+(pb-p0).normalized()).normalized();
	auto B = t.cross(principal_n);

	// Reflect the edge direction through the oscullating plane
	Eigen::RowVector3d edge_dir = (p0-V.row(v1)).normalized();
	std::cout << "edge_dir.dot(B) = " << edge_dir.dot(B) << std::endl;
	edge_dir = edge_dir- 2*edge_dir.dot(B)*B;
	Eigen::RowVector3d dest_pos = p0+len2*edge_dir;

	
	bc << dest_pos(0),dest_pos(1),dest_pos(2);
}

void ReflectionFoldConstraintsBuilder::add_fold(const Dog& dog, int curve_idx, int edge_idx, std::vector<bool>& submeshes_set) {
	ReflectionFold refF(dog, curve_idx, edge_idx, submeshes_set);
	folds.push_back(refF);
}

void ReflectionFoldConstraintsBuilder::get_folds_constraint_indices(const Dog& dog, Eigen::VectorXi& b) const {
	int fold_n = folds.size(); b.resize(3*fold_n);
	int cnt = 0; // The order is x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < folds.size(); i++) {
			Eigen::VectorXi bFold; EdgePoint epFold; folds[i].get_constraint_indices(dog, bFold, epFold);
			b(cnt++) = bFold(coord_i);
		}	
	}
}

void ReflectionFoldConstraintsBuilder::get_folds_constraint_coords(const Dog& dog, Eigen::VectorXd& bc) const {
	int fold_n = folds.size(); bc.resize(3*fold_n);
	int cnt = 0; // The order is x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < folds.size(); i++) {
			Eigen::VectorXd bcFold; Eigen::MatrixXd edgeCoord; folds[i].get_constraint_coords(dog, bcFold, edgeCoord);
			bc(cnt++) = bcFold(coord_i);
		}	
	}
}