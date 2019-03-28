#include "Dog.h"

#include <igl/boundary_loop.h>

using namespace std;

int DogEdgeStitching::get_vertex_edge_point_deg(Edge& edge) const {
	int mult_edge_index = edge_to_duplicates.at(edge);
	int mult_edge_const_start = multiplied_edges_start[mult_edge_index]; 
	int mult_edges_num = multiplied_edges_num[mult_edge_index];
	return mult_edges_num;
}

Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, DogEdgeStitching edgeStitching, 
		const Eigen::MatrixXd& V_ren, const Eigen::MatrixXi& F_ren, std::vector<int> submeshVSize, std::vector<int> submeshFSize,
		const std::vector< std::vector<int> >& submesh_adjacency) :
				V(V),F(F),edgeStitching(edgeStitching),V_ren(V_ren), F_ren(F_ren), submeshVSize(submeshVSize), submeshFSize(submeshFSize),
				submesh_adjacency(submesh_adjacency) {

	// set mesh_min_max_i;
	vi_to_submesh.resize(V.rows());
	int min_idx = 0;
	for (int sub_i = 0; sub_i < submeshVSize.size(); sub_i++) {
		int submesh_vn = submeshVSize[sub_i]; int max_idx = min_idx +submesh_vn-1;
		for (int subm_vi = min_idx; subm_vi <= max_idx; subm_vi++) vi_to_submesh[subm_vi] = sub_i;
		min_idx += submesh_vn;
	}
	quad_topology(V,F,quadTop);
	setup_stitched_curves_initial_l_angles_length();
}
 
Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F),V_ren(V), F_ren(Fsqr_to_F(F)) {
	submeshVSize.push_back(V.rows()); submeshFSize.push_back(F.rows());
	vi_to_submesh.assign(V.rows(),0);
	quad_topology(V,F,quadTop);
}

Dog::Dog(const Dog& d) : V(d.V),F(d.F),quadTop(d.quadTop),V_ren(d.V_ren), F_ren(d.F_ren), submeshVSize(d.submeshVSize), submeshFSize(d.submeshFSize),
						submesh_adjacency(d.submesh_adjacency), edgeStitching(d.edgeStitching),
						stitched_curves_l(d.stitched_curves_l),stitched_curves_angles(d.stitched_curves_angles), stitched_curves_curvature(d.stitched_curves_curvature),
						vi_to_submesh(d.vi_to_submesh) {
	// empty on purpose
}

void Dog::setup_stitched_curves_initial_l_angles_length() {
	const DogEdgeStitching& eS = edgeStitching;
	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		std::vector<double> curve_l;
		// Go through all the inner vertices in a curve
		for (int edge_idx = 0; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			auto pt1 = foldingCurve[edge_idx].getPositionInMesh(V); auto pt2 = foldingCurve[edge_idx+1].getPositionInMesh(V);
			curve_l.push_back((pt1-pt2).norm());
		}
		std::vector<double> curve_a,curve_k;
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			auto pt = foldingCurve[edge_idx].getPositionInMesh(V); auto pt_b = foldingCurve[edge_idx-1].getPositionInMesh(V); auto pt_f = foldingCurve[edge_idx+1].getPositionInMesh(V);
			auto e1 = pt_f-pt; auto e2 = pt-pt_b;
			double edges_dot = (e1).dot(e2); double l_normalization = e1.norm()*e2.norm();
			double angle = acos(edges_dot/l_normalization);
			double k = sin(angle)/(pt_f-pt_b).norm();

			curve_a.push_back(angle);
			curve_k.push_back(k);
		}

		stitched_curves_l.push_back(curve_l);
		stitched_curves_angles.push_back(curve_a);
		stitched_curves_curvature.push_back(curve_k);
	}
}

Dog* Dog::get_submesh(int submesh_i) {
	if (submesh_i >= get_submesh_n()) return NULL;
	int submesh_v_min_i, submesh_v_max_i;
	get_submesh_min_max_i(submesh_i, submesh_v_min_i, submesh_v_max_i, true);
	Eigen::MatrixXd submeshV = V.block(submesh_v_min_i,0,submesh_v_max_i-submesh_v_min_i+1,3);

	int submesh_f_min_i, submesh_f_max_i;
	get_submesh_min_max_i(submesh_i, submesh_f_min_i, submesh_f_max_i, false);
	Eigen::MatrixXi submeshF = F.block(submesh_f_min_i,0,submesh_f_max_i-submesh_f_min_i+1,4);
	// substract submesh_f_min_i from all F faces
	for (int i = 0; i < submeshF.rows(); i++) {submeshF.row(i).array() -= submesh_v_min_i;}


	Dog* submeshDog = new Dog(submeshV, submeshF);
	return submeshDog;
}

void Dog::update_submesh_V(int submesh_i, const Eigen::MatrixXd& submeshV) {
	int submesh_v_min_i, submesh_v_max_i;
	get_submesh_min_max_i(submesh_i, submesh_v_min_i, submesh_v_max_i, true);
	for (int i = 0; i < submeshV.rows(); i++) {V.row(submesh_v_min_i + i) = submeshV.row(i);}
	//update_rendering_v();
	update_Vren2();
}

void Dog::update_rendering_v() {
	//V_ren_from_V_and_const(V,edgeStitching,V_ren);
	update_Vren2();
}

void Dog::get_2_submeshes_vertices_from_edge(const Edge& edge, int &v1_out, int &v2_out, int &w1_out, int& w2_out) const {
	int mult_edge_start = edgeStitching.multiplied_edges_start[edgeStitching.edge_to_duplicates.at(edge)];
	Edge e1 = edgeStitching.edge_const_1[mult_edge_start], e2 = edgeStitching.edge_const_2[mult_edge_start];
	v1_out = e1.v1; v2_out = e1.v2;
	w1_out = e2.v1; w2_out = e2.v2;
}

void Dog::get_2_inner_vertices_from_edge(const Edge& edge, int &v1_out, int &v2_out) const {
	int mult_edge_start = edgeStitching.multiplied_edges_start[edgeStitching.edge_to_duplicates.at(edge)];
	Edge e1 = edgeStitching.edge_const_1[mult_edge_start], e2 = edgeStitching.edge_const_2[mult_edge_start];
	if (!quadTop.is_bnd_v[e1.v1]) v1_out = e1.v1; else v1_out = e1.v2;
	if (!quadTop.is_bnd_v[e2.v1]) v2_out = e2.v1; else v2_out = e2.v2;
}

void Dog::V_ren_from_V_and_const(const Eigen::MatrixXd& V, const DogEdgeStitching& fC, Eigen::MatrixXd& V_ren) {
	int consts_num = fC.edge_coordinates.size();
	Eigen::MatrixXd V_folds_polygons(consts_num,3);
	for (int const_i = 0; const_i < consts_num; const_i++) {
		double t = fC.edge_coordinates[const_i];
		V_folds_polygons.row(const_i) = t*V.row(fC.edge_const_1[const_i].v1) + (1-t)*V.row(fC.edge_const_1[const_i].v2);
	}
	V_ren.resize(V.rows()+consts_num,3);
	V_ren << V,V_folds_polygons;
}

void Dog::update_Vren2() {
	int subm_n = edgeStitching.submesh_to_edge_pt.size();
	// Check if there is only 1 submesh
	if (subm_n < 2) {
		V_ren = V; return;
	}
	// go through the entire submeshes
	int vi_cnt = 0;
	for (int subi = 0; subi < subm_n; subi++) {
		// First add the normal vertices, that are on the DOG but are not on the crease points
		int subm_min,subm_max; get_submesh_min_max_i(subi, subm_min, subm_max);
		for (int vi = subm_min; vi <= subm_max; vi++) V_ren.row(vi_cnt++) = V.row(vi);
		// Then add the crease points
		for (auto ei : edgeStitching.submesh_to_edge_pt[subi]) {
			// Get the vertex
			double t = edgeStitching.edge_coordinates[ei];
			V_ren.row(vi_cnt++) = t*V.row(edgeStitching.edge_const_1[ei].v1) + (1-t)*V.row(edgeStitching.edge_const_1[ei].v2);
		}
	}
}

void Dog::get_submesh_min_max_i(int submesh_i, int& submesh_min_i, int& submesh_max_i, bool vertices) {
	std::vector<int> idx_list;
	if (vertices) 
		idx_list = submeshVSize;
	else
		idx_list = submeshFSize; 

	int sub_i = 0; submesh_min_i = 0; int submesh_vn = idx_list[0];
	for (int sub_i = 1; sub_i <= submesh_i; sub_i++) {
		submesh_min_i += submesh_vn;
		submesh_vn = idx_list[sub_i];
	}
	submesh_max_i = submesh_min_i+submesh_vn-1;
}

int Dog::v_ren_idx_to_v_idx(int v_idx) const {
	auto coords = V_ren.row(v_idx);
	for (int i = 0; i < V.rows(); i++) {
		if ((V.row(i)-coords).norm() == 0) 
			return i;
	}
	// did not find matching point
	return -1;
}

int Dog::v_ren_idx_to_edge(int v_idx, EdgePoint& edgePt) const {
	auto coords = V_ren.row(v_idx);
	for (auto curve_v: edgeStitching.stitched_curves) {
		for (auto ep : curve_v) {
			if ((ep.getPositionInMesh(V)-coords).norm() == 0) {
				edgePt = ep;
				return 0;
			}
		}
	}
	return -1;
	std::vector<std::vector<EdgePoint>> stitched_curves;
}

bool Dog::is_crease_vertex_flat(int curve_i, int edge_i) const {
	double flat_point_curvature_tolerance = 1e-5;
	// First verify the range
	if ( (curve_i >=0 ) && (curve_i <= edgeStitching.stitched_curves.size()) ) {
		if ( (edge_i >= 1) && (edge_i <= edgeStitching.stitched_curves[curve_i].size()-1) ) {
			//cout << "abs(stitched_curves_curvature[curve_i][edge_i-1] = " << abs(stitched_curves_curvature[curve_i][edge_i-1]) << endl;
			return (abs(stitched_curves_curvature[curve_i][edge_i-1]) <= flat_point_curvature_tolerance);
		}
	}
	return true;
}
