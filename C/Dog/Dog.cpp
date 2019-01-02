#include "Dog.h"

#include <igl/boundary_loop.h>

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
}

Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F),V_ren(V), F_ren(Fsqr_to_F(F)) {
	submeshVSize.push_back(V.rows()); submeshFSize.push_back(F.rows());
	vi_to_submesh.assign(V.rows(),0);
	quad_topology(V,F,quadTop);
}

Dog::Dog(const Dog& d) : V(d.V),F(d.F),quadTop(d.quadTop),V_ren(d.V_ren), F_ren(d.F_ren), submeshVSize(d.submeshVSize), submeshFSize(d.submeshFSize),
						submesh_adjacency(d.submesh_adjacency), edgeStitching(d.edgeStitching),
						vi_to_submesh(d.vi_to_submesh) {
	// empty on purpose
}

Dog* Dog::get_submesh(int submesh_i) {
	if (submesh_i >= get_submesh_n()) return NULL;
	int submesh_v_min_i, submesh_v_max_i;
	get_submesh_min_max_i(submesh_i, submesh_v_min_i, submesh_v_max_i, true);
	Eigen::MatrixXd submeshV = V.block(submesh_v_min_i,0,submesh_v_max_i-submesh_v_min_i+1,3);

	int submesh_f_min_i, submesh_f_max_i;
	get_submesh_min_max_i(submesh_i, submesh_f_min_i, submesh_f_max_i, false);
	Eigen::MatrixXi submeshF = F.block(submesh_f_min_i,0,submesh_f_max_i-submesh_f_min_i+1,4);
	Dog* submeshDog = new Dog(submeshV, submeshF);
	return submeshDog;
}

void Dog::update_rendering_v() {
	V_ren_from_V_and_const(V,edgeStitching,V_ren);
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
