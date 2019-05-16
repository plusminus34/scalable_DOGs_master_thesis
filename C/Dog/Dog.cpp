#include "Dog.h"

#include <igl/boundary_loop.h>
#include <igl/edges.h>

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
				V(V),F(F),flatV(V), edgeStitching(edgeStitching),V_ren(V_ren), F_ren(F_ren), submeshVSize(submeshVSize), submeshFSize(submeshFSize),
				submesh_adjacency(submesh_adjacency) {
	cout << "Dog::Dog" << endl;
	// set mesh_min_max_i;
	vi_to_submesh.resize(V.rows());
	int min_idx = 0;
	for (int sub_i = 0; sub_i < submeshVSize.size(); sub_i++) {
		int submesh_vn = submeshVSize[sub_i]; int max_idx = min_idx +submesh_vn-1;
		for (int subm_vi = min_idx; subm_vi <= max_idx; subm_vi++) vi_to_submesh[subm_vi] = sub_i;
		min_idx += submesh_vn;
	}
	cout << "setting neighborhood structures" << endl;
	quad_topology(V,F,quadTop);
	cout << "setting curve stuff" << endl;
	setup_stitched_curves_initial_l_angles_length();
	cout << "setting wireframe edges" << endl;
	setup_rendered_wireframe_edges_from_planar();
	cout << "setting up uv" << endl;
	setup_uv_and_texture();
	cout << "DOG setup complete" << endl;
}
 
Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F),flatV(V), V_ren(V), F_ren(Fsqr_to_F(F)) {
	submeshVSize.push_back(V.rows()); submeshFSize.push_back(F.rows());
	vi_to_submesh.assign(V.rows(),0);
	quad_topology(V,F,quadTop);
	setup_rendered_wireframe_edges_from_planar();
}

Dog::Dog(const Dog& d) : V(d.V),F(d.F),flatV(d.flatV),quadTop(d.quadTop),V_ren(d.V_ren), F_ren(d.F_ren), rendered_wireframe_edges(d.rendered_wireframe_edges),
						submeshVSize(d.submeshVSize), submeshFSize(d.submeshFSize),
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
			//std::cout << "i = " << edge_idx << " with angle = " << angle << ", curvature " << k << " l1 = " << e1.norm() << " l2 = " << e2.norm() << std::endl;

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
	update_Vren();
}

void Dog::update_rendering_v() {
	V_ren = V;
	//V_ren_from_V_and_const(V,edgeStitching,V_ren);
	//update_Vren();
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

void Dog::update_Vren() {
	int subm_n = submeshVSize.size();
	// go through the entire submeshes
	int vi_cnt = 0;
	for (int subi = 0; subi < subm_n; subi++) {
		// First add the normal vertices, that are on the DOG but are not on the crease points
		int subm_min,subm_max; get_submesh_min_max_i(subi, subm_min, subm_max);
		for (int vi = subm_min; vi <= subm_max; vi++) V_ren.row(vi_cnt++) = V.row(vi);
		// Then add the crease points
		if (edgeStitching.submesh_to_edge_pt.size()){
			for (auto ei : edgeStitching.submesh_to_edge_pt[subi]) {
				// Get the vertex
				double t = edgeStitching.edge_coordinates[ei];
				V_ren.row(vi_cnt++) = t*V.row(edgeStitching.edge_const_1[ei].v1) + (1-t)*V.row(edgeStitching.edge_const_1[ei].v2);
			}
			for (int j = 0; j < edgeStitching.submesh_to_bnd_edge[subi].size(); j++) {
				EdgePoint ep = edgeStitching.submesh_to_bnd_edge[subi][j];
				V_ren.row(vi_cnt++) = ep.getPositionInMesh(V);
			//std::cout << "Added row = " << V_ren_list[i].row(submeshVList[i].rows() + eS.submesh_to_edge_pt[i].size() + j)  << std::endl;
			}
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
	double eps = 1e-10;
	auto coords = V_ren.row(v_idx);
	for (int i = 0; i < V.rows(); i++) {
		if ((V.row(i)-coords).norm() < eps) 
			return i;
	}
	// did not find matching point
	std::cout << "could not find point!" << std::endl;
	return -1;
}

int Dog::v_ren_idx_to_edge(int v_idx, EdgePoint& edgePt) const {
	double eps = 1e-10;
	auto coords = V_ren.row(v_idx);
	for (auto curve_v: edgeStitching.stitched_curves) {
		for (auto ep : curve_v) {
			if ((ep.getPositionInMesh(V)-coords).norm() < eps) {
				edgePt = ep;
				return 0;
			}
		}
	}
	return -1;
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

void Dog::setup_rendered_wireframe_edges_from_planar() {
	// Hack but will work for now
	double eps = 1e-10;
	Eigen::MatrixXi E; igl::edges(F_ren,E);
	for (int i = 0; i < E.rows(); i++) {
		// make sure the edge is an 'x' or 'y' edge
		std::cout << "checking edge i = " << i << " out of " << E.rows() << std::endl;
		if ( abs(V_ren(E(i,0),0)-V_ren(E(i,1),0)) < eps ) {
			rendered_wireframe_edges.push_back(std::pair<int,int>(E(i,0),E(i,1)));	
		} else if ( abs(V_ren(E(i,0),1)-V_ren(E(i,1),1)) < eps ) {
			rendered_wireframe_edges.push_back(std::pair<int,int>(E(i,0),E(i,1)));
		}
	}
}

bool Dog::is_rectangular() {
	return (edgeStitching.submesh_to_edge_pt.size() == 0);
}

void Dog::setup_boundary_curves_indices() {
	double min_x = flatV.col(0).minCoeff(); double max_x = flatV.col(0).maxCoeff();
	double min_y = flatV.col(1).minCoeff(); double max_y = flatV.col(1).maxCoeff();
	//std::cout << "min_x = " << min_x << " max_x = " << max_x << " min_y = " << min_y << " max_y = " << max_y << std::endl;
	int left_lower_i = find_v_idx(flatV,Eigen::RowVector3d(min_x,min_y,0));
	int right_lower_i = find_v_idx(flatV,Eigen::RowVector3d(max_x,min_y,0));
	int left_upper_i = find_v_idx(flatV,Eigen::RowVector3d(min_x,max_y,0));
	//std::cout << "At indices " << left_lower_i << "," << right_lower_i << "," << left_upper_i << std::endl;

	//left_bnd,right_bnd,lower_bnd,upper_bnd;
	get_all_curves_on_parameter_line(left_lower_i, Eigen::RowVector3d(0,1,0), left_bnd);
	get_all_curves_on_parameter_line(right_lower_i, Eigen::RowVector3d(0,1,0), right_bnd);
	get_all_curves_on_parameter_line(left_lower_i, Eigen::RowVector3d(1,0,0), lower_bnd);
	get_all_curves_on_parameter_line(left_upper_i, Eigen::RowVector3d(1,0,0), upper_bnd);
	/*
	cout << "left ";for (auto idx: left_bnd ) cout << idx <<","; cout << endl;
	cout << "right "; for (auto idx: right_bnd ) cout << idx <<","; cout << endl;
	cout << "down "; for (auto idx: lower_bnd ) cout << idx <<","; cout << endl;
	cout << "up "; for (auto idx: upper_bnd ) cout << idx <<","; cout << endl;
	*/
}

int Dog::find_v_idx(Eigen::MatrixXd& Vertices, Eigen::RowVector3d v) {
	double eps = 1e-12;
	for (int i = 0; i < Vertices.rows(); i++) {
		if ((Vertices.row(i)-v).norm() < eps) return i;
	}
	return -1;
}

int Dog::find_other_v_idx(Eigen::MatrixXd& Vertices, int other_v_i, Eigen::RowVector3d v) {
	double eps = 1e-12;
	for (int i = 0; i < Vertices.rows(); i++) {
		if (i == other_v_i) continue;
		if ((Vertices.row(i)-v).norm() < eps) return i;
	}
	return -1;
}

void Dog::get_all_curves_on_parameter_line(int v_idx, const Eigen::RowVector3d& direction, std::vector<int>& indices) {
	double eps = 1e-10;
	int prev_idx = -1; int cur_idx = v_idx;
	indices.push_back(cur_idx);
	bool more_to_go = true;
	while (more_to_go) {
		more_to_go = false;
		for (int i = 0; i < quadTop.A[cur_idx].size(); i++) {
			Eigen::RowVector3d edge_dir = V.row(quadTop.A[cur_idx][i])-V.row(cur_idx);
			if (edge_dir.dot(direction) > eps) {
				prev_idx = cur_idx;
				cur_idx = quadTop.A[cur_idx][i];
				indices.push_back(cur_idx);
				more_to_go = true; 
				continue;
			}
		}
		if (!more_to_go) {
			// Found the end of the curve, maybe we can go backwards and get the same vertex from another connected component
			// Find another vertex with the same flattened location, but another component
			int other_v_i = find_other_v_idx(flatV, prev_idx, flatV.row(prev_idx));
			if (other_v_i != -1) {
				prev_idx = cur_idx;
				cur_idx = other_v_i;
				indices.push_back(cur_idx);
				more_to_go = true;
			}
		}
	}
	
}

void Dog::setup_uv_and_texture() {
	std::cout << "here" << std::endl;

	uv.resize(V.rows(),2); uv.col(0) = V.col(0); uv.col(1) = V.col(1);
	auto mesh_bb_size = uv.colwise().maxCoeff()-uv.colwise().maxCoeff();
	double mesh_W = mesh_bb_size[0], mesh_H = mesh_bb_size[1];

	// now go component by component and set the uv values accordingly

	// First set every vertex in the uv as its location + some offset of the x value of the mesh
	int submesh_n = get_submesh_n(); int subm_v_start = 0;
	for (int subm_i = 0; subm_i < submesh_n; subm_i++) {
		int subm_size = submeshVSize[subm_i];
		double x_offset = subm_i*mesh_W + 1;
		for (int v_idx = subm_v_start; v_idx < subm_v_start + subm_size; v_idx++) {
			uv.row(v_idx) << x_offset + V(v_idx,0), V(v_idx,1);
		}
		subm_v_start += subm_size;
	}

	// Now create a texture. Set the resolution for 1024Xsubmesh_n*(1024+1)
	// So height is 1024, and width is the paccking of everything
	int y_resolution = 1024, x_resolution = (1+1024)*submesh_n;
	// Set everything to be invisible, and then mark those that are inside
	text_R.resize(y_resolution,x_resolution);text_R.setZero();
  	text_G.resize(y_resolution,x_resolution);text_G.setZero();
  	text_B.resize(y_resolution,x_resolution);text_B.setZero();
  	text_A.resize(y_resolution,x_resolution);text_A.setZero();

  	// This should zero out half of the y things
  	for (int i = 0; i < y_resolution/2; i++) {
  		for (int j = 0; j < x_resolution; j++) {
  			text_A(i,j) = 255;
  		}
  	}

  	
  	// Now scale the uv. Scale it such that the biggest axis of the texture will be from 0 to 1
	// while the other one be at that ratio
	// move the minimal coordinates x to be 0, and same for y
	Eigen::VectorXd max_c = uv.colwise().maxCoeff();
	Eigen::VectorXd min_c = uv.colwise().minCoeff();
  	// move the minimal coordinates x to be 0, and same for y
  	double t_x = -1*min_c[0]; double t_y = -1*min_c[1];
  	//cout << "t_x = " << t_x << " t_y = " << t_y << endl;
  	// scale it such that the maximum 'x' distance will be 1, and same for y
  	double x_rad = max_c[0]-min_c[0]; double y_rad = max_c[1]-min_c[1];

  	double x_s,y_s;
	if (x_resolution >= y_resolution){
		x_s = 1./x_rad, y_s = (double(y_resolution)/x_resolution)*1./y_rad;
	} else {
		y_s = 1./y_rad, x_s = (double(x_resolution)/y_resolution)*1./x_rad;
	}
  	// For now scale it by the smaller factor so it will always fit to a 1x1 box
  	double scale = min(1./x_rad,1./y_rad);
  	for (int i = 0; i < uv.rows(); i++) {
    	uv.row(i) << x_s*(V(i,0)+t_x),y_s*(V(i,1)+t_y);
  	}
}

const Eigen::MatrixXd& Dog::getTexture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& text_Ri,
					Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& text_Gi,
					Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& text_Bi,
					Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& text_Ai) const {
	text_Ri = text_R; text_Gi = text_G; text_Bi = text_B;text_Ai = text_A;
}