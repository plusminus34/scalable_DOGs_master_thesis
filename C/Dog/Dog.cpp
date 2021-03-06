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
	setup_uv();
	cout << "DOG setup complete" << endl;
}

Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F),flatV(V), V_ren(V), F_ren(Fsqr_to_F(F)) {
	submeshVSize.push_back(V.rows()); submeshFSize.push_back(F.rows());
	vi_to_submesh.assign(V.rows(),0);
	quad_topology(V,F,quadTop);
	setup_rendered_wireframe_edges_from_planar();
}

Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, DogEdgeStitching edgeStitching) :
				V(V), F(F), flatV(V), edgeStitching(edgeStitching), V_ren(V),  F_ren(Fsqr_to_F(F)) {
	quad_topology(V,F,quadTop);
	//find connected components / assign vertices to submeshes
	int num_submeshes = 0;
	vi_to_submesh.resize(V.rows(), -1);
	for(int i=0; i<V.rows(); ++i){
		if(vi_to_submesh[i] > -1) continue;
		int submesh = num_submeshes++;
		vector<int> tocheck(1,i);
		while(tocheck.size()>0){
			int v = tocheck[0];
			tocheck[0] = tocheck[tocheck.size()-1];
			tocheck.pop_back();
			vi_to_submesh[v] = submesh;
			for(int j=0; j<quadTop.VF[v].size(); ++j){
				int face = quadTop.VF[v][j];
				for(int k=0; k<4; ++k){
					int w = F(face, k);
					if(vi_to_submesh[w]==-1){
						vi_to_submesh[w]=-2;
						tocheck.push_back(w);
					}
				}
			}
		}
	}
	// set submesh sizes
	submeshVSize.resize(num_submeshes, 0);
	submeshFSize.resize(num_submeshes, 0);
	submesh_adjacency.resize(num_submeshes);
	for(int i=0; i<V.rows(); ++i) ++submeshVSize[vi_to_submesh[i]];
	for(int i=0; i<F.rows(); ++i) ++submeshFSize[vi_to_submesh[F(i,0)]];
	// and submesh_adjacency
	Eigen::MatrixXi submesh_A(num_submeshes, num_submeshes); submesh_A.setConstant(0);
	for(int i=0; i<edgeStitching.edge_const_1.size(); ++i){
		int submesh_1 = vi_to_submesh[edgeStitching.edge_const_1[i].v1];
		int submesh_2 = vi_to_submesh[edgeStitching.edge_const_2[i].v1];
		submesh_A(submesh_1, submesh_2) = 1;
	}
	for(int i=0;i<num_submeshes;++i){
		for(int j=i+1;j<num_submeshes;++j){
			if(submesh_A(i,j)==1){
				submesh_adjacency[i].push_back(j);
				submesh_adjacency[j].push_back(i);
			}
		}
	}
	// rest like in full constructor
	setup_stitched_curves_initial_l_angles_length();
	setup_rendered_wireframe_edges_from_planar();
	setup_uv();
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
	if (submesh_i == -1) return this;
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

Dog* Dog::get_curve_submesh(int curve_i) {
	if (curve_i < 0 || curve_i >= edgeStitching.stitched_curves.size()) return NULL;
	// find all vertices/faces along curve
	vector<int> v_to_local(V.rows(), -1);
	vector<int> f_to_local(F.rows(), -1);
	vector<EdgePoint> curve = edgeStitching.stitched_curves[curve_i];
	for(int i=0; i<curve.size(); ++i){
		int where = edgeStitching.edge_to_duplicates.at(curve[i].edge);
		Edge e1 = edgeStitching.edge_const_1[where];
		Edge e2 = edgeStitching.edge_const_2[where];
		for(int j=0; j<quadTop.VF[e1.v1].size(); ++j){
			for(int k=0; k<quadTop.VF[e1.v2].size(); ++k){
				if(quadTop.VF[e1.v1][j] == quadTop.VF[e1.v2][k]){
					f_to_local[quadTop.VF[e1.v1][j]] = -2;
				}
			}
		}
		for(int j=0; j<quadTop.VF[e2.v1].size(); ++j){
			for(int k=0; k<quadTop.VF[e2.v2].size(); ++k){
				if(quadTop.VF[e2.v1][j] == quadTop.VF[e2.v2][k]){
					f_to_local[quadTop.VF[e2.v1][j]] = -2;
				}
			}
		}
	}
	for(int i=0;i<F.rows();++i){
		if(f_to_local[i] == -2)
			for(int j=0;j<4;++j) v_to_local[F(i,j)] = -2;
	}
	//fill them into curve_V and curve_F
	int count = 0;
	for(int i=0; i<V.rows(); ++i){
		if(v_to_local[i] == -2) v_to_local[i] = count++;
	}
	Eigen::MatrixXd curve_V(count,3);
	count = 0;
	for(int i=0; i<F.rows(); ++i){
		if(f_to_local[i] == -2) f_to_local[i] = count++;
	}
	Eigen::MatrixXi curve_F(count,4);
	for(int i=0; i<V.rows(); ++i){
		if(v_to_local[i] > -1) curve_V.row(v_to_local[i]) = V.row(i);
	}
	for(int i=0; i<F.rows(); ++i){
		if(f_to_local[i] > -1) {
			for(int j=0; j<4; ++j) curve_F(f_to_local[i], j) = v_to_local[F(i,j)];
		}
	}

	// curve_eS is basically edgeStitching, but remove everything but curve i
	DogEdgeStitching curve_eS;
	curve_eS.stitched_curves.push_back(curve);
	curve_eS.multiplied_edges_start.resize(curve.size());
	curve_eS.multiplied_edges_num.resize(curve.size(), 1);
	for(int i=0; i<curve.size(); ++i){
		int where = edgeStitching.edge_to_duplicates.at(curve[i].edge);
		Edge e1 = edgeStitching.edge_const_1[where];
		Edge e2 = edgeStitching.edge_const_2[where];
		Edge local_e1(v_to_local[e1.v1], v_to_local[e1.v2]);
		Edge local_e2(v_to_local[e2.v1], v_to_local[e2.v2]);
		curve_eS.edge_to_duplicates[local_e1] = i;
		curve_eS.edge_to_duplicates[local_e2] = i;
		curve_eS.edge_const_1.push_back(local_e1);
		curve_eS.edge_const_2.push_back(local_e2);
		curve_eS.edge_coordinates.push_back(curve[i].t);
		curve_eS.multiplied_edges_start[i] = i;
	}
	//missing other less important stuff
	Dog* submeshDog = new Dog(curve_V, curve_F, curve_eS);
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
	//V_ren = V;
	//V_ren_from_V_and_const(V,edgeStitching,V_ren);
	update_Vren();
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
	V_ren *= render_scale;
}

void Dog::get_submesh_min_max_i(int submesh_i, int& submesh_min_i, int& submesh_max_i, bool vertices) const {
	std::vector<int> idx_list;
	if (vertices)
		idx_list = submeshVSize;
	else
		idx_list = submeshFSize;

	submesh_min_i = 0; int submesh_vn = idx_list[0];
	for (int sub_i = 1; sub_i <= submesh_i; sub_i++) {
		submesh_min_i += submesh_vn;
		submesh_vn = idx_list[sub_i];
	}
	submesh_max_i = submesh_min_i+submesh_vn-1;
}

int Dog::v_in_submesh(int v) const {
	int submesh_i = vi_to_submesh[v];
	for(int i=1; i<=submesh_i; ++i) v -= submeshVSize[i-1];
	return v;
}

int Dog::v_ren_idx_to_v_idx(int v_idx) const {
	double eps = 1e-5;
	auto coords = V_ren.row(v_idx);
	for (int i = 0; i < V.rows(); i++) {
		if ((V.row(i)/render_scale-coords).norm() < eps)
			return i;
	}
	// did not find matching point
	std::cout << "could not find point!" << std::endl;
	return -1;
}

int Dog::v_ren_idx_to_edge(int v_idx, EdgePoint& edgePt) const {
	double eps = 1e-10;
	auto coords = V_ren.row(v_idx) / render_scale;
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
		//std::cout << "checking edge i = " << i << " out of " << E.rows() << std::endl;
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

void Dog::setup_uv() {
	uv.resize(V.rows(),2); uv.col(0) = V.col(0); uv.col(1) = V.col(1);
	Eigen::VectorXd max_c = uv.colwise().maxCoeff();
	Eigen::VectorXd min_c = uv.colwise().minCoeff();
  	// move the minimal coordinates x to be 0, and same for y
  	double t_x = -1*min_c[0]; double t_y = -1*min_c[1];
  	//cout << "t_x = " << t_x << " t_y = " << t_y << endl;
  	// scale it such that the maximum 'x' distance will be 1, and same for y
  	double x_rad = max_c[0]-min_c[0]; double y_rad = max_c[1]-min_c[1];
  	// For now scale it by the smaller factor so it will always fit to a 1x1 box
  	double scale = min(1./x_rad,1./y_rad);
  	for (int i = 0; i < uv.rows(); i++) {
    	uv.row(i) << scale*(V(i,0)+t_x),scale*(V(i,1)+t_y);
  	}
}
