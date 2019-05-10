#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <Eigen/Dense>
#include <vector>

#include "igl/serialize.h"

#include "../QuadMesh/Quad.h"

struct DogEdgeStitching  : public igl::Serializable {
	// Stores edge constraints (if an edge is duplicated n times, than this edge involves n-1 linearly independent constraints)
	// The duplicated constraints always needs to be stored consecuitvely in this list
	std::vector<Edge> edge_const_1, edge_const_2;
	std::vector<double> edge_coordinates;

	// A list of all the duplicated edge point constraints size (usually there are 2 edge points, but a vertex can have more)
	// The start is their starting offset in std:vector<Edge> edge_const_1/2 and the num is the number of duplicated constraints
	// If there are n duplicated edges, then there are n-1 duplicated edge constraints
	// Each one of the edges will map to multplied edges start and the num will be n-1
	std::vector<int> multiplied_edges_start; std::vector<int> multiplied_edges_num;
	std::map<Edge, int> edge_to_duplicates; // A map between an edge and its index in multiplied_edges_start

	// The folds polylines (holds an arbitrary edge point, and not all the duplicated)
	// Can use any EdgePoint to map to other equal edges with edge_to_duplicates to get an index in multiplied_edges_start
	std::vector<std::vector<EdgePoint>> stitched_curves;

	// Use for cases when it's important to have a precise representation (usually it doesn't)
	std::vector<CGAL::Exact_predicates_exact_constructions_kernel::FT> edge_coordinates_precise;

	// For each submesh, hold the indices to edge-points (crease points) that are on that submesh. The indices point to edge_const_1 and edge_coordinates vectors
	std::vector<std::vector<int>> submesh_to_edge_pt;
	std::vector<std::vector<EdgePoint>> submesh_to_bnd_edge;

	int get_vertex_edge_point_deg(Edge& edge) const;

	// For now we don't serialize the exact CGAL coordinates
	void InitSerialization() {
      Add(edge_const_1,std::string("edge_const_1"));
      Add(edge_const_2,std::string("edge_const_2"));
      Add(edge_coordinates,std::string("edge_coordinates"));
      Add(multiplied_edges_start,std::string("multiplied_edges_start"));
      Add(multiplied_edges_num,std::string("multiplied_edges_num"));
      Add(edge_to_duplicates,std::string("edge_to_duplicates"));
      Add(stitched_curves,std::string("stitched_curves"));
      Add(submesh_to_edge_pt,std::string("submesh_to_edge_pt"));
      Add(submesh_to_bnd_edge,std::string("submesh_to_bnd_edge"));
    }
};

// Encapsulate a dog mesh, the stitching constraints for multiple components and its rendered mesh
class Dog : public igl::Serializable {
public:
	Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, DogEdgeStitching edgeStitching, const Eigen::MatrixXd& V_ren, const Eigen::MatrixXi& F_ren,
		std::vector<int> submeshVSize, std::vector<int> submeshFSize, const std::vector< std::vector<int> >& submesh_adjacency);
	Dog(const Dog& dog);
	Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	Dog(){/*Needed for deserilization*/};

	void update_Vren();
	void update_V(const Eigen::MatrixXd& V_new) {V = V_new; update_rendering_v();}
	void update_V_vector(const Eigen::VectorXd& x) {vec_to_mat2(x,V); update_rendering_v();}
	void update_submesh_V(int submesh_i, const Eigen::MatrixXd& submeshV);
	std::vector<std::vector <int> > get_submesh_adjacency() {return submesh_adjacency;}

	Dog* get_submesh(int submesh_i);

	int get_v_num() {return V.rows();}

	// If true, gives the min/max in terms of submesh vertices, otherwise faces
	void get_submesh_min_max_i(int submesh_i, int& submesh_min_i, int& submesh_max_i, bool vertices = true);
	int get_submesh_n() const { return submeshVSize.size();}
	int get_submesh_i_size(int submesh_i) const {return submeshVSize[submesh_i];}
	int v_to_submesh_idx(int v) const {return vi_to_submesh[v];}

	bool has_creases() const {return (edgeStitching.edge_const_1.size()>0);}

	bool is_crease_vertex_flat(int curve_i, int edge_i) const;

	const DogEdgeStitching& getEdgeStitching() const {return edgeStitching;}
	const Eigen::MatrixXi& getF() const {return F;}
	const Eigen::MatrixXd& getV() const {return V;}
	//const Eigen::MatrixXd& getFlatV() const {return flatV;} somehow nto working now
	const QuadTopology& getQuadTopology() const {return quadTop;}
	Eigen::MatrixXd& getVMutable() {return V;}
	Eigen::VectorXd getV_vector() const {Eigen::VectorXd x; mat2_to_vec(V,x); return x;}
	Eigen::MatrixXi getFTriangular() const {return Fsqr_to_F(F);} // useful for the editor who needs a triangular mesh (still different then the rendering)
	const Eigen::MatrixXi& getFrendering() const {return F_ren;}
	const Eigen::MatrixXd& getVrendering() const {return V_ren;}

	void get_2_submeshes_vertices_from_edge(const Edge& edge, int &v1_out, int &v2_out, int &w1_out, int& w2_out) const;
	void get_2_inner_vertices_from_edge(const Edge& edge, int &v1_out, int &v2_out) const;

	// Helper functinon to go from V to V_ren (rendering) vice-versa (used for the GUI and editor vertex picking)
	static void V_ren_from_V_and_const(const Eigen::MatrixXd& V, const DogEdgeStitching& edgeStitching, Eigen::MatrixXd& V_ren);
	int v_ren_idx_to_v_idx(int v_idx) const;
	int v_ren_idx_to_edge(int v_idx, EdgePoint& edgePt) const;

	bool is_rectangular();
	void setup_rendered_wireframe_edges_from_planar();
	const std::vector< std::pair<int,int> >& get_rendered_wireframe_edges() const {return rendered_wireframe_edges;}
	void setup_boundary_curves_indices();

	void InitSerialization() {
      Add(V,std::string("_V"));
      Add(F,std::string("_F"));
      Add(F,std::string("_flatV"));
      Add(left_bnd,std::string("_left_bnd"));
      Add(right_bnd,std::string("_right_bnd"));
      Add(lower_bnd,std::string("_lower_bnd"));
      Add(upper_bnd,std::string("_upper_bnd"));
      Add(quadTop,std::string("_quadTop"));
      Add(V_ren,std::string("_V_ren"));
      Add(F_ren,std::string("_F_ren"));
      Add(rendered_wireframe_edges,std::string("_rendered_wireframe_edgesn"));
      Add(submeshVSize,std::string("_submeshVSize"));
      Add(submeshFSize,std::string("_submeshFSize"));
      Add(submesh_adjacency,std::string("_submesh_adjacency"));
      Add(edgeStitching,std::string("_edgeStitching"));
      Add(stitched_curves_l,std::string("stitched_curves_l"));
      Add(stitched_curves_angles,std::string("stitched_curves_angles"));
      Add(stitched_curves_curvature,std::string("stitched_curves_curvature"));
      Add(vi_to_submesh,std::string("_vi_to_submesh"));
    }
	
	// The initial length/angles/curvatures of the initial stitched curves
	std::vector<std::vector<double>> stitched_curves_l; std::vector<std::vector<double>> stitched_curves_angles; std::vector<std::vector<double>> stitched_curves_curvature;
	std::vector<int> left_bnd,right_bnd,lower_bnd,upper_bnd;

private:
	void update_rendering_v();
	void setup_stitched_curves_initial_l_angles_length();
	void get_all_curves_on_parameter_line(int v_idx, const Eigen::RowVector3d& direction, std::vector<int>& indices);
	static int find_v_idx(Eigen::MatrixXd& Vertices, Eigen::RowVector3d v);
	static int find_other_v_idx(Eigen::MatrixXd& Vertices, int other_v_i, Eigen::RowVector3d v);

	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	Eigen::MatrixXd flatV;
	// Indices of boundary curves (also when there are creases). Only relevant for square patches and used for wallpapers.
	QuadTopology quadTop;
	// The initial rendered (triangular) mesh
	Eigen::MatrixXd V_ren; Eigen::MatrixXi F_ren;
	// The rendered quad wireframe
	std::vector< std::pair<int,int> > rendered_wireframe_edges;

	// Edge stitching along multiple connected components in the DOG. Used to represent a piecewise developable mesh and in particular allow for folds.
	std::vector<int> submeshVSize; std::vector<int> submeshFSize;
	std::vector< std::vector<int> > submesh_adjacency;
	DogEdgeStitching edgeStitching;
	std::vector<int> vi_to_submesh;
};
