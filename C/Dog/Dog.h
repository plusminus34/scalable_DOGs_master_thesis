#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <Eigen/Dense>
#include <vector>

#include "igl/serialize.h"

#include "../QuadMesh/Quad.h"

struct DogEdgeStitching  : public igl::Serializable {
	std::vector<Edge> edge_const_1, edge_const_2;
	std::vector<double> edge_coordinates;

	// The folds polylines
	std::vector<std::vector<EdgePoint>> stitched_curves;

	// Use for cases when it's important to have a precise representation (usually it doesn't)
	std::vector<CGAL::Exact_predicates_exact_constructions_kernel::FT> edge_coordinates_precise;

	// For now we don't serialize the exact CGAL coordinates
	void InitSerialization() {
      Add(edge_const_1,std::string("edge_const_1"));
      Add(edge_const_2,std::string("edge_const_2"));
      Add(edge_coordinates,std::string("edge_coordinates"));
      Add(stitched_curves,std::string("stitched_curves"));
    }
};

// Encapsulate a dog mesh, the stitching constraints for multiple components and its rendered mesh
class Dog : public igl::Serializable {
public:
	Dog(Eigen::MatrixXd V, Eigen::MatrixXi F, DogEdgeStitching edgeStitching, Eigen::MatrixXd V_ren, Eigen::MatrixXi F_ren,
		std::vector<int> submeshVSize);
	Dog(const Dog& dog);
	Dog(Eigen::MatrixXd V, Eigen::MatrixXi F);
	Dog(){/*Needed for deserilization*/};

	//void get_rendering_mesh(Eigen::MatrixXd& Vi, Eigen::MatrixXi& Fi) {Vi = V_ren; Fi = F_ren;}
	//void get_rendering_mesh(Eigen::MatrixXd& Vi) {Vi = V_ren;}
	
	void update_V(const Eigen::MatrixXd& V_new) {V = V_new; update_rendering_v();}
	void update_V_vector(const Eigen::VectorXd& x) {vec_to_mat2(x,V); update_rendering_v();}

	void get_submesh_min_max_i(int submesh_i, int& submesh_min_i, int& submesh_max_i);
	int get_submesh_n() const { return submesh_min_max_i.size();}
	int get_submesh_i_size(int submesh_i) const {return 1+submesh_min_max_i[submesh_i].second-submesh_min_max_i[submesh_i].first;}
	int v_to_submesh_idx(int v) const {return vi_to_submesh[v];}

	bool has_creases() const {return (edgeStitching.edge_const_1.size()>0);}

	const DogEdgeStitching& getEdgeStitching() const {return edgeStitching;}
	const Eigen::MatrixXi& getF() const {return F;}
	const Eigen::MatrixXd& getV() const {return V;}
	Eigen::MatrixXd& getVMutable() {return V;}
	Eigen::VectorXd getV_vector() const {Eigen::VectorXd x; mat2_to_vec(V,x); return x;}
	const Eigen::MatrixXi& getFrendering() const {return F_ren;}
	const Eigen::MatrixXd& getVrendering() const {return V_ren;}

	static void V_ren_from_V_and_const(const Eigen::MatrixXd& V, const DogEdgeStitching& edgeStitching, Eigen::MatrixXd& V_ren);

	void InitSerialization() {
      Add(V,std::string("_V"));
      Add(F,std::string("_F"));
      Add(V_ren,std::string("_V_ren"));
      Add(F_ren,std::string("_F_ren"));
      Add(edgeStitching,std::string("_edgeStitching"));
    }
	
private:
	void update_rendering_v();

	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	// The initial rendered (triangular) mesh
	Eigen::MatrixXd V_ren; Eigen::MatrixXi F_ren;

	// Edge stitching along multiple connected components in the DOG. Used to represent a piecewise developable mesh and in particular allow for folds.
	DogEdgeStitching edgeStitching;
	std::vector<std::pair<int,int>> submesh_min_max_i;
	std::vector<int> vi_to_submesh;
};


