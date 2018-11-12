#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <Eigen/Dense>
#include <vector>

#include "../QuadMesh/Quad.h"

struct DogEdgeStitching {
	std::vector<Edge> edge_const_1, edge_const_2;
	std::vector<double> edge_coordinates;
	// Use for cases when it's important to have a precise representation (usually it doesn't)
	std::vector<CGAL::Exact_predicates_exact_constructions_kernel::FT> edge_coordinates_precise;
};

class Dog {
public:
	Dog(Eigen::MatrixXd V, Eigen::MatrixXi F, DogEdgeStitching edgeStitching, Eigen::MatrixXd V_ren, Eigen::MatrixXi F_ren);
	Dog(const Dog& dog);

	void get_rendering_mesh(Eigen::MatrixXd& Vi, Eigen::MatrixXi& Fi) {Vi = V_ren; Fi = F_ren;}
	void get_rendering_mesh(Eigen::MatrixXd& Vi) {Vi = V_ren;}
	
	void update_rendering_v();

	static void V_ren_from_V_and_const(const Eigen::MatrixXd& V, const DogEdgeStitching& edgeStitching, Eigen::MatrixXd& V_ren);
	
private:
	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	// The initial rendered (triangular) mesh
	Eigen::MatrixXd V_ren; Eigen::MatrixXi F_ren;

	// Edge stitching along multiple connected components in the DOG. Used to represent a piecewise developable mesh and in particular allow for folds.
	DogEdgeStitching edgeStitching;
};