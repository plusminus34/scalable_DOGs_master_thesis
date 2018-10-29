#include "PlanarArrangement.h"
#include "arr_print.h"

#include <igl/combine.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <boost/range/irange.hpp>

void PlanarArrangement::add_polylines(std::vector<Polyline_2>& polylines) {
	insert(arr, polylines.begin(), polylines.end());
}

void PlanarArrangement::add_polyline(Polyline_2& polyline) {
	insert(arr, polyline);
}

void PlanarArrangement::get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd& colors) {
	int fn = get_face_n();

	std::vector<Eigen::MatrixXd> V_list; std::vector<Eigen::MatrixXi> F_list;

	// Build faces polygons
	Arrangement_2::Face_const_iterator fit;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		Eigen::MatrixXd Vk; Eigen::MatrixXi Fk;
    	get_face_vertices(fit, Vk);
    	auto range = boost::copy_range<std::vector<int>>(boost::irange(0, int(Vk.rows())));
    	std::vector<std::vector<int> > vlist_for_tri; vlist_for_tri.push_back(range);
    	igl::polygon_mesh_to_triangle_mesh(vlist_for_tri,Fk);

    	V_list.push_back(Vk);
    	F_list.push_back(Fk);
	}
	igl::combine(V_list,F_list, V, F);

	// Color every polygonal face (now represented as many triangles for rendering) in another color
	Eigen::VectorXd components(F.rows());
	int c = 0;
	for (int i = 0; i < F_list.size();i++) {
		for (int j = 0; j < F_list[i].rows(); j++) {
			components[c] = i;
			c++;
		}
	}
	igl::jet(components,true,colors);
}

int PlanarArrangement::get_face_n() {
	return arr.number_of_faces();
}

int PlanarArrangement::get_vertices_n() {
	return arr.number_of_vertices();
}

void PlanarArrangement::get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p) {
	typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
	typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;

	// count number of vertices
	int v_num = 1;
	do {++curr; v_num++;} while (curr != circ);
	p.resize(v_num,3);

	curr = circ;
 	do {
		p.row(0) << CGAL::to_double(curr->source()->point().x()),CGAL::to_double(curr->source()->point().y()),0;
		++curr;
	} while (curr != circ);
}