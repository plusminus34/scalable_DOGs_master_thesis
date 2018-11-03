#include "PlanarArrangement.h"

#include <igl/combine.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/jet.h>
#include <igl/triangle/triangulate.h>

#include <boost/range/irange.hpp>

void PlanarArrangement::add_segments(const std::vector<Segment_2>& segments) {
	insert(arr, segments.begin(), segments.end());
}

void PlanarArrangement::add_segment(const Segment_2& segment) {
	add_segments({segment});
}

void PlanarArrangement::add_polylines(const std::vector<Polyline_2>& polylines) {
	// we need to insert the whole polyline together since we are using an arrangement with history
	insert(arr, polylines.begin(), polylines.end());
	/*
	//for (auto p: polylines) { std::cout << "p = " << p << std::endl << "ok?" << std::endl;}
	std::vector<Segment_2> segments;
	for (auto polyline : polylines) {
		for (auto it = polyline.subcurves_begin(); it != polyline.subcurves_end(); it++) {
			segments.push_back(*it);
		}	
	}
	add_segments(segments);
	*/
}

void PlanarArrangement::add_polyline(const Polyline_2& polyline) {
	insert(arr, polyline);
}

bool PlanarArrangement::locate_point_on_vertex(const Point_2& pt, Vertex_const_handle& v) {
	Point_location pl(arr);
	Vertex_const_handle* v_ptr;
	typename CGAL::Arr_point_location_result<Arrangement_2>::Type obj = pl.locate(pt);
	// Check if it is on a vertex
	if ( (v_ptr = boost::get<Vertex_const_handle>(&obj) ) ) {
        v = *v_ptr;
		return true;
	} else {
		return false;
	}
}

void PlanarArrangement::get_visualization_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors) {
	std::vector<Eigen::MatrixXd> V_list; std::vector<Eigen::MatrixXi> F_list;

	// Build faces polygons
	Arrangement_2::Face_const_iterator fit;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (fit->is_unbounded()) continue; // skip the unbounded face

		Eigen::MatrixXd Vk; Eigen::MatrixXi Fk;
    	get_face_vertices(fit, Vk);
    	auto range = boost::copy_range<std::vector<int>>(boost::irange(0, int(Vk.rows())));
    	std::vector<std::vector<int> > vlist_for_tri; vlist_for_tri.push_back(range);
    	//igl::polygon_mesh_to_triangle_mesh(vlist_for_tri,Fk);
    	Eigen::MatrixXi Ek(Vk.rows(),2); for (int i = 0; i < Vk.rows(); i++) {Ek.row(i) << i,(i+1) % Vk.rows();}
    	Eigen::MatrixXi H;
    	Eigen::MatrixXd newVk;
    	igl::triangle::triangulate(Vk,Ek, H, "p",newVk,Fk);
    	Eigen::MatrixXd V3d(newVk.rows(),3); V3d.setZero(); V3d.col(0) = Vk.col(0); V3d.col(1) = Vk.col(1);
    	V_list.push_back(V3d);
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

int PlanarArrangement::get_faces_n() {
	return arr.number_of_faces();
}

int PlanarArrangement::get_vertices_n() {
	return arr.number_of_vertices();
}

void PlanarArrangement::get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p) {
	typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
	typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
	// switch edge direction if not ordered as the segments
	//if (curr->curve().subcurves_begin()->source()!= curr->source()->point()) {curr = curr->twin();}

	// count number of vertices
	int v_num = 0;
	do {
		// count also polyline edges (that are not necessarily ones in the graph)
		//v_num += curr->curve().subcurves_end()-curr->curve().subcurves_begin();
		for (auto it = curr->curve().subcurves_begin(); it != curr->curve().subcurves_end(); it++) {
			v_num++;
		}
		curr++;
	} while (curr != circ);
	//std::cout << "vnum = "  << v_num << std::endl;
	// Fill up p with the vertices
	//p.resize(v_num,3);
	p.resize(v_num,2);


	int ri = 0; curr = circ;
	do {
		//p.row(ri) << CGAL::to_double(curr->source()->point().x()),CGAL::to_double(curr->source()->point().y()),0;
		// Go through segments in the polyline (not every vertex in the polyline is an actual vertex in the graph)
		//std::cout << "vertex starting at " << curr->source()->point() << std::endl;
		bool flipped_order = false;
		if (curr->curve().subcurves_begin()->source()!= curr->source()->point()) {flipped_order = true;}
		for (auto it = curr->curve().subcurves_begin(); it != curr->curve().subcurves_end(); it++) {
			Point_2 next_pt = flipped_order ? it->target() : it->source();
			//std::cout << "Segment: " << *it << std::endl << " and adding " << next_pt << std::endl;
			p.row(ri) << CGAL::to_double(next_pt.x()),CGAL::to_double(next_pt.y());
			ri++;
		}
		curr++;
	} while (curr != circ);
	//std::cout << "p = " << p << std::endl;
	//exit(1);
}

void PlanarArrangement::get_boundary_polyline_pts(std::vector<Point_2>& pts) {
	pts.clear();
	// Build faces polygons
	Arrangement_2::Face_const_iterator fit;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (fit->is_unbounded()) {
			Eigen::MatrixXd p;
			std::cout << "found unbounded face!" << std::endl;
			get_face_vertices(fit,p);
			std::cout << "p = " << p << std::endl;
			pts.resize(p.rows());
			for (int i = 0; i < p.rows(); i++) pts[i] = Point_2(p(i,0),p(i,1));
		}
	}
}

void get_multiple_arrangements_visualization_mesh(std::vector<PlanarArrangement*> arrangements, double spacing,
							Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& colors) {
	// Visualize all
	std::vector<Eigen::MatrixXd> V_list; std::vector<Eigen::MatrixXi> F_list; std::vector<Eigen::MatrixXd> F_colors_list;
	int cnt = 0;
	for (auto arr: arrangements) {
		Eigen::MatrixXd Vk,Ck; Eigen::MatrixXi Fk;
		arr->get_visualization_mesh(Vk,Fk,Ck);
		Vk.rowwise() += Eigen::RowVector3d(cnt*spacing,0,0);
		V_list.push_back(Vk);
		F_list.push_back(Fk);
		F_colors_list.push_back(Ck);
		cnt++;
	}
	igl::combine(V_list,F_list, V, F);
	colors.resize(F.rows(),3); int kv = 0;
	for (auto Ck : F_colors_list) {
		int ni = Ck.rows();
		colors.block(kv,0,ni,3) = Ck;
		kv+=ni;
	}
}