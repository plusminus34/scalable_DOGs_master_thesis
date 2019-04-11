#include "PlanarArrangement.h"

#include <igl/combine.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/jet.h>
#include <igl/triangle/triangulate.h>

#include <boost/range/irange.hpp>
#include <random>

void PlanarArrangement::add_segments(const std::vector<Segment_2>& segments) {
	insert(arr, segments.begin(), segments.end());
}

void PlanarArrangement::add_segment(const Segment_2& segment) {
	add_segments({segment});
}

void PlanarArrangement::add_polylines(const std::vector<Polyline_2>& polylines) {
	// we need to insert the whole polyline together since we are using an arrangement with history
	insert(arr, polylines.begin(), polylines.end());
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
	// shuffle colors

	//std::random_device rd; std::mt19937 g(rd()); 
	//auto color_permute = boost::copy_range<std::vector<int>>(boost::irange(0, int(F_list.size())));
	//std::shuffle(color_permute.begin(), color_permute.end(), g);
	Eigen::VectorXd components(F.rows());
	int c = 0;
	for (int i = 0; i < F_list.size();i++) {
		for (int j = 0; j < F_list[i].rows(); j++) {
			components[c] = i;//color_permute[i];
			c++;
		}
	}
	if (components.maxCoeff() < 15) {
		igl::jet(components,true,colors);
	} else {
		colors.resize(F.rows(),3);
		// too many colors, better to set everything to white
		for (int i = 0; i < F.rows(); i++) colors.row(i) << 1,1,1;
	}
}

int PlanarArrangement::get_faces_n() {
	return arr.number_of_faces();
}

int PlanarArrangement::get_vertices_n() {
	return arr.number_of_vertices();
}

void PlanarArrangement::get_face_vertices(Arrangement_2::Face_const_handle f, Eigen::MatrixXd& p) const {
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
	p.resize(v_num,2);
	int ri = 0; curr = circ;
	do {
		//std::cout << "Edge from " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
		std::vector<Segment_2> polyline_segments(curr->curve().subcurves_begin(),curr->curve().subcurves_end());
		bool flipped_order = false;
		if (curr->curve().subcurves_begin()->source()!= curr->source()->point()) {flipped_order = true;}
		if (flipped_order) {
			std::reverse(std::begin(polyline_segments), std::end(polyline_segments));
		}
		for (auto seg: polyline_segments) {
			if (!flipped_order) {
				p.row(ri) << CGAL::to_double(seg.source().x()),CGAL::to_double(seg.source().y());
			} else {
				p.row(ri) << CGAL::to_double(seg.target().x()),CGAL::to_double(seg.target().y());
			}
			
			ri++;
		}
		curr++;
	} while (curr != circ);
}

void PlanarArrangement::get_face_vertices_from_circulator_iter(Arrangement_2::Ccb_halfedge_const_circulator circ, 
						std::vector<Point_2>& p) {
	//typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
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
	p.resize(v_num);
	int ri = 0; curr = circ;
	do {
		//std::cout << "Edge from " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
		std::vector<Segment_2> polyline_segments(curr->curve().subcurves_begin(),curr->curve().subcurves_end());
		bool flipped_order = false;
		if (curr->curve().subcurves_begin()->source()!= curr->source()->point()) {flipped_order = true;}
		if (flipped_order) {
			std::reverse(std::begin(polyline_segments), std::end(polyline_segments));
		}
		for (auto seg: polyline_segments) {
			if (!flipped_order) {
				p[ri] = seg.source();
			} else {
				p[ri] = seg.target();
			}
			
			ri++;
		}
		curr++;
	} while (curr != circ);
}

void PlanarArrangement::get_faces_pts(std::vector<std::vector<Point_2>>& pts) const {
	pts.clear();
	// Build faces polygons
	Arrangement_2::Face_const_iterator fit;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (!fit->is_unbounded()) {
			Eigen::MatrixXd p; get_face_vertices(fit,p);
			std::vector<Point_2> face_pts(p.rows());
			for (int i = 0; i < p.rows(); i++) face_pts[i] = Point_2(p(i,0),p(i,1));
			pts.push_back(face_pts);
		}
	}
}

void PlanarArrangement::get_faces_adjacency_list(std::vector<std::vector<int>>& A) const {
	A.clear();
	std::map<Arrangement_2::Face_const_iterator, int> face_to_id;
	// First number each face
	Arrangement_2::Face_const_iterator fit; int cnt = 0;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (!fit->is_unbounded()) {
			face_to_id[fit] = cnt; cnt++;
		}
	}
	A.resize(cnt);

	int face_idx = 0;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (!fit->is_unbounded()) {
			//face_to_id[fit] = cnt; cnt++;
			//std::cout << "face idx = " << face_idx << " and map value is " << face_to_id[fit] << std::endl;


			// Go through all of the vertices in the face, and add their neighbouring faces
			typename Arrangement_2::Ccb_halfedge_const_circulator circ = fit->outer_ccb();
			typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;

			/*
			// This includes faces sharing a vertex
			do {
				auto v_handle = curr->source();

				Arrangement_2::Halfedge_around_vertex_const_circulator first_v_edge, curr_v_edge;
				first_v_edge = curr_v_edge = v_handle->incident_halfedges();
				//std::cout << "The neighbors of the vertex (" << v->point() << ") are:";
				do {
					// Note that the current halfedge is directed from u to v:
					//int f1 = face_to_id[curr_v_edge->face()], f2 = face_to_id[curr_v_edge->twin()->face()];
					Arrangement_2::Face_const_iterator f1 = curr_v_edge->face(); Arrangement_2::Face_const_iterator f2 = curr_v_edge->twin()->face();

					if ((!f1->is_unbounded()) && (face_to_id[f1] != face_idx) ) A[face_idx].push_back(face_to_id[f1]);
					if ((!f2->is_unbounded()) && (face_to_id[f2] != face_idx) ) A[face_idx].push_back(face_to_id[f2]);
					//std::cout << "Connected to " << f1 << " and " << f2 << std::endl;
				} while (++curr_v_edge != first_v_edge);
				curr++;
			} while (curr != circ);
			*/

			do {
				auto nb_face = curr->twin()->face();
				if (!nb_face->is_unbounded()) {
					int nb = face_to_id[curr->twin()->face()];
					A[face_idx].push_back(nb);
				}
				curr++;
			} while (curr != circ);
			
			// remove duplicates
			std::sort(A[face_idx].begin(), A[face_idx].end());
    		auto last = std::unique(A[face_idx].begin(), A[face_idx].end());
    		A[face_idx].erase(last, A[face_idx].end()); 
			face_idx++;
		}
			
	}
}

void PlanarArrangement::get_visualization_edges(Eigen::MatrixXd& bnd_pts1, Eigen::MatrixXd& bnd_pts2) {
	/*
  std::vector<std::vector<Point_2>> all_faces_pts; 
  get_faces_pts(all_faces_pts);
  std::vector<Point_2> bnd_pts = all_faces_pts[1];
  bnd_pts1.resize(bnd_pts.size(),3); bnd_pts2.resize(bnd_pts.size(),3);
  
  for (int i = 0; i < bnd_pts.size();i++) {
    double x = CGAL::to_double(bnd_pts[i].x()), y = CGAL::to_double(bnd_pts[i].y());
    int next_i = (i+1)%bnd_pts.size();
    double next_x = CGAL::to_double(bnd_pts[next_i].x()), next_y = CGAL::to_double(bnd_pts[next_i].y());
    bnd_pts1.row(i) << x,y,0;
    bnd_pts2.row(i) << next_x,next_y,0;
  }
  */
}
void PlanarArrangement::get_faces_polygons(std::vector<Polygon_2>& polygons) const {
	std::vector<std::vector<Point_2>> facePts; get_faces_pts(facePts);
	polygons.resize(facePts.size());
	for (int i = 0; i < facePts.size();i++) {polygons[i] = Polygon_2(facePts[i].begin(), facePts[i].end());}
}

void PlanarArrangement::get_faces_polygons_with_holes(std::vector<Polygon_with_holes_2>& polygons) const {
	// Build faces polygons
	Arrangement_2::Face_const_iterator fit;
	for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
		if (!fit->is_unbounded()) {
			// Get outer face points
			std::vector<Point_2> face_pts; get_face_vertices_from_circulator_iter(fit->outer_ccb(), face_pts);
			//Polygon_with_holes_2 polygon(Polygon_2(face_pts.begin(), face_pts.end()), fit->holes_begin(),fit->holes_end());
			Polygon_2 outer_bnd(face_pts.begin(), face_pts.end());
			std::vector<Polygon_2> holes;	
			for (auto hole_fit  = fit->holes_begin(); hole_fit != fit->holes_end(); hole_fit++) {
				std::vector<Point_2> hole_pts; get_face_vertices_from_circulator_iter(*hole_fit, hole_pts);
				holes.push_back(Polygon_2(hole_pts.begin(), hole_pts.end()));
			}
			Polygon_with_holes_2 polygon_with_holes(outer_bnd,holes.begin(),holes.end());
			polygons.push_back(polygon_with_holes);
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

void get_multiple_arrangements_visualization_edges(std::vector<PlanarArrangement*> arrangements, double spacing, 
  Eigen::MatrixXd& bnd_pts1, Eigen::MatrixXd& bnd_pts2) {

  std::vector<Point_2> pts1,pts2;
  for (int i = 0; i < arrangements.size(); i++) {
  	std::vector<std::vector<Point_2>> faces_pts; arrangements[i]->get_faces_pts(faces_pts);
  	for (auto f_pts: faces_pts) {
  		// add spacing
  		for(auto& pt : f_pts) pt = Point_2(pt.x()+spacing*i,pt.y());
  		// add to list of edges
  		pts1.insert(pts1.end(),f_pts.begin(),f_pts.end());
  		pts2.insert(pts2.end(),f_pts.begin()+1,f_pts.end());pts2.push_back(f_pts[0]);
  	}
  }
  bnd_pts1.resize(pts1.size(),3); bnd_pts2.resize(pts2.size(),3);
  for (int i = 0; i < pts1.size();i++) {
    bnd_pts1.row(i) << CGAL::to_double(pts1[i].x()),CGAL::to_double(pts1[i].y()),0;
    bnd_pts2.row(i) << CGAL::to_double(pts2[i].x()),CGAL::to_double(pts2[i].y()),0;
  }
}
void sort_segments(const std::vector<Segment_2>& unsorted_seg, const Point_2& firstPoint,
						std::vector<Segment_2>& sorted_seg) {
	sorted_seg.resize(unsorted_seg.size());
	std::vector<Segment_2> available_segs = unsorted_seg;
	
	bool found_next = false;
	Point_2 next_pt = firstPoint;

	if (unsorted_seg.size()>10) {
		std::cout << "init! next point = " << next_pt << std::endl;	
	}
	
	for (auto seg: available_segs) {std::cout << "seg = " << seg << std::endl;}
	int seg_n = available_segs.size();
	for (int seg_i = 0; seg_i < seg_n; seg_i++) {
		for (int i = 0; i < available_segs.size(); i++) {
			if (available_segs[i].source() == next_pt) {
				sorted_seg[seg_i] = available_segs[i];
				found_next = true;
			}
			if (available_segs[i].target() == next_pt) {
				sorted_seg[seg_i] = Segment_2(available_segs[i].target(),available_segs[i].source());
				found_next = true;
			}
			if (found_next) {
				available_segs.erase(available_segs.begin() + i);
				next_pt = sorted_seg[seg_i].target();
				if (seg_i < 10) {
					std::cout << "found segment " << sorted_seg[seg_i] << std::endl;
					std::cout << "next pt = " << next_pt << std::endl;	
				}
				break;
			}
		}
	}
	



}