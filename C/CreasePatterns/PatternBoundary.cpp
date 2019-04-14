
#include "PatternBoundary.h"


PatternBoundary::PatternBoundary(const std::vector<Polygon_2>& boundary_polygons) {
	std::cout << std::endl << "Polygon: " << std::endl << boundary_polygons[0] << std::endl; //exit(1);
	std::vector<Number_type> areas; Number_type max_area(0); int boundary_poly_i = 0;
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (!boundary_polygons[i].is_simple()) {

			//std::cout << std::endl << std::endl << "Polygon:" << std::endl << boundary_polygons[i] << std::endl;
			std::cout << "Error at PatternBoundary: Polygon " << i << " is not simple" << std::endl;
			exit(1);
			
		}
		if (abs(boundary_polygons[i].area()) > max_area) {
			max_area = abs(boundary_polygons[i].area());
			boundary_poly_i = i;
		}
	}
	//std::cout << "max_area = " << max_area << std::endl; int wait; std::cin >> wait;
	outer_boundary = boundary_polygons[boundary_poly_i];
	for (int i = 0; i < boundary_polygons.size(); i++) {
		if (i!= boundary_poly_i) {holes.push_back(Polygon_2(boundary_polygons[i]));}
	}

	// Save all of the segments (needed for snapping vertices to edges)
	for (auto poly: boundary_polygons) {
		for (auto ei = poly.edges_begin(); ei != poly.edges_end(); ei++) {
			all_polygons_edges.push_back(*ei);
		}
	}
}

PatternBoundary::PatternBoundary(const std::vector<Polyline_2>& boundary_polylines) : 
	PatternBoundary(polylines_to_polygons(boundary_polylines)) {
	// empty on purpose
}

std::vector<Point_2> PatternBoundary::get_all_boundary_points() const {
	std::vector<Point_2> boundary_vertices(outer_boundary.vertices_begin(),outer_boundary.vertices_end());
	for (auto h : holes) boundary_vertices.insert(boundary_vertices.end(), h.vertices_begin(), h.vertices_end());
	std::unique(boundary_vertices.begin(),boundary_vertices.end());
	return boundary_vertices;
}

Polyline_2 PatternBoundary::filter_and_snap(Polyline_2& polyline, const Number_type& squared_dist_threshold,
		Point_2& firstPt, Point_2& lastPt, bool& snappedFirst, bool& snappedLast) {
	std::vector<Point_2> pts; polyline_to_points(polyline, pts);
	std::vector<Point_2> filtered_pts;
	// special case for one segment
	if (pts.size() > 2 ) {
		for (auto pt: pts) {
			if (inside(pt)) filtered_pts.push_back(pt);
			else {std::cout << "filtered pt " << pt << std::endl;}
			//filtered_pts.push_back(pt);
		}
	} else {
		filtered_pts = pts;
	}
	filtered_pts[0] = snap_pt(filtered_pts[0], squared_dist_threshold, snappedFirst);
	filtered_pts.back() = snap_pt(filtered_pts.back(), squared_dist_threshold, snappedLast);
	Geom_traits_2 traits;
    Geom_traits_2::Construct_curve_2 polyline_construct = traits.construct_curve_2_object();

    firstPt = filtered_pts[0]; lastPt = filtered_pts.back();

    return polyline_construct(filtered_pts.begin(), filtered_pts.end());
}

bool PatternBoundary::inside(const Point_2& pt) {
	// True if the point is inside the outer boundary but outside the holes (and including the boundary of both)
	if (outer_boundary.bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE) return false;
	for (int i = 0; i < holes.size(); i++) if (holes[i].bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) return false;
	return true;
}

Point_2 PatternBoundary::snap_pt(const Point_2& pt, const Number_type& squared_dist_threshold, bool& has_snapped) {
	has_snapped = false;
	if (outer_boundary.bounded_side(pt) == CGAL::ON_BOUNDARY) return pt;
	for (auto poly: holes) if (poly.bounded_side(pt) == CGAL::ON_BOUNDARY) return pt;

	Point_2 snappedPt = pt;
	int closest_seg = 0; Number_type min_dist = CGAL::squared_distance(pt, all_polygons_edges[0].source());
	for (auto seg: all_polygons_edges) {
		Point_2 proj_pt; Number_type squared_dist;
		proj_pt_to_segment(pt,seg,proj_pt,squared_dist);
		if (squared_dist < min(squared_dist_threshold,min_dist)) {
			min_dist = squared_dist;
			snappedPt = proj_pt;
			has_snapped = true;
		}
	}
	return snappedPt;
}

std::vector<Point_2> PatternBoundary::get_vertical_and_horizontal_intersections(Point_2& pt) {
	std::vector<Point_2> pt_intersections;
	// Create a segment whose length depends on the bounding box of the boundary polygon
	// Compute the intersection of this segment with the entire segments of the polylines
	// 	CGAL::compute_intersection_points(int_segments.begin(), int_segments.end(), std::back_inserter(intersections), false, geom_traits_2);	
	//Segment_2 seg_x(pt, CGAL::Direction_2<Kernel>(1,0)), seg_y(pt, CGAL::Direction_2<Kernel>(0,1));
	//CGAL::Ray_2<Kernel> ray_x(pt, CGAL::Direction_2<Kernel>(1,0)), ray_y(pt, CGAL::Direction_2<Kernel>(0,1));

	Geom_traits_2 geom_traits_2;
	std::vector<Segment_2> segments_and_x = all_polygons_edges; segments_and_x.push_back(Segment_2(Point_2(pt.x()-10000,pt.y()),Point_2(pt.x()+10000,pt.y())));
	CGAL::compute_intersection_points(segments_and_x.begin(), segments_and_x.end(),
                                    std::back_inserter(pt_intersections), false, geom_traits_2);

	std::vector<Segment_2> segments_and_y = all_polygons_edges; segments_and_y.push_back(Segment_2(Point_2(pt.x(),pt.y()-10000),Point_2(pt.x(),pt.y()+10000)));
	CGAL::compute_intersection_points(segments_and_y.begin(), segments_and_y.end(),
                                    std::back_inserter(pt_intersections), false, geom_traits_2);

	return pt_intersections;
}

std::vector<Polygon_2> PatternBoundary::polylines_to_polygons(const std::vector<Polyline_2>& boundary_polylines) {
	std::vector<Polygon_2> polygons;
	for (auto poly:boundary_polylines) {
		std::vector<Point_2> points; polyline_to_points(poly, points);
		// Make sure the first and last point are not equal
		if (points[0] == points.back()) points.pop_back();
		polygons.push_back(Polygon_2(points.begin(),points.end()));
	}
	return polygons;
}

void PatternBoundary::proj_pt_to_segment(const Point_2& pt, const Segment_2& seg_i, Point_2& proj_pt, Number_type& squared_dist) {
	CGAL::Segment_2<Kernel> seg(seg_i.source(),seg_i.target());
	Point_2 line_proj = seg.supporting_line().projection(pt);
	if (seg.has_on(line_proj)) proj_pt = line_proj;
	else proj_pt =  CGAL::compare_distance_to_point(pt, seg.source(), seg.target()) == CGAL::SMALLER ? seg.source():seg.target(); 

	squared_dist = CGAL::squared_distance(pt,proj_pt);
}

void PatternBoundary::polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points) {
	int seg_n = poly.subcurves_end()-poly.subcurves_begin();
	points.push_back(poly.subcurves_begin()->source());
	for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
		points.push_back(seg_i->target());
	}
}