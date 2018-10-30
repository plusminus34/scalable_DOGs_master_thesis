#include "OrthogonalGrid.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


OrthogonalGrid::OrthogonalGrid(Bbox_2 bbox2, int x_res, int y_res, const std::vector<Point_2>& additional_grid_points) :
										bbox2(bbox) {
	set_up_initial_grid(x_res, y_res);
	for (auto pt : additional_grid_points) {
		subdivide_grid_at_pt(pt);
	}
}

void OrthogonalGrid::set_up_initial_grid(int x_res, int y_res) {
	create_spaced_range(bbox2.xmin(), bbox2.xmax(), x_res, x_coords);
	create_spaced_range(bbox2.ymin(), bbox2.ymax(), y_res, y_coords);
}

void OrthogonalGrid::subdivide_grid_at_pt(const Point_2& pt) {
	if (is_point_on_grid(pt)) return;
	// Todo: check what is better x or y
	if (1) {
		x_coords.push_back(pt.x());
		std::sort(x_coords.begin(),x_coords.end());
	} else {
		x_coords.push_back(pt.y());
		std::sort(y_coords.begin(),y_coords.end());
	}	
}

bool OrthogonalGrid::is_point_on_grid(const Point_2& pt) {
	// TODO: implement
	return false;
}

void OrthogonalGrid::create_spaced_range(const Number_type min, const Number_type max, const int num_points, std::vector<Number_type>& range) {
	Number_type spacing = (max-min)/(num_points-1);
	for (int i = 0; i < num_points; i++) {range.push_back(min+i*spacing);}
}