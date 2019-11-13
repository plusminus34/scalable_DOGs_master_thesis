#pragma once

#include "ArrangementDefs.h"
#include <string>

void read_svg_crease_pattern(const std::string& path, CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines,
				std::vector<Polyline_2>& boundary_polylines, std::vector<Eigen::MatrixXd>& polylines_data,
			  std::vector<Eigen::MatrixXd>& boundary_polylines_data);

//Polyline_2 points_to_polylines_snapped_at_start_end(const Eigen::MatrixXd& p, const CGAL::Bbox_2& bbox);
Polyline_2 points_to_polylines(const Eigen::MatrixXd& p);

Point_2 snap_pt_to_bbox(const Point_2& pt, const CGAL::Bbox_2& bbox);
void round_if_close(Number_type& pt, const Number_type& floor_pt);
