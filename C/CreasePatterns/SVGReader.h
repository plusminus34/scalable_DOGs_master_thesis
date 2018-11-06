#include "ArrangementDefs.h"
#include <string>

void read_svg_crease_pattern(const std::string& path, CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines);

// Used for constructors
CGAL::Bbox_2 bbox_from_svg_crease_pattern(const std::string& path);
std::vector<Polyline_2> polylines_from_svg_crease_pattern(const std::string& path);

Polyline_2 points_to_polylines_snapped_at_start_end(const Eigen::MatrixXd& p, const CGAL::Bbox_2& bbox);

Point_2 snap_pt_to_bbox(const Point_2& pt, const CGAL::Bbox_2& bbox);
void round_if_close(Number_type& pt, const Number_type& floor_pt);