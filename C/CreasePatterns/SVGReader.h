#include "ArrangementDefs.h"
#include <string>

void read_svg_crease_pattern(const std::string& path, CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines);

Polyline_2 eigen_to_polyline_2(const Eigen::MatrixXd& p);