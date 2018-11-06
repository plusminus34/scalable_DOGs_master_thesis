#include "DogBuilder.h"

#include "SVGReader.h"

DogBuilder::DogBuilder(const std::string& svgPath, int x_res, int y_res, bool snap_rounding) : 
		creasePattern(bbox_from_svg_crease_pattern(svgPath), polylines_from_svg_crease_pattern(svgPath), x_res, y_res, snap_rounding) {
	initialize();
}

DogBuilder::DogBuilder(const CGAL::Bbox_2& bbox, std::vector<Polyline_2>& polylines, int x_res, int y_res, bool snap_rounding) :
		creasePattern(bbox, polylines, x_res, y_res, snap_rounding) {
	initialize();
}


void DogBuilder::get_submeshes(std::vector<Eigen::MatrixXd>& V, std::vector<Eigen::MatrixXi>& F) {
	// Create faces polygons and find their intersections with the faces

	// Call get_faces_polygons on the clipped arrangement
	// Iterate over the polygons and add faces that intersect
	// Implement a function that creates a mesh given a grid and a boolean flag for each face (whether it's in or out)
	// Count the occurence of each faces, those with one occurences are inner faces and are usefull later for rendering
	// Return the meshes
}

void DogBuilder::initialize() {
	init_grid_polygons();	
}
void DogBuilder::init_grid_polygons() {
	// Squares num are the faces num -1 for the exterior face
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	const std::vector<Number_type>& gx_coords(orthGrid.get_x_coords()), gy_coords(orthGrid.get_y_coords());
	int squares_num = (gx_coords.size()-1)*(gy_coords.size()-1);
	gridPolygons.resize(squares_num);
	
	for (int y_i = 0; y_i < gy_coords.size()-1; y_i++) {
		for (int x_i = 0; x_i < gx_coords.size()-1; x_i++) {
			Polygon_2 P;
			// Add points in counter clockwise ordering
  			P.push_back(Point_2(gx_coords[x_i],gy_coords[y_i])); // 0,0
  			P.push_back(Point_2(gx_coords[x_i+1],gy_coords[y_i])); // 1,0
  			P.push_back(Point_2(gx_coords[x_i+1],gy_coords[y_i+1])); // 1,1
  			P.push_back(Point_2(gx_coords[x_i],gy_coords[y_i+1])); // 0,1
  			gridPolygons[y_i*gx_coords.size()+x_i] = P;
		}
	}
}