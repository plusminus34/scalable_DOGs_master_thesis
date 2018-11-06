#include "DogBuilder.h"

#include "SVGReader.h"

DogBuilder::DogBuilder(const DogCreasePattern& i_dogCreasePattern) : creasePattern(i_dogCreasePattern) {
	initialize();
}

void DogBuilder::set_sqr_in_polygon() {
	// Get the faces' polygons and find their intersections with the faces
	std::vector<Polygon_2> facePolygons; 
	creasePattern.get_clipped_arrangement().get_faces_polygons(facePolygons);

	sqr_in_polygon.resize(facePolygons.size());
	// Iterate over the polygons and add faces that intersect
	int face_i = 0;
	for (auto poly: facePolygons) {
		sqr_in_polygon[face_i].resize(gridPolygons.size());
		std::cout << "Polygon number " << face_i << " with " << poly.size() << " vertices" << std::endl;
		for (int f_i = 0; f_i < gridPolygons.size(); f_i++) {
			sqr_in_polygon[face_i][f_i] = CGAL::do_intersect(poly, gridPolygons[f_i]);
			std::cout << "face " << f_i << " in polygon = " << sqr_in_polygon[face_i][f_i] << std::endl;
		}
		face_i++;
	}
	// Implement a function that creates a mesh given a grid and a boolean flag for each face (whether it's in or out)
	// Count the occurence of each faces, those with one occurences are inner faces and are usefull later for rendering
	// Return the meshes
}

void DogBuilder::initialize() {
	init_grid_polygons();
	set_sqr_in_polygon();
}

void DogBuilder::init_grid_polygons() {
	// Squares num are the faces num -1 for the exterior face
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	const std::vector<Number_type>& gx_coords(orthGrid.get_x_coords()), gy_coords(orthGrid.get_y_coords());
	int squares_num = (gx_coords.size()-1)*(gy_coords.size()-1);
	std::cout << "squares_num = " << squares_num << std::endl;
	gridPolygons.resize(squares_num);
	
	for (int y_i = 0; y_i < gy_coords.size()-1; y_i++) {
		for (int x_i = 0; x_i < gx_coords.size()-1; x_i++) {
			Polygon_2 P;
			// Add points in counter clockwise ordering
  			P.push_back(Point_2(gx_coords[x_i],gy_coords[y_i])); // 0,0
  			P.push_back(Point_2(gx_coords[x_i+1],gy_coords[y_i])); // 1,0
  			P.push_back(Point_2(gx_coords[x_i+1],gy_coords[y_i+1])); // 1,1
  			P.push_back(Point_2(gx_coords[x_i],gy_coords[y_i+1])); // 0,1

  			gridPolygons[y_i*(gx_coords.size()-1)+x_i] = P;
  			//std::cout << "y_i*(gx_coords.size()-1)+x_i = " <<  << std::endl;
		}
	}
}