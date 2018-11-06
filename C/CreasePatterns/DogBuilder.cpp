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
		sqr_in_polygon[face_i] = std::vector<bool>(gridPolygons.size(), false);
		std::cout << "Polygon number " << face_i << " with " << poly.size() << " vertices" << std::endl;

		for (int f_i = 0; f_i < gridPolygons.size(); f_i++) {
			bool face_intersection = CGAL::do_intersect(poly, gridPolygons[f_i]);

			// NOTE: Minor inaccuracies (in CGAL??) cause vertex intersections to sometime return a polygon with a very small area (1e-28)
			// If we do intersect, we filter those
			if (face_intersection) {
				Polygon_set R;
				CGAL::intersection(poly, gridPolygons[f_i], std::back_inserter(R));
				bool all_areas_are_zero = true;
  				for (auto rit = R.begin(); rit != R.end(); ++rit) {
					auto outer_bnd_poly = rit->outer_boundary();
					all_areas_are_zero = all_areas_are_zero & (outer_bnd_poly.area() < 1e-20);
  				}
  				//std::cout << " }" << std::endl;
  				//std::cout << "all_areas_are_zero = " << all_areas_are_zero << std::endl;
  				face_intersection = !all_areas_are_zero;
			}
			sqr_in_polygon[face_i][f_i] = face_intersection;
			std::cout << "face " << f_i << " in polygon = " << sqr_in_polygon[face_i][f_i] << std::endl;
		}
		face_i++;
	}
}

void DogBuilder::generate_mesh() {
	submesh_n = gridPolygons.size();
	submeshV.resize(submesh_n); submeshF.resize(submesh_n);

	Eigen::MatrixXd gridV; Eigen::MatrixXi gridF; 
	init_mesh_vertices_and_faces_from_grid(gridV, gridF);
	
	for (auto submesh_flags: sqr_in_polygon ) {
		Eigen::MatrixXd submeshV; Eigen::MatrixXi submeshF;
		int f_cnt = std::count(submesh_flags.begin(), submesh_flags.end(), true);
		std::cout << "f_cnt = " << f_cnt << std::endl;
		submeshF.resize(f_cnt,4); int cnt = 0;
		for (int fi = 0; fi < submesh_flags.size(); fi++) {
			if (submesh_flags[fi]) submeshF.row(cnt++) << gridF.row(fi);
		}
	}
}

void DogBuilder::init_mesh_vertices_and_faces_from_grid(Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF) {
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	const std::vector<Number_type>& gx_coords(orthGrid.get_x_coords()), gy_coords(orthGrid.get_y_coords());
	gridV.resize(gx_coords.size()*gy_coords.size(),3); gridF.resize((gx_coords.size()-1)*(gy_coords.size()-1),4);
	int fcnt = 0;
	for (int y_i = 0; y_i < gy_coords.size(); y_i++) {
		for (int x_i = 0; x_i < gx_coords.size(); x_i++) {
			gridV.row(y_i*(gx_coords.size())+x_i) << CGAL::to_double(gx_coords[x_i]),CGAL::to_double(gy_coords[y_i]),0;

			if ((x_i < gx_coords.size()) && (y_i < gy_coords.size())) {
				// In DDG index notation add the face [F,F_1,F_12,F_2]
				gridF.row(fcnt) << y_i*gx_coords.size()+x_i, y_i*gx_coords.size()+x_i+1,(y_i+1)*gx_coords.size()+x_i+1,(y_i+1)*gx_coords.size()+x_i;
				fcnt++;
			}
		}
	}
}

void DogBuilder::initialize() {
	init_grid_polygons();
	set_sqr_in_polygon();
	generate_mesh();
	//generate_rendered_mesh();
	//generate_constraints();

	//Go through the polyline vertices, create a unique list of them, locate them in the global mesh.
	//Each occurence has to be in a different submesh, and we need to constrain couples of these occurnces (so k occurence means k-1 "circular" constraints)
	//generate_constraints();
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