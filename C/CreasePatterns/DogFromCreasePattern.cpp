#include "DogFromCreasePattern.h"

#include <igl/combine.h>
#include <igl/remove_unreferenced.h>

Dog dog_from_crease_pattern(const CreasePattern& creasePattern) {
	std::vector<Polygon_2> gridPolygons;
	init_grid_polygons(creasePattern, gridPolygons);
	// Per polygon contains a flag (whether face 'i' intersects that polygon)
	std::vector<std::vector<bool>> sqr_in_polygon;
	set_sqr_in_polygon(creasePattern, gridPolygons, sqr_in_polygon);

	// Generated mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F; DogFoldingConstraints foldingConstraints;
	generate_mesh(creasePattern, gridPolygons, sqr_in_polygon, V, F);

	generate_constraints(foldingConstraints);
	Eigen::MatrixXi F_ren = generate_rendered_mesh_faces();

	return Dog(V,F,foldingConstraints,F_ren);
}

void set_sqr_in_polygon(const CreasePattern& creasePattern, std::vector<Polygon_2>& gridPolygons, 
						std::vector<std::vector<bool>>& sqr_in_polygon) {
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
			//std::cout << "face " << f_i << " in polygon = " << sqr_in_polygon[face_i][f_i] << std::endl;
		}
		face_i++;
	}
}

void generate_mesh(const CreasePattern& creasePattern, const std::vector<Polygon_2>& gridPolygons, const std::vector<std::vector<bool>>& sqr_in_polygon,
					Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	int submesh_n = gridPolygons.size();
	std::vector<Eigen::MatrixXd> submeshVList(submesh_n); std::vector<Eigen::MatrixXi> submeshFList(submesh_n);
	
	Eigen::MatrixXd gridV; Eigen::MatrixXi gridF; 
	init_mesh_vertices_and_faces_from_grid(creasePattern, gridV, gridF);
	int poly_idx = 0;
	for (auto submesh_flags: sqr_in_polygon ) {
		Eigen::MatrixXd submeshV; Eigen::MatrixXi submeshF;
		int f_cnt = std::count(submesh_flags.begin(), submesh_flags.end(), true);
		submeshF.resize(f_cnt,4); int cnt = 0;
		for (int fi = 0; fi < submesh_flags.size(); fi++) {
			if (submesh_flags[fi]) {
				submeshF.row(cnt++) << gridF.row(fi);	
			} 
		}
		Eigen::MatrixXi IV;
		igl::remove_unreferenced(gridV,submeshF,submeshV,submeshF,IV);
		submeshVList[poly_idx] = submeshV;
		submeshFList[poly_idx] = submeshF;
		poly_idx++;
	}
	igl::combine(submeshVList,submeshFList, V, F);
}


void generate_constraints(DogFoldingConstraints& foldingConstraints) {
	// TODO
}

Eigen::MatrixXi generate_rendered_mesh_faces() {
	Eigen::MatrixXi F_ren;
	// The rendered mesh should have the vertices of V as well as polyline vertices.
	// So V_ren = [V,V_p], where V_p is a (not necessarily unique) list of the polyline vertices, who'se values we can take from a constraints list.
	// The faces can be obtained from the clipped arrangement, which contains only these vertices. 
	// They will be written in coordinates, so we will first need to create a mapping between coordinates and indices in V_ren.
	// We can first save the polygonal faces, and call igl::triangulate which will work perfectly as everything is convex.
	// The only thing that will then be needed for rendering is to update V_p, which could be done with the constraints list.

	return F_ren;
}


void init_mesh_vertices_and_faces_from_grid(const CreasePattern& creasePattern, Eigen::MatrixXd& gridV, Eigen::MatrixXi& gridF) {
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	const std::vector<Number_type>& gx_coords(orthGrid.get_x_coords()), gy_coords(orthGrid.get_y_coords());
	gridV.resize(gx_coords.size()*gy_coords.size(),3); gridF.resize((gx_coords.size()-1)*(gy_coords.size()-1),4);
	int fcnt = 0;
	for (int y_i = 0; y_i < gy_coords.size(); y_i++) {
		for (int x_i = 0; x_i < gx_coords.size(); x_i++) {
			gridV.row(y_i*(gx_coords.size())+x_i) << CGAL::to_double(gx_coords[x_i]),CGAL::to_double(gy_coords[y_i]),0;
			if ((x_i < gx_coords.size()-1) && (y_i < gy_coords.size()-1)) {
				// In DDG index notation add the face [F,F_1,F_12,F_2]
				gridF.row(fcnt) << y_i*gx_coords.size()+x_i, y_i*gx_coords.size()+x_i+1,(y_i+1)*gx_coords.size()+x_i+1,(y_i+1)*gx_coords.size()+x_i;
				fcnt++;
			}
		}
	}
}

void init_grid_polygons(const CreasePattern& creasePattern,std::vector<Polygon_2>& gridPolygons) {
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