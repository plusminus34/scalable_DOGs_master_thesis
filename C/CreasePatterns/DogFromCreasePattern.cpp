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
	std::vector<Eigen::MatrixXd> submeshVList; std::vector<Eigen::MatrixXi> submeshFList;
	Eigen::MatrixXd V; Eigen::MatrixXi F; DogFoldingConstraints foldingConstraints;
	generate_mesh(creasePattern, gridPolygons, sqr_in_polygon, submeshVList, submeshFList, V, F);

	generate_constraints(creasePattern, submeshVList, submeshFList, foldingConstraints);
	Eigen::MatrixXi F_ren = generate_rendered_mesh_faces(creasePattern, V, foldingConstraints);

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
		//std::cout << "Polygon number " << face_i << " with " << poly.size() << " vertices" << std::endl;

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
					std::vector<Eigen::MatrixXd>& submeshVList, std::vector<Eigen::MatrixXi>& submeshFList,
					Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	int submesh_n = gridPolygons.size();
	submeshVList.resize(submesh_n); submeshFList.resize(submesh_n);
	
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


void generate_constraints(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
						const std::vector<Eigen::MatrixXi>& submeshFList, DogFoldingConstraints& foldingConstraints) {
	// Get all the polylines unique points.
	const std::vector<Polyline_2>& polyline_pts = creasePattern.get_clipped_polylines();
	std::set<Point_2> constrained_pts;
	for (auto poly : polyline_pts) {
		for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
			constrained_pts.insert(seg_i->source()); constrained_pts.insert(seg_i->target());
		} 
	}
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	// For each pt, perform a query on the orthogoanl grid arrangement. It can be on a vertex or an edge.
	for (auto pt: constrained_pts) {
		std::pair<Point_2,Point_2> edge_pts; double t;
		if (!orthGrid.get_pt_edge_coordinates(pt, edge_pts,t)) {
			std::cout << "Error, got a point pt = " << pt << " that is not on the grid " << std::endl;
			exit(1); // Should never get here, and if so all is lost
		}
		//std::cout << "point p = " << pt << " lies between " << edge_pts.first << " and " << edge_pts.second << " with t = " << t << std::endl;
		// Now find the indices of both points and add them as constraints
		// For every point, find all submeshes that contain it. We need to have both points for a submesh to count.
		
		Eigen::RowVector3d pt1(CGAL::to_double(edge_pts.first.x()),CGAL::to_double(edge_pts.first.y()),0);
		Eigen::RowVector3d pt2(CGAL::to_double(edge_pts.second.x()),CGAL::to_double(edge_pts.second.y()),0);
		std::vector<std::pair<int,int>> global_edge_indices;
		// Go through every submesh
		int global_idx_base = 0;
		for (int sub_i = 0; sub_i < submeshVList.size(); sub_i++) {
			// Find the indices of the pts in the submesh (if they exist)
			int pt1_idx = -1, pt2_idx = -1;
			for (int v_i = 0; v_i < submeshVList[sub_i].rows(); v_i++) {
				if (submeshVList[sub_i].row(v_i) == pt1) {pt1_idx = v_i;}
				if (submeshVList[sub_i].row(v_i) == pt2) {pt2_idx = v_i;}
				if ((pt1_idx!=-1)&& (pt2_idx != -1)) break;
			}
			// Found the edge (both points) on the submesh
			if ((pt1_idx!=-1)&& (pt2_idx != -1)) {
				// Insert the (global V) indices
				global_edge_indices.push_back(std::pair<int,int>(global_idx_base+pt1_idx,global_idx_base+pt2_idx));
			}
			global_idx_base += submeshVList[sub_i].rows();
		}
		//std::cout << "global_edge_indices.size() = " << global_edge_indices.size() << std::endl;
		// We then need to map it to the global V which just have the submeshes concatenated
		if (!global_edge_indices.size()){
			std::cout << "Error, got an edge that is not in a submesh, with pt1 = " << pt1 << " and pt2 = " << pt2 << std::endl;
			exit(1); // Should not get here, and if so there's really nothing to do but debug the crease pattern
		}
		
		// We got 'n' different vertex pairs, hence we need n-1 (circular) constraints
		for (int const_i = 0; const_i < global_edge_indices.size()-1; const_i++) {
			foldingConstraints.edge_const_1.push_back(global_edge_indices[const_i]);
			foldingConstraints.edge_const_2.push_back(global_edge_indices[const_i+1]);
			foldingConstraints.edge_coordinates.push_back(t);
		}
	}
	// 
}

Eigen::MatrixXi generate_rendered_mesh_faces(const CreasePattern& creasePattern, const Eigen::MatrixXd& V, const DogFoldingConstraints& fC) {
	Eigen::MatrixXi F_ren;
	// The rendered mesh should have the vertices of V as well as polyline vertices.
	// So V_ren = [V,V_p], where V_p is a (not necessarily unique) list of the polyline vertices, who'se values we can take from a constraints list.
	// The faces can be obtained from the clipped arrangement with the grid, which contains only these vertices. 
	// They will be written in coordinates, so we will first need to create a mapping between coordinates and indices in V_ren.
	// We can first save the polygonal faces, and call igl::triangulate which will work perfectly as everything is convex.
	// The only thing that will then be needed for rendering is to update V_p, which could be done with the constraints list.
	
	Eigen::MatrixXd V_ren; Dog::get_V_ren(V, fC, V_ren);
	std::vector<Polygon_2> faces_polygons; 
	// Get orth grid and add it the polylines
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	PlanarArrangement grid_with_snapped(orthGrid);
	grid_with_snapped.add_polylines(creasePattern.get_clipped_polylines());
	grid_with_snapped.get_faces_polygons(faces_polygons);
	
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
	//std::cout << "squares_num = " << squares_num << std::endl;
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