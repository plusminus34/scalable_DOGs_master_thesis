#include "DogFromCreasePattern.h"

#include <igl/combine.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/remove_unreferenced.h>

#include <CGAL/Polygon_2_algorithms.h>

#include "../QuadMesh/Quad.h"

Dog dog_from_crease_pattern(const CreasePattern& creasePattern) {
	std::vector<Polygon_2> gridPolygons;
	init_grid_polygons(creasePattern, gridPolygons);
	// Per polygon contains a flag (whether face 'i' intersects that polygon)
	std::vector<std::vector<bool>> sqr_in_polygon;
	set_sqr_in_polygon(creasePattern, gridPolygons, sqr_in_polygon);

	// Generated mesh
	std::vector<Eigen::MatrixXd> submeshVList; std::vector<Eigen::MatrixXi> submeshFList;
	std::vector<std::vector<bool>> submeshV_is_inner;
	Eigen::MatrixXd V; Eigen::MatrixXi F; DogEdgeStitching edgeStitching;
	generate_mesh(creasePattern, gridPolygons, sqr_in_polygon, submeshVList, submeshFList, submeshV_is_inner, V, F);

	std::vector<Point_2> constrained_pts_non_unique;
	generate_constraints(creasePattern, submeshVList, submeshFList, edgeStitching, constrained_pts_non_unique,V);

	std::vector<SubmeshPoly> submesh_polygons;
	get_faces_partitions_to_submeshes(creasePattern, submesh_polygons);
	Eigen::MatrixXd V_ren; Dog::V_ren_from_V_and_const(V, edgeStitching, V_ren);
	Eigen::MatrixXi F_ren = generate_rendered_mesh_faces(creasePattern, submesh_polygons, submeshVList, V_ren, constrained_pts_non_unique);

	std::vector<int> submeshVSize(submeshVList.size());
	for (int i = 0; i < submeshVList.size(); i++) {submeshVSize[i] = submeshVList[i].rows();}

	return Dog(V,F,edgeStitching,V_ren,F_ren, submeshVSize);
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
			// This could be do the inaccuracies converting CGAL to double? (I think there's no conversion here so this is weird)
			// If we do intersect, we filter those numerical errors
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
					std::vector<Eigen::MatrixXd>& submeshVList, std::vector<Eigen::MatrixXi>& submeshFList, std::vector<std::vector<bool>>& submeshV_is_inner,
					Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	int submesh_n = gridPolygons.size();
	submeshVList.resize(submesh_n); submeshFList.resize(submesh_n); submeshV_is_inner.resize(submesh_n);

	
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
		Eigen::MatrixXi I,J; // J is s.t slice(V,J,1,NV);
		igl::remove_unreferenced(gridV,submeshF,submeshV,submeshF,I,J);
		submeshVList[poly_idx] = submeshV;
		submeshFList[poly_idx] = submeshF;

		std::cout << "poly_idx = " << poly_idx << std::endl;
		submeshV_is_inner[poly_idx].resize(submeshV.rows());
		
		poly_idx++;
	}
	igl::combine(submeshVList,submeshFList, V, F);
}



// Another way to go about it is to get the vertices of the )orth grid + polylines) graph that share more than one face
void generate_constraints(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
						const std::vector<Eigen::MatrixXi>& submeshFList, DogEdgeStitching& edgeStitching,
						std::vector<Point_2>& constrained_pts_non_unique, const Eigen::MatrixXd& V) {
	// Get all the polylines unique points (vertices will appear twice with each polyline)
	std::set<Point_2> constrained_pts;
	const std::vector<Polyline_2>& polylines = creasePattern.get_clipped_polylines();
	for (auto poly : polylines) {
		//std::cout << "new poly"<<std::endl;
		std::vector<Point_2> pts; polyline_to_points(poly,pts);
		constrained_pts.insert(pts.begin(),pts.end());
	}
	// For each pt, perform a query on the orthogoanl grid arrangement. It can be on a vertex or an edge.
	for (auto pt: constrained_pts) {

		Number_type edge_t; std::vector<Edge> edge_v_indices;
		pt_to_edge_coordiantes(pt, creasePattern, submeshVList, edge_v_indices, edge_t);

		// We got 'n' different vertex pairs, hence we need n-1 (circular) constraints
		for (int const_i = 0; const_i < edge_v_indices.size()-1; const_i++) {
			edgeStitching.edge_const_1.push_back(edge_v_indices[const_i]);
			edgeStitching.edge_const_2.push_back(edge_v_indices[const_i+1]);
			edgeStitching.edge_coordinates.push_back(CGAL::to_double(edge_t));
			edgeStitching.edge_coordinates_precise.push_back(edge_t);
			constrained_pts_non_unique.push_back(pt);
		}
	}
	
	// set stitched_curves (once per edge, even if it is duplicated)
	edgeStitching.stitched_curves.resize(polylines.size());
	
	for (int i = 0; i < polylines.size(); i++) {
		std::vector<Point_2> pts; polyline_to_points(polylines[i],pts);
		edgeStitching.stitched_curves[i].resize(pts.size());
		for (int j = 0; j < pts.size(); j++) {
			Number_type edge_t; std::vector<Edge> edge_v_indices;
			pt_to_edge_coordiantes(pts[j], creasePattern, submeshVList, edge_v_indices, edge_t);

			auto edge_v = CGAL::to_double(edge_t)*V.row(edge_v_indices[0].v1)+(1-CGAL::to_double(edge_t))*V.row(edge_v_indices[0].v2);
			// Choose one of the edges, don't need the duplicates and they are all (approximately) equal anyhow
			edgeStitching.stitched_curves[i][j] = EdgePoint(edge_v_indices[0],CGAL::to_double(edge_t));
		}
	}
}

void get_faces_partitions_to_submeshes(const CreasePattern& creasePattern, std::vector<SubmeshPoly>& submesh_polygons) {
	std::vector<Polygon_2> faces_polygons; std::vector<int> faces_to_submesh;

	// Get orth grid and add it the polylines
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	PlanarArrangement grid_with_snapped(orthGrid);
	grid_with_snapped.add_polylines(creasePattern.get_clipped_polylines());
	grid_with_snapped.get_faces_polygons(faces_polygons);


	std::vector<Polygon_2> submeshBnd; 
	creasePattern.get_clipped_arrangement().get_faces_polygons(submeshBnd);

	// The rendered mesh faces will be the polygons in grid_with_snapped (meaning faces_polygons) after triangulation.
	// Here we'll find which submesh has them, and then translate the point indices to the correct point in something
	std::vector<int> face_to_submesh(faces_polygons.size());
	for (int f_i = 0; f_i < faces_polygons.size(); f_i++){
		auto face_polygon = faces_polygons[f_i];
		// Now go through each submeshPolygon and check if its there
		face_to_submesh[f_i] = -1;
		//std::cout << "----- Checking face " << face_polygon << " ---------" << std::endl;
		int submesh_i = 0;
		while ( (submesh_i < submeshBnd.size()) && (face_to_submesh[f_i] == -1) ) {
			bool is_in_submesh = true;
			//std::cout << "checking if its in polygon = " << submeshBnd[submesh_i] << std::endl;
			for (auto vptr = face_polygon.vertices_begin(); vptr != face_polygon.vertices_end(); vptr++) {
				auto pt_in_face = ! (submeshBnd[submesh_i].bounded_side(*vptr) == CGAL::ON_UNBOUNDED_SIDE);
				//std::cout << "pt " << *vptr << " in face = " << pt_in_face << std::endl;
				is_in_submesh = is_in_submesh & (pt_in_face);
			}
			if (is_in_submesh) face_to_submesh[f_i] = submesh_i;
			submesh_i++;
		}

		if (face_to_submesh[f_i] == -1) {std::cout << "Error, couldn't find a submesh for face = " << face_polygon << std::endl; exit(1); /*nothing to do but debug*/}
		//std::cout << "face in submesh = " << face_to_submesh[f_i] << std::endl;

	}
	submesh_polygons.resize(faces_polygons.size());
	for (int i = 0; i < submesh_polygons.size(); i++) {submesh_polygons[i] = SubmeshPoly(face_to_submesh[i], faces_polygons[i]);}
	// sort it by the submeshes indices
	std::sort(submesh_polygons.begin(), submesh_polygons.end(), [](const SubmeshPoly& p1, const SubmeshPoly& p2){return p1.first < p2.first;});
	//for (auto sp: submesh_polygons) {std::cout << "submesh = " << sp.first << " polygon = " << sp.second << std::endl;}
}

Eigen::MatrixXi generate_rendered_mesh_faces(const CreasePattern& creasePattern, std::vector<SubmeshPoly>& submesh_polygons,
			const std::vector<Eigen::MatrixXd>& submeshVList, const Eigen::MatrixXd& V_ren, const std::vector<Point_2>& constrained_pts_non_unique) {
	typedef std::pair<double,double> PointDouble;
	std::vector<std::vector<int> > polygons_v_ren_indices(submesh_polygons.size());


	// IMPORTANT:
	// To make things a (bit) faster, this method uses a mapping from Points (in double) to a given submesh
	// The function then iterates on CGAL exact kernel Point_2, convert them to double with CGAL::to_double() and then use this map
	//	This should always work for points on the grid (since the mesh points in double with the same CGAL::to_double routing)
	//	The only problem might be edge points which are in the exact form t*v1+(1-t)*v2.
	//	V_ren has the points in the from to_double(t)*to_double(v1) + to_double(1-t)*to_double(v2) which is different than converting it all together, i.e:
	//			to_double( t*v1 + (1-t)*v2  )
	//  Therefore the solution is to use the original polygon points which we saved at constrained_pts_non_unique
	
	// set [min_sub_i, max_sub_i) as the (half-closed) range of indices of the submesh in subPoly inside V_ren 
	int sub_i = 0; int min_sub_i = 0; int max_sub_i = min_sub_i + submeshVList[sub_i].rows();

	// Save a mapping between constrained fold points to V_ren indices
	std::map<PointDouble, int> poly_fold_pt_to_V_ren; int v_ren_start_idx = 0;
	for (auto subV: submeshVList) v_ren_start_idx+=subV.rows();
	for (int const_i = 0; const_i < constrained_pts_non_unique.size(); const_i++) {
		auto const_pt = constrained_pts_non_unique[const_i];
		PointDouble pt(CGAL::to_double(const_pt.x()),CGAL::to_double(const_pt.y()));
		poly_fold_pt_to_V_ren[pt] = v_ren_start_idx+const_i;
	}
	std::map<PointDouble, int> submesh_pt_to_V_ren;
	// create a map from point coordinates to index in V_ren (which should be a vertex in the correct submesh)
	for (int ri = min_sub_i; ri < max_sub_i; ri++) {submesh_pt_to_V_ren[PointDouble(V_ren(ri,0),V_ren(ri,1))] = ri;}
	std::cout <<"sub_i = " << sub_i << std::endl;

	for (int poly_i = 0; poly_i < submesh_polygons.size(); poly_i++) {
		auto submesh_i = submesh_polygons[poly_i].first; auto submesh_poly = submesh_polygons[poly_i].second;
		// Update the interval [min_sub_i, max_sub_i) to the new submesh vertices indices in V_ren
		if (submesh_i != sub_i) {
			sub_i++;
			min_sub_i = max_sub_i;
			max_sub_i = min_sub_i + submeshVList[sub_i].rows();
			
			submesh_pt_to_V_ren.clear();
			// update the mapping from points to coordinates in v_ren for the new submesh
			for (int ri = min_sub_i; ri < max_sub_i; ri++) {submesh_pt_to_V_ren[PointDouble(V_ren(ri,0),V_ren(ri,1))] = ri;}
			std::cout <<"sub_i = " << sub_i << std::endl;
		}
		// every pt should have some cashed indices to V_ren
		polygons_v_ren_indices[poly_i].resize(submesh_poly.size()); int poly_vi = 0;
		for (auto vptr = submesh_poly.vertices_begin(); vptr != submesh_poly.vertices_end(); vptr++) {
			PointDouble pt(CGAL::to_double(vptr->x()),CGAL::to_double(vptr->y()));
			
			if (submesh_pt_to_V_ren.count(pt)) {
				std::cout << "pt has key in submesh "  << std::endl;
				polygons_v_ren_indices[poly_i][poly_vi] = submesh_pt_to_V_ren[pt];
			} else if (poly_fold_pt_to_V_ren.count(pt)) {
				std::cout << "pt has key in fold map" << std::endl;
				polygons_v_ren_indices[poly_i][poly_vi] = poly_fold_pt_to_V_ren[pt];
			} else {
				std::cout << "Error! Could not map point to polygon " << std::endl;

				exit(1); // No solution to this other than debugging
			}
			poly_vi++;
		}
	}

	// Generate a triangular mesh from the polygonal faces (which are all convex, otherwise this routine is problematic)
	Eigen::MatrixXi F_ren;
	igl::polygon_mesh_to_triangle_mesh(polygons_v_ren_indices,F_ren);
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

void pt_to_edge_coordiantes(const Point_2& pt, const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
				std::vector<Edge>& edge_v_indices, Number_type& edge_t) {

	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	std::pair<Point_2,Point_2> edge_pts; 
	if (!orthGrid.get_pt_edge_coordinates(pt, edge_pts,edge_t)) {
		std::cout << "Error, got a point pt = " << pt << " that is not on the grid " << std::endl;
		exit(1); // Should never get here, and if so all is lost
	}
	//std::cout << "point p = " << pt << " lies between " << edge_pts.first << " and " << edge_pts.second << " with t = " << t << std::endl;
	// Now find the indices of both points and add them as constraints
	// For every point, find all submeshes that contain it. We need to have both points for a submesh to count.
	
	Eigen::RowVector3d pt1(CGAL::to_double(edge_pts.first.x()),CGAL::to_double(edge_pts.first.y()),0);
	Eigen::RowVector3d pt2(CGAL::to_double(edge_pts.second.x()),CGAL::to_double(edge_pts.second.y()),0);
	
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
			edge_v_indices.push_back(Edge(global_idx_base+pt1_idx,global_idx_base+pt2_idx));
		}
		global_idx_base += submeshVList[sub_i].rows();
	}
	//std::cout << "global_edge_indices.size() = " << global_edge_indices.size() << std::endl;
	// We then need to map it to the global V which just have the submeshes concatenated
	if (!edge_v_indices.size()){
		std::cout << "Error, got an edge that is not in a submesh, with pt1 = " << pt1 << " and pt2 = " << pt2 << std::endl;
		exit(1); // Should not get here, and if so there's really nothing to do but debug the crease pattern
	}	
}

bool is_closed_polyline(const Polyline_2& poly) {
	int seg_n = poly.subcurves_end()-poly.subcurves_begin();
	auto first_pt = poly.subcurves_begin()->source(), last_pt = (poly.subcurves_begin()+(seg_n-1))->target();
	bool is_closed = (first_pt == last_pt);
	return is_closed;
}

void polyline_to_points(const Polyline_2& poly, std::vector<Point_2>& points) {
	// handle closed and open curves differently
	int seg_n = poly.subcurves_end()-poly.subcurves_begin();
	int points_n = (is_closed_polyline(poly)) ? seg_n : seg_n+1;
	
	points.resize(points_n); int cnt = 0;
	points[cnt++] = poly.subcurves_begin()->source();
	for (auto seg_i = poly.subcurves_begin(); seg_i!= poly.subcurves_end(); seg_i++) {
		points[cnt++] = seg_i->target();
	}
}