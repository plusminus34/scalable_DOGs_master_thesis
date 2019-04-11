#include "DogFromCreasePattern.h"

#include <igl/combine.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/remove_unreferenced.h>

#include <CGAL/Polygon_2_algorithms.h>

#include "../QuadMesh/Quad.h"

using namespace std;

Dog dog_from_crease_pattern(const CreasePattern& creasePattern) {
	std::vector<bool> is_polygon_hole = submesh_is_hole(creasePattern);
	std::vector<Polygon_2> gridPolygons;
	init_grid_polygons(creasePattern, gridPolygons);
	cout << "init grid polygons done" << endl;
	// Per polygon contains a flag (whether face 'i' intersects that polygon)
	std::vector<std::vector<bool>> sqr_in_polygon;
	set_sqr_in_polygon(creasePattern, is_polygon_hole, gridPolygons, sqr_in_polygon);
	cout << "set sqr in polygone done" << endl;

	// Generated mesh
	std::vector<Eigen::MatrixXd> submeshVList; std::vector<Eigen::MatrixXi> submeshFList;
	std::vector<std::vector<bool>> submeshV_is_inner;
	Eigen::MatrixXd V; Eigen::MatrixXi F; DogEdgeStitching edgeStitching;
	generate_mesh(creasePattern, gridPolygons, sqr_in_polygon, submeshVList, submeshFList, submeshV_is_inner, V, F);
	std::vector<SubmeshPoly> submesh_polygons;
	get_faces_partitions_to_submeshes(creasePattern, is_polygon_hole, submesh_polygons);
	std::cout << "partitioned faces to submeshes" << std::endl;

	std::vector<Point_2> constrained_pts_non_unique;
	generate_constraints(creasePattern, submeshVList, submeshFList, edgeStitching, constrained_pts_non_unique,V);
	save_submesh_bnd_edge_points(creasePattern, submeshVList, edgeStitching);
	cout << "saved submesh bnd edge points" << endl;
	std::vector<Eigen::MatrixXd> V_ren_list; generate_V_ren_list(V, submeshVList,edgeStitching,V_ren_list);
	Eigen::MatrixXd V_ren2; Eigen::MatrixXi F_ren2; generate_rendered_mesh_vertices_and_faces(creasePattern, submesh_polygons, 
									V_ren_list, edgeStitching, V_ren2, F_ren2);
	std::cout << "Generated constraints" << std::endl;

	std::vector<int> submeshVSize(submeshVList.size()); std::vector<int> submeshFSize(submeshFList.size());
	for (int i = 0; i < submeshVList.size(); i++) {
		submeshVSize[i] = submeshVList[i].rows();
		submeshFSize[i] = submeshFList[i].rows();
	}
	cout << "Have set submeshVSize and Fsize" << endl;
	vector< vector<int> > submesh_adjacency;
	creasePattern.get_clipped_arrangement().get_faces_adjacency_list(submesh_adjacency);
	cout << "got submesh_adjacency" << endl;

	//return Dog(V,F,edgeStitching,V_ren,F_ren, submeshVSize, submeshFSize, submesh_adjacency);
	return Dog(V,F,edgeStitching,V_ren2,F_ren2, submeshVSize, submeshFSize, submesh_adjacency);
}

void set_sqr_in_polygon(const CreasePattern& creasePattern, std::vector<bool>& is_polygon_hole, std::vector<Polygon_2>& gridPolygons, 
						std::vector<std::vector<bool>>& sqr_in_polygon) {
	// Get the faces' polygons and find their intersections with the faces
	std::vector<Polygon_with_holes_2> facePolygons; creasePattern.get_submeshes_faces_polygons(facePolygons);
	std::cout << "facePolygons.size() = " << facePolygons.size() << std::endl;

	sqr_in_polygon.resize(facePolygons.size());
	// Iterate over the polygons and add faces that intersect
	int face_i = 0;
	for (auto poly: facePolygons) {
		if (is_polygon_hole[face_i]) {face_i++; continue;}
		sqr_in_polygon[face_i] = std::vector<bool>(gridPolygons.size(), false);
		//std::cout << "Polygon number " << face_i << " with " << poly.size() << " vertices" << std::endl;

		for (int f_i = 0; f_i < gridPolygons.size(); f_i++) {
			//std::cout << "before" << std::endl;
			//std::cout << "poly = " << poly << std::endl;
			//std::cout << "gridPolygons[f_i] = " << gridPolygons[f_i] << std::endl;
			bool face_intersection = CGAL::do_intersect(poly, gridPolygons[f_i]);
			//std::cout << "after" << std::endl;
			//exit(1);

			// NOTE: Minor inaccuracies (in CGAL??) cause vertex intersections to sometime return a polygon with a very small area (1e-28)
			// This could be do the inaccuracies converting CGAL to double? (I think there's no conversion here so this is weird)
			// If we do intersect, we filter those numerical errors
			if (face_intersection) {
				cout << "got face intersection" << endl;
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
	int grid_polygons_n = gridPolygons.size();
	int submesh_n = sqr_in_polygon.size();
	submeshVList.resize(submesh_n); submeshFList.resize(submesh_n); submeshV_is_inner.resize(grid_polygons_n);

	
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
	const std::vector<Polyline_2>& polylines = creasePattern.get_clipped_fold_polylines();
	for (auto poly : polylines) {
		//std::cout << "new poly"<<std::endl;
		std::vector<Point_2> pts; polyline_to_points(poly,pts);
		constrained_pts.insert(pts.begin(),pts.end());
	}
	edgeStitching.submesh_to_edge_pt.resize(submeshFList.size());
	// For each pt, perform a query on the orthogoanl grid arrangement. It can be on a vertex or an edge.
	int pt_const_i = 0; int edge_constraints_cnt = 0;
	for (auto pt: constrained_pts) {

		Number_type edge_t; std::vector<Edge> edge_v_indices; std::vector<int> submeshes_with_pt;
		pt_to_edge_coordinates(pt, creasePattern, submeshVList, edge_v_indices, edge_t, submeshes_with_pt);

		edgeStitching.multiplied_edges_start.push_back(edge_constraints_cnt);
		// We got 'n' different vertex pairs, hence we need n-1 (circular) constraints
		double double_t = CGAL::to_double(edge_t);
		for (int const_i = 0; const_i < edge_v_indices.size()-1; const_i++) {
			edgeStitching.edge_const_1.push_back(edge_v_indices[const_i]);
			edgeStitching.edge_const_2.push_back(edge_v_indices[const_i+1]);
			edgeStitching.edge_coordinates.push_back(double_t);
			edgeStitching.edge_coordinates_precise.push_back(edge_t);
			constrained_pts_non_unique.push_back(pt);
			edge_constraints_cnt++;
		}
		// Save the number of duplicated edges
		edgeStitching.multiplied_edges_num.push_back(edge_v_indices.size()-1);
		// Add all the duplicated edge points
		for (int dup_edge_i = 0; dup_edge_i < edge_v_indices.size(); dup_edge_i++) {
			edgeStitching.edge_to_duplicates[edge_v_indices[dup_edge_i]] = pt_const_i;
		}
		int index_to_one_of_the_edge_pts = edge_constraints_cnt-1; // doesn't matter which one who pick, they are equal
		for (auto subm_i : submeshes_with_pt) {
			edgeStitching.submesh_to_edge_pt[subm_i].push_back(index_to_one_of_the_edge_pts); // push one of the vertices
		}
		
		pt_const_i++;
	}
	
	// set stitched_curves (once per edge, even if it is duplicated)
	edgeStitching.stitched_curves.resize(polylines.size());
	
	for (int i = 0; i < polylines.size(); i++) {
		std::vector<Point_2> pts; polyline_to_points(polylines[i],pts);
		edgeStitching.stitched_curves[i].resize(pts.size());
		for (int j = 0; j < pts.size(); j++) {
			Number_type edge_t; std::vector<Edge> edge_v_indices; std::vector<int> submeshes_with_pt;
			pt_to_edge_coordinates(pts[j], creasePattern, submeshVList, edge_v_indices, edge_t, submeshes_with_pt);

			auto edge_v = CGAL::to_double(edge_t)*V.row(edge_v_indices[0].v1)+(1-CGAL::to_double(edge_t))*V.row(edge_v_indices[0].v2);
			// Choose one of the edges, don't need the duplicates and they are all (approximately) equal anyhow
			edgeStitching.stitched_curves[i][j] = EdgePoint(edge_v_indices[0],CGAL::to_double(edge_t));
		}
	}
}

void save_submesh_bnd_edge_points(const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList,
				DogEdgeStitching& edgeStitching) {
	edgeStitching.submesh_to_bnd_edge.resize(submeshVList.size());

	std::vector<Point_2> boundary_vertices = creasePattern.boundary()->get_all_boundary_points();
	for (auto pt: boundary_vertices) {
		Number_type edge_t; std::vector<Edge> edge_v_indices; std::vector<int> submeshes_with_pt;
		pt_to_edge_coordinates(pt, creasePattern, submeshVList, edge_v_indices, edge_t, submeshes_with_pt);
		//for (auto subm_i : submeshes_with_pt) {
		for (int i = 0; i < submeshes_with_pt.size(); i++) {
			int subm_i = submeshes_with_pt[i];
			EdgePoint edgePoint(edge_v_indices[i], CGAL::to_double(edge_t));
			edgeStitching.submesh_to_bnd_edge[subm_i].push_back(edgePoint); // push one of the vertices
			//std::cout << "t = " << edge_t << std::endl;
		}
	}
	//int wait; std::cin >> wait;
}

void get_faces_partitions_to_submeshes(const CreasePattern& creasePattern, std::vector<bool>& is_polygon_hole, std::vector<SubmeshPoly>& submesh_polygons) {
	std::vector<Polygon_2> faces_polygons; std::vector<int> faces_to_submesh;

	// Get orth grid and add it the polylines
	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	PlanarArrangement grid_with_snapped(orthGrid);
	grid_with_snapped.add_polylines(creasePattern.get_clipped_fold_polylines()); 
	grid_with_snapped.add_polylines(creasePattern.get_clipped_bnd_polylines());
	grid_with_snapped.get_faces_polygons(faces_polygons);


	std::vector<Polygon_with_holes_2> submeshBnd;
	creasePattern.get_submeshes_faces_polygons(submeshBnd);
	//std::cout << "submeshBnd.size() = " << submeshBnd.size() << std::endl; int wait; std::cin >> wait;

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
			if (is_polygon_hole[submesh_i]) {submesh_i++; continue;}
			bool is_in_submesh = true;
			//std::cout << "checking if its in polygon = " << submeshBnd[submesh_i] << std::endl;
			for (auto vptr = face_polygon.vertices_begin(); vptr != face_polygon.vertices_end(); vptr++) {
				auto pt_in_face = pt_inside_polygon(submeshBnd[submesh_i],*vptr);
				//std::cout << "pt " << *vptr << " in face = " << pt_in_face << std::endl;
				is_in_submesh = is_in_submesh & (pt_in_face);
			}
			if (is_in_submesh) face_to_submesh[f_i] = submesh_i;
			submesh_i++;
		}

		//if (face_to_submesh[f_i] == -1) {std::cout << "Error, couldn't find a submesh for face = " << face_polygon << std::endl; exit(1); /*nothing to do but debug*/}
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
				//std::cout << "pt has key in submesh "  << std::endl;
				polygons_v_ren_indices[poly_i][poly_vi] = submesh_pt_to_V_ren[pt];
			} else if (poly_fold_pt_to_V_ren.count(pt)) {
				//std::cout << "pt has key in fold map" << std::endl;
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

void pt_to_edge_coordinates(const Point_2& pt, const CreasePattern& creasePattern, const std::vector<Eigen::MatrixXd>& submeshVList, 
				std::vector<Edge>& edge_v_indices, Number_type& edge_t, std::vector<int>& submeshes_with_pt) {

	const OrthogonalGrid& orthGrid(creasePattern.get_orthogonal_grid());
	std::pair<Point_2,Point_2> edge_pts; 
	if (!orthGrid.get_pt_edge_coordinates(pt, edge_pts,edge_t)) {
		std::cout << "Error, got a point pt = " << pt << " that is not on the grid " << std::endl;
		exit(1); // Should never get here, and if so all is lost
	}
	std::vector<Polygon_with_holes_2> submeshBnd; 
	creasePattern.get_submeshes_faces_polygons(submeshBnd);
	//std::cout << "point p = " << pt << " lies between " << edge_pts.first << " and " << edge_pts.second << " with t = " << t << std::endl;
	// Now find the indices of both points and add them as constraints
	// For every point, find all submeshes that contain it. We need to have both points for a submesh to count.
	
	Eigen::RowVector3d pt1(CGAL::to_double(edge_pts.first.x()),CGAL::to_double(edge_pts.first.y()),0);
	Eigen::RowVector3d pt2(CGAL::to_double(edge_pts.second.x()),CGAL::to_double(edge_pts.second.y()),0);
	
	// Go through every submesh
	int global_idx_base = 0;
	for (int sub_i = 0; sub_i < submeshVList.size(); sub_i++) {
		bool pt_in_submesh = pt_inside_polygon(submeshBnd[sub_i],pt);

		// Find the indices of the pts in the submesh (if they exist)
		int pt1_idx = -1, pt2_idx = -1;
		for (int v_i = 0; v_i < submeshVList[sub_i].rows(); v_i++) {
			if (submeshVList[sub_i].row(v_i) == pt1) {pt1_idx = v_i;}
			if (submeshVList[sub_i].row(v_i) == pt2) {pt2_idx = v_i;}
			if ((pt1_idx!=-1) && (pt2_idx != -1)) break;
		}
		// Found the edge (both points) on the submesh
		if ((pt1_idx!=-1)&& (pt2_idx != -1)) {
			Number_type min_edge = std::min(orthGrid.get_avg_x_edge(),orthGrid.get_avg_y_edge());
			double step_size = 0.01*CGAL::to_double(min_edge);

			Point_2 offset_xplus = Point_2(pt.x()+Number_type(step_size),pt.y());
			Point_2 offset_xminus = Point_2(pt.x()-Number_type(step_size),pt.y());
			Point_2 offset_yplus = Point_2(pt.x(),pt.y()+Number_type(step_size));
			Point_2 offset_yminus = Point_2(pt.x()-Number_type(step_size),pt.y()-+Number_type(step_size));
			
			bool pos_eps_insidex = pt_inside_polygon(submeshBnd[sub_i],offset_xplus);//(submeshBnd[sub_i].bounded_side(offset_xplus) == CGAL::ON_BOUNDED_SIDE);
			bool pos_m_eps_insidex = pt_inside_polygon(submeshBnd[sub_i],offset_xminus);//(submeshBnd[sub_i].bounded_side(offset_xminus) == CGAL::ON_BOUNDED_SIDE);
			bool pos_eps_insidey = pt_inside_polygon(submeshBnd[sub_i],offset_yplus);//(submeshBnd[sub_i].bounded_side(offset_yplus) == CGAL::ON_BOUNDED_SIDE);
			bool pos_m_eps_insidey = pt_inside_polygon(submeshBnd[sub_i],offset_yminus);//(submeshBnd[sub_i].bounded_side(offset_yminus) == CGAL::ON_BOUNDED_SIDE);

			pt_in_submesh = pt_in_submesh | pos_eps_insidex | pos_m_eps_insidex | pos_eps_insidey | pos_m_eps_insidey;

			//std::cout << "inside and pt_in_submesh = " << pt_in_submesh << std::endl;
			
			if (pt_in_submesh) {
				// Insert the (global V) indices
				edge_v_indices.push_back(Edge(global_idx_base+pt1_idx,global_idx_base+pt2_idx));
				submeshes_with_pt.push_back(sub_i);
			}
			
		}
		global_idx_base += submeshVList[sub_i].rows();
	}
	//std::cout << "global_edge_indices.size() = " << global_edge_indices.size() << std::endl;
	// We then need to map it to the global V which just have the submeshes concatenated
	/*
	if (!edge_v_indices.size()){
		std::cout << "Error, got an edge that is not in a submesh, with pt1 = " << pt1 << " and pt2 = " << pt2 << std::endl;
		exit(1); // Should not get here, and if so there's really nothing to do but debug the crease pattern
	}
	*/	
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

void generate_V_ren_list(Eigen::MatrixXd& V, std::vector<Eigen::MatrixXd>& submeshVList,DogEdgeStitching& eS, std::vector<Eigen::MatrixXd>& V_ren_list) {
	V_ren_list.resize(submeshVList.size());
	for (int i = 0; i < submeshVList.size(); i++) {
		V_ren_list[i].resize(submeshVList[i].rows() + eS.submesh_to_edge_pt[i].size()+eS.submesh_to_bnd_edge[i].size(),3);
		// add normal points
		for (int j = 0; j < submeshVList[i].rows(); j++) {V_ren_list[i].row(j) = submeshVList[i].row(j);}
		// add crease points
		for (int j = 0; j < eS.submesh_to_edge_pt[i].size(); j++) {
			int ei = eS.submesh_to_edge_pt[i][j];
			auto t = eS.edge_coordinates_precise[ei];
			auto coord_x =  t*V(eS.edge_const_1[ei].v1,0) + (1-t)*V(eS.edge_const_1[ei].v2,0);
			auto coord_y =  t*V(eS.edge_const_1[ei].v1,1) + (1-t)*V(eS.edge_const_1[ei].v2,1);
			V_ren_list[i].row(submeshVList[i].rows() + j) << CGAL::to_double(coord_x), CGAL::to_double(coord_y),0;
			//V_ren_list[i].row(submeshVList[i].rows() + j) = 
			//double t = eS.edge_coordinates[ei];
			//V_ren_list[i].row(submeshVList[i].rows() + j) = t*V.row(eS.edge_const_1[ei].v1) + (1-t)*V.row(eS.edge_const_1[ei].v2);
		}
		for (int j = 0; j < eS.submesh_to_bnd_edge[i].size(); j++) {
			EdgePoint ep = eS.submesh_to_bnd_edge[i][j];
			V_ren_list[i].row(submeshVList[i].rows() + eS.submesh_to_edge_pt[i].size() + j) = ep.getPositionInMesh(V);
			//std::cout << "Added row = " << V_ren_list[i].row(submeshVList[i].rows() + eS.submesh_to_edge_pt[i].size() + j)  << std::endl;
		}
	}
}


void generate_rendered_mesh_vertices_and_faces(const CreasePattern& creasePattern, 
		const std::vector<SubmeshPoly>& submesh_polygons,
		std::vector<Eigen::MatrixXd>& V_ren_list, DogEdgeStitching& edgeStitching,
		Eigen::MatrixXd& V_ren,
		Eigen::MatrixXi& F_ren) {
	typedef std::pair<double,double> PointDouble;
	// First create a list of F_ren
	std::vector<Eigen::MatrixXi> F_ren_list(V_ren_list.size());
	for (int subi = 0; subi < V_ren_list.size(); subi++) {
		std::vector<std::vector<int> > polygons_v_ren_indices; int poly_n = 0;
		// Create faces for every submesh
		// We need to convert the submeshpolygons indices from coordinates to the ones in the Vlist
		// We can search the closest to each one, but its better to just cache all points before
		std::map<PointDouble, int> point_to_sub_V_ren;
		// map normal points
		for (int i = 0; i < V_ren_list[subi].rows(); i++) {
			point_to_sub_V_ren[PointDouble(V_ren_list[subi](i,0),V_ren_list[subi](i,1))] = i;
		}
		// map points
		// Now go through the submesh polygons
		for (int poly_i = 0; poly_i < submesh_polygons.size(); poly_i++) {
			auto submesh_i = submesh_polygons[poly_i].first; auto submesh_poly = submesh_polygons[poly_i].second;
			// Update the interval [min_sub_i, max_sub_i) to the new submesh vertices indices in V_ren
			if (submesh_i != subi) continue;
			// every pt should have some cached indices to point_to_sub_V_ren
			polygons_v_ren_indices.push_back(std::vector<int>(submesh_poly.size())); int poly_vi = 0;
			for (auto vptr = submesh_poly.vertices_begin(); vptr != submesh_poly.vertices_end(); vptr++) {
				PointDouble pt(CGAL::to_double(vptr->x()),CGAL::to_double(vptr->y()));
				
				if (point_to_sub_V_ren.count(pt)) {
					polygons_v_ren_indices[poly_n][poly_vi] = point_to_sub_V_ren[pt];
					//std::cout << "Mapped point to polygon at: (" << pt.first << "," << pt.second << ")" << std::endl;
				} else {
					//std::cout << "Could not map point to polygon at: (" << pt.first << "," << pt.second << ")" << std::endl;
					// find the closest point to it on Vsub. Just a bit slower but there's some rounding error
					double min_dist = std::numeric_limits<double>::infinity();
					Eigen::RowVector3d pt_3d(pt.first,pt.second,0); int closest_idx = 0;
					for (int k = 0; k < V_ren_list[subi].rows(); k++) {
						if ((pt_3d-V_ren_list[subi].row(k)).norm() < min_dist) {
							min_dist = (pt_3d-V_ren_list[subi].row(k)).norm();
							closest_idx = k;
						}
					}
					//std::cout << "found point at idx = " << closest_idx << std::endl;
					polygons_v_ren_indices[poly_n][poly_vi] = closest_idx;
				}
				poly_vi++;
			}
			poly_n++;
		}
		// Create the submesh
		igl::polygon_mesh_to_triangle_mesh(polygons_v_ren_indices,F_ren_list[subi]);
	}
	igl::combine(V_ren_list,F_ren_list, V_ren, F_ren);
}

bool pt_inside_polygon(const Polygon_with_holes_2& poly, const Point_2& pt) {
	//if (poly.bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE) return false;
	if (poly.outer_boundary().bounded_side(pt) == CGAL::ON_UNBOUNDED_SIDE) return false;
	for (auto hole = poly.holes_begin(); hole != poly.holes_end(); hole++) if (hole->bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) return false;
	return true;
}

Number_type bbox_diff(const CGAL::Bbox_2& bbox1, const CGAL::Bbox_2& bbox2) {
	auto max_diff_x = max(CGAL::abs(bbox1.xmax()-bbox2.xmax()),CGAL::abs(bbox1.xmin()-bbox2.xmin()));
	auto max_diff_y = max(CGAL::abs(bbox1.ymax()-bbox2.ymax()),CGAL::abs(bbox1.ymin()-bbox2.ymin()));
	return max(max_diff_x, max_diff_y);
}

Number_type bbox_max_edge(const CGAL::Bbox_2& bbox) {
	return max(CGAL::abs(bbox.xmax()-bbox.xmin()),CGAL::abs(bbox.ymax()-bbox.ymin()));
}

std::vector<bool> submesh_is_hole(const CreasePattern& creasePattern) {
	auto holes = creasePattern.boundary()->get_holes();
	std::vector<Polygon_with_holes_2> facePolygons; creasePattern.get_submeshes_faces_polygons(facePolygons);	
	std::vector<bool> submesh_is_hole(facePolygons.size(), false); int poly_i = 0;
	for (auto poly : facePolygons) {
		auto bbox_error_threshold = Number_type(1e-2)*bbox_max_edge(poly.bbox());
		for (auto hole : holes) {
			if (bbox_diff(poly.outer_boundary().bbox(), hole.bbox()) < bbox_error_threshold) {
				submesh_is_hole[poly_i] = true;
			}
		}
		if (submesh_is_hole[poly_i]) std::cout << "found a hole" << std::endl;
		else std::cout << "not a hole" << std::endl;
		poly_i++;
	}
	return submesh_is_hole;
}
