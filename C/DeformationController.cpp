#include "DeformationController.h"

#include <queue>
using namespace std;

DeformationController::DeformationController() : dogEditor(NULL), globalDog(NULL),
			 editedSubmesh(NULL), dogSolver(NULL), curveConstraintsBuilder(NULL) {
	// empty on piurpose
}

void DeformationController::setup_optimization_measurements(std::string log_file_name) {
	opt_measurements_log = new ofstream(log_file_name);
}

void DeformationController::init_from_new_dog(Dog& dog) {
	if (globalDog) delete globalDog;
	if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = -1; // Editing the global dog

	dogEditor = new DogEditor(*viewer, *globalDog, edit_mode, select_mode, 
								has_new_constraints,b,bc,paired_vertices,edgePoints,edgeCoords, z_only_editing);
	
	init_x0 = dog.getV_vector();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles, 
				mvTangentCreaseAngleParams, mv_cos_angles, paired_vertices, matching_curve_pts_y, matching_curve_pts_x, opt_measurements_log);

	foldingDihedralAngleConstraintsBuilder = new FoldingDihedralAngleConstraintsBuilder(*globalDog, deformation_timestep);
	mvFoldingDihedralAngleConstraintsBuilder = new MVFoldingDihedralAngleConstraintsBuilder(*globalDog, deformation_timestep);
}

bool DeformationController::has_constraints() {
	// positional constraints
	if ( (b.rows() + edgePoints.size()) > 0) return true;
	if (is_curve_constraint) return true;
	if (dihedral_constrained.size()) return true;
	return false;
}

void DeformationController::single_optimization() {
	if ((is_time_dependent_deformation) && (deformation_timestep < 1) ) {
		deformation_timestep+=deformation_timestep_diff;
	}
	if (has_new_constraints) reset_dog_solver();
	if (is_curve_constraint) update_edge_curve_constraints();
	update_dihedral_constraints();
	dogSolver->update_point_coords(bc);
	dogSolver->update_edge_coords(edgeCoords);
	dogSolver->single_iteration(constraints_deviation, objective);
}

void DeformationController::apply_new_editor_constraint() {
	if (edit_mode == DogEditor::VERTEX_PAIRS) {
		if ( (dogEditor->pair_vertex_1!= -1) && (dogEditor->pair_vertex_2!= -1) ) {
			int vnum = globalDog->get_v_num();
			for (int i = 0; i < 3; i++) {
				paired_vertices.push_back(std::pair<int,int>(i*vnum+dogEditor->pair_vertex_1,i*vnum+dogEditor->pair_vertex_2));	
			}
		}
	} else if (edit_mode == DogEditor::DIHEDRAL_ANGLE) {
		if (dogEditor->picked_edge.t !=-1) {
			foldingDihedralAngleConstraintsBuilder->add_constraint(dogEditor->picked_edge, dst_dihedral_angle);
			foldingDihedralAngleConstraintsBuilder->get_edge_angle_pairs(edge_angle_pairs);
			foldingDihedralAngleConstraintsBuilder->get_edge_angle_constraints(edge_cos_angles);
			dihedral_constrained.push_back(dogEditor->picked_edge);
			is_time_dependent_deformation = true;
		}
	} else if (edit_mode == DogEditor::MV_DIHEDRAL_ANGLE) {
		if (dogEditor->picked_edge.t !=-1) {
			mvFoldingDihedralAngleConstraintsBuilder->add_constraint(dogEditor->picked_edge, dst_dihedral_angle);
			mvFoldingDihedralAngleConstraintsBuilder->get_mv_tangent_crease_folds(mvTangentCreaseAngleParams);
			mvFoldingDihedralAngleConstraintsBuilder->get_edge_angle_constraints(mv_cos_angles);
			//dihedral_constrained.push_back(dogEditor->picked_edge);
			is_time_dependent_deformation = true;
		}
	} else if (edit_mode == DogEditor::EDGES_ANGLE) {
		if ((dogEditor->edge_angle_v1 != -1) && (dogEditor->edge_angle_v2 != -1) && (dogEditor->edge_angle_center != -1) ) {
			cout << "adding constraint with edge = " << dogEditor->edge_angle_v1 << "," << dogEditor->edge_angle_center << 
						"," << dogEditor->edge_angle_v2 << endl;
			edge_angle_pairs.push_back(pair<Edge,Edge>(Edge(dogEditor->edge_angle_v1,dogEditor->edge_angle_center),
														Edge(dogEditor->edge_angle_center,dogEditor->edge_angle_v2)));
			edge_cos_angles.push_back(cos(dst_dihedral_angle));
		}
	}
	has_new_constraints = true;
	reset_new_editor_constraint();
}

void DeformationController::setup_curve_constraints() {
	if (curveConstraintsBuilder) delete curveConstraintsBuilder;
	if (globalDog->has_creases()) {
		curveConstraintsBuilder = new CurveInterpolationConstraintsBuilder(globalDog->getV(), 
															globalDog->getEdgeStitching(), deformed_curve_idx, deformation_timestep,
															curve_k_translation, curve_k_mult, curve_t_addition, max_curve_points);
	} else {
		std::vector<int> v_indices; int n = globalDog->getV().rows();
		for (int i = 0; i < sqrt(n); i++) v_indices.push_back(i);
		std::cout << "v_indices.size() = " << v_indices.size() << std::endl;
		curveConstraintsBuilder = new CurveInterpolationConstraintsBuilder(globalDog->getV(), 
															v_indices, deformation_timestep,
															curve_k_translation, curve_k_mult, curve_t_addition);
	}
	SurfaceCurve surfaceCurve; Eigen::MatrixXd edgeCoords;
	curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	is_curve_constraint = true;
	is_time_dependent_deformation = true;
	add_edge_point_constraints(surfaceCurve.edgePoints,edgeCoords);
}

void DeformationController::update_edge_curve_constraints() {
	if (curveConstraintsBuilder) {
		SurfaceCurve surfaceCurve;
		curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);	
	}
}

void DeformationController::update_dihedral_constraints() {
	foldingDihedralAngleConstraintsBuilder->get_edge_angle_constraints(edge_cos_angles);
	dogSolver->update_edge_angles(edge_cos_angles);

	mvFoldingDihedralAngleConstraintsBuilder->get_edge_angle_constraints(mv_cos_angles);
	dogSolver->update_mv_cos_angles(mv_cos_angles);
}

void DeformationController::update_edge_coords(Eigen::MatrixXd& edgeCoords_i) {
	edgeCoords = edgeCoords_i; dogSolver->update_edge_coords(edgeCoords);
}

void DeformationController::update_point_coords(Eigen::VectorXd& bc_i) {
	bc = bc_i; dogSolver->update_point_coords(bc);
}

void DeformationController::add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords) {
	edgePoints.insert(edgePoints.end(), new_edgePoints.begin(), new_edgePoints.end());
	Eigen::MatrixXd old_edgeCoords = edgeCoords; edgeCoords.resize(old_edgeCoords.rows()+new_edgeCoords.rows(), old_edgeCoords.cols());
	if (old_edgeCoords.rows()) edgeCoords << old_edgeCoords, new_edgeCoords; else edgeCoords = new_edgeCoords; // Eigen's concatenate crashes if one of them is empty

	has_new_constraints = true;
}

void DeformationController::add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc) {
	Eigen::VectorXi old_b = b; Eigen::VectorXd old_bc = bc;
	b.resize(old_b.rows()+new_b.rows()); bc.resize(b.rows());
	if (old_b.rows()) {
		b << old_b,new_b; bc << old_bc,new_bc;
	} else {
		b = new_b; bc = new_bc; // Eigen's concatenate crashes if one of them is empty
	}
	
	has_new_constraints = true;
}

void DeformationController::add_edge_point_constraint(const EdgePoint& new_edgePoint, const Eigen::RowVector3d& new_edgeCoords) {
	std::vector<EdgePoint> new_edgePoints = {new_edgePoint}; Eigen::MatrixXd newCoords(1,3); newCoords.row(0) = new_edgeCoords;
	add_edge_point_constraints(new_edgePoints, newCoords);
}

void DeformationController::add_pair_vertices_constraints(const std::vector<std::pair<int,int>>& new_pair_vertices) {
	paired_vertices.insert(paired_vertices.end(), new_pair_vertices.begin(), new_pair_vertices.end());
	has_new_constraints = true;
}

void DeformationController::add_pair_vertices_constraint(int v1, int v2) {
	std::vector<std::pair<int,int>> new_pair_vertices{std::pair<int,int>(v1,v2)}; 
	add_pair_vertices_constraints(new_pair_vertices);
}


void DeformationController::reset_constraints() {
	b.resize(0);
	bc.resize(0); 
	paired_vertices.clear(); 
	edgePoints.clear();
	dihedral_constrained.clear();
	edge_angle_pairs.clear(); 
	edge_cos_angles.clear();
	mv_cos_angles.clear();
	mvTangentCreaseAngleParams.clear();
	edgeCoords.resize(0,3); 
	dogEditor->clearHandles(); 
	reset_dog_solver(); 
	is_curve_constraint = false;
}

// t = 0.5 in the edge constraint means it is equally spaced
EdgePoint DeformationController::find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int &min_edge) {
	auto eS = globalDog->getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
	int curve_v_n = foldingCurve.size();
	min_edge = 1; double min_dist_from_equal = abs(0.5-foldingCurve[1].t);
	for (int ei = 1; ei < foldingCurve.size()-1; ei++) {
		double dist_from_equal = abs(0.5-foldingCurve[ei].t);
		if ( dist_from_equal < min_dist_from_equal) {
			min_edge = ei;
			min_dist_from_equal = dist_from_equal;
		}
	}
	return foldingCurve[min_edge];
}

void DeformationController::reset_dog_solver() {
	Dog& dog = dogSolver->getDog();
	if (dogSolver) delete dogSolver;
	cout << "reseting dog solver" << endl;
	std::cout << "matching_curve_pts_y.size() = " << matching_curve_pts_y.first.size() << std::endl;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles,
		 mvTangentCreaseAngleParams, mv_cos_angles, paired_vertices, matching_curve_pts_y, matching_curve_pts_x, opt_measurements_log);
	//cout << "edge_cos_angles.size() = "<< edge_cos_angles.size() << endl;
	//int wait; cin >> wait;
	has_new_constraints = false;
}

void DeformationController::set_wallpaper_constraints() {
	// Set Y curve constraints
	auto left_curve = dogSolver->getDog().left_bnd; auto right_curve = dogSolver->getDog().right_bnd;
	if ((wallpaperType == XUY) || (wallpaperType == UXUY) ) {std::reverse(right_curve.begin(),right_curve.end());}
	auto lower_curve = dogSolver->getDog().lower_bnd; auto upper_curve = dogSolver->getDog().upper_bnd;
	if ((wallpaperType == UXY) || (wallpaperType == UXUY) ) {std::reverse(upper_curve.begin(),upper_curve.end());}

	//matching_curve_pts_y.first = left_curve; matching_curve_pts_y.second = right_curve;
	matching_curve_pts_x.first = lower_curve; matching_curve_pts_x.second = upper_curve;

	reset_dog_solver();
}