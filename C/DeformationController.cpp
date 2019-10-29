#include "DeformationController.h"

#include <queue>
using namespace std;

DeformationController::DeformationController() : dogEditor(NULL), globalDog(NULL),
			 editedSubmesh(NULL), dogSolver(NULL), curveConstraintsBuilder(NULL) {
	// empty on purpose
}

void DeformationController::setup_optimization_measurements(std::string log_file_name) {
	opt_measurements_log = new ofstream(log_file_name);
}

void DeformationController::init_from_new_dog(Dog& dog) {
	//if (globalDog) delete globalDog;
	//if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = -1; // Editing the global dog

	dogEditor = new DogEditor(*viewer, *globalDog, edit_mode, select_mode,
								has_new_constraints,b,bc,paired_vertices,edgePoints,edgeCoords, z_only_editing);

	init_x0 = dog.getV_vector();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles,
				mvTangentCreaseAngleParams, mv_cos_angles, paired_vertices, bnd_vertices_pairs, opt_measurements_log);

	//std::cout << "setting up boundary curves!" << std::endl; dogSolver->getDog().setup_boundary_curves_indices();

	foldingDihedralAngleConstraintsBuilder = new FoldingDihedralAngleConstraintsBuilder(*globalDog, deformation_timestep);
	mvFoldingDihedralAngleConstraintsBuilder = new MVFoldingDihedralAngleConstraintsBuilder(*globalDog, deformation_timestep);
}

bool DeformationController::has_constraints() {
	// positional constraints
	if ( (b.rows() + edgePoints.size()) > 0) return true;
	if (is_curve_constraint) return true;
	if (dihedral_constrained.size()) return true;
	if (paired_vertices.size()) return true;
	return false;
}

void DeformationController::change_submesh(int submesh_i){
	if(submesh_i < -1 || submesh_i >= globalDog->get_submesh_n() ){
		//cycle through available submeshes
		submesh_i = editedSubmeshI + 1;
		if( submesh_i < -1 || submesh_i >= globalDog->get_submesh_n() ) submesh_i = -1;
	}
	editedSubmeshI = submesh_i;
	editedSubmesh = globalDog->get_submesh(submesh_i);
}

void DeformationController::single_optimization() {
	if ((is_time_dependent_deformation) && (deformation_timestep < 1) ) {
		deformation_timestep+=deformation_timestep_diff;
		p.pair_weight = deformation_timestep*0.1*p.soft_pos_weight;
		// the default of paired_boundary_bending_weight_mult is 1, so we just interpolate the weight from 0 to the bending weight
		p.paired_boundary_bending_weight = deformation_timestep*paired_boundary_bending_weight_mult*p.bending_weight;
		//p.paired_boundary_bending_weight = deformation_timestep*paired_boundary_bending_weight_mult*1;//p.bending_weight;
	}
	if (has_new_constraints) reset_dog_solver();
	if (is_curve_constraint) update_edge_curve_constraints();
	update_dihedral_constraints();
	dogSolver->set_solver_mode(solver_mode);
	dogSolver->update_point_coords(bc);
	dogSolver->update_edge_coords(edgeCoords);

	Eigen::VectorXd data_i;
	if(current_iteration < stored_iterations) data_i = dogSolver->get_obj_parts();
	dogSolver->single_iteration(constraints_deviation, objective);

	if(current_iteration < stored_iterations){
		obj_data.block(current_iteration,0, 1,3) << data_i[0], data_i[1], data_i[2];
		if(current_iteration+1 < stored_iterations){
			obj_data(current_iteration+1,3) = constraints_deviation;
			obj_data(current_iteration+1,4) = objective;
		}
	}
	if(current_iteration == stored_iterations - 1){
		cout << "Writing output\n";
		std::ofstream outfile;
		outfile.open ("output.txt");

		outfile << "Output: "<< stored_iterations << " iterations" << endl;

		outfile << "iteration" << endl;
		outfile << 0;
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << i;
		outfile << endl;

		outfile << "bending" << endl;
		outfile << obj_data(0,0);
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << obj_data(i,0);
		outfile << endl;

		outfile << "isometry" << endl;
		outfile << obj_data(0,1);
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << obj_data(i,1);
		outfile << endl;

		outfile << "objective (global mesh)" << endl;
		outfile << obj_data(0,2);
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << obj_data(i,2);
		outfile << endl;

		outfile << "constraints deviation (from iteration)" << endl;
		outfile << obj_data(0,3);
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << obj_data(i,3);
		outfile << endl;

		outfile << "objective (from iteration)" << endl;
		outfile << obj_data(0,4);
		for(int i=1; i<stored_iterations; ++i) outfile << ", " << obj_data(i,4);
		outfile << endl;

		outfile.close();
		cout << "Written to output.txt\n";
	}
	++current_iteration;

}

void DeformationController::apply_new_editor_constraint() {
	if (edit_mode == DogEditor::VERTEX_PAIRS) {
		if ( (dogEditor->pair_vertex_1!= -1) && (dogEditor->pair_vertex_2!= -1) ) {
			int vnum = globalDog->get_v_num();
			for (int i = 0; i < 3; i++) {
					paired_vertices.push_back(std::pair<int,int>(i*vnum+dogEditor->pair_vertex_1,i*vnum+dogEditor->pair_vertex_2));
			}
			is_time_dependent_deformation = true;
			/*
			if (z_only_editing) {
				paired_vertices.push_back(std::pair<int,int>(2*vnum+dogEditor->pair_vertex_1,2*vnum+dogEditor->pair_vertex_2));
			} else {
				for (int i = 0; i < 3; i++) {
					paired_vertices.push_back(std::pair<int,int>(i*vnum+dogEditor->pair_vertex_1,i*vnum+dogEditor->pair_vertex_2));
				}
			}
			*/
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
	bnd_vertices_pairs.clear();
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
	auto vars = dogSolver->get_opt_vars();
	if (dogSolver) delete dogSolver;
	cout << "resetting dog solver" << endl;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles,
		 mvTangentCreaseAngleParams, mv_cos_angles, paired_vertices, bnd_vertices_pairs, opt_measurements_log);
	//cout << "edge_cos_angles.size() = "<< edge_cos_angles.size() << endl;
	//int wait; cin >> wait;
	dogSolver->set_opt_vars(vars);
	has_new_constraints = false;
}

void DeformationController::set_cylindrical_boundary_constraints() {
	if ((!dogSolver->getDog().left_bnd.size()) || (!dogSolver->getDog().right_bnd.size()) ) dogSolver->getDog().setup_boundary_curves_indices();
	// Set Y curve constraints
	auto left_curve = dogSolver->getDog().left_bnd; auto right_curve = dogSolver->getDog().right_bnd;
	// Add curves pair constraints and curves boundary smoothness constraints
	for (int pair_i = 0; pair_i < left_curve.size(); pair_i++) {
		int vnum = globalDog->get_v_num();
		int v1(left_curve[pair_i]),v2(right_curve[pair_i]);
		for (int axis = 0; axis < 3; axis++) {
			paired_vertices.push_back(std::pair<int,int>(axis*vnum+v1,axis*vnum+v2));
		}
		bnd_vertices_pairs.push_back(std::pair<int,int>(v1,v2));

	}
	is_time_dependent_deformation = true;
	reset_dog_solver();
}

void DeformationController::store_data(int num_iterations){
	stored_iterations = num_iterations;
	obj_data = Eigen::MatrixXd::Zero(stored_iterations,5);
}
