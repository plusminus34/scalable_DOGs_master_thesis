#include "DeformationController.h"

#include <queue>
using namespace std;

void DeformationController::init_from_new_dog(Dog& dog) {
	if (globalDog) delete globalDog;
	if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = -1; // Editing the global dog
	dogEditor.init_from_new_dog(dog);
}

void DeformationController::update_fold_constraints() {
	Eigen::VectorXd bc_MV; Eigen::MatrixXd edgeCoords;
	mvFoldingConstraintsBuilder.get_folds_constraint_coords(fold_dihedral_angle, *globalDog, bc_MV, edgeCoords);
	Eigen::VectorXd bc_ref;
	refFoldingConstrainsBuilder.get_folds_constraint_coords(*globalDog, bc_ref);

	Eigen::VectorXd bc(bc_MV.rows()+bc_ref.rows());
	int cnt = 0;
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < bc_MV.rows()/3; i++) {
			bc(cnt) = bc_MV(coord_i*bc_MV.rows()/3+i);
			cnt++;
		}
		for (int i = 0; i < bc_ref.rows()/3; i++) {
			bc(cnt) = bc_ref(coord_i*bc_ref.rows()/3+i);
			cnt++;
		}
	}
	dogEditor.update_edge_coords(edgeCoords);
	dogEditor.update_point_coords(bc);
}
void DeformationController::single_optimization() {
	if (is_folding()) update_fold_constraints();
	if (is_curve_constraint) update_edge_curve_constraints();
	return dogEditor.single_optimization();
}

void DeformationController::setup_curve_constraints() {
	if (curveConstraintsBuilder) delete curveConstraintsBuilder;
	curveConstraintsBuilder = new CurveInterpolationConstraintsBuilder(globalDog->getV(), 
															globalDog->getEdgeStitching(), curve_timestep);
	SurfaceCurve surfaceCurve; Eigen::MatrixXd edgeCoords;
	curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	is_curve_constraint = true;
	dogEditor.add_edge_point_constraints(surfaceCurve.edgePoints,edgeCoords);
}

void DeformationController::update_edge_curve_constraints() {
	SurfaceCurve surfaceCurve; Eigen::MatrixXd edgeCoords;
	curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	dogEditor.update_edge_coords(edgeCoords);
}

void DeformationController::setup_reflection_fold_constraints() {
	cout << "Setting up fold constraints!" << endl;

	// Setting up on M/V fold
	int vnum = globalDog->getV().rows();
	auto eS = globalDog->getEdgeStitching(); auto curves = eS.stitched_curves;
	int fold_curve_idx = 0; int e_idx; 
	EdgePoint edgePoint = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx);
	e_idx=curves[0].size()/2;

	bool is_mountain = true; bool keep_rigid_motion = true;
	mvFoldingConstraintsBuilder.add_fold(*globalDog, fold_curve_idx, e_idx, is_mountain, keep_rigid_motion);

	int submesh_n = globalDog->get_submesh_n(); int curves_n = curves.size();
	vector<bool> passed_on_submesh(submesh_n, false); vector<bool> passed_on_curves(curves_n, false);
	passed_on_curves[0] = true;
	
	// 0) Implement a method called submesh to curves. Goes through every curve points, then get some edge points and get all submeshes related to it - Done
	auto subm_to_curves = submeshes_to_curves(*globalDog);
	// 1) Mark one of the 2 submeshes as "passed"
	int sub1_v,sub2_v; globalDog->get_2_inner_vertices_from_edge(curves[0][e_idx].edge,sub1_v,sub2_v);
	int submesh1 =  globalDog->v_to_submesh_idx(sub1_v); int submesh2 =  globalDog->v_to_submesh_idx(sub2_v); 
	passed_on_submesh[submesh1] = true; passed_on_submesh[submesh2] = true; 
	int cur_submesh = submesh1;
	//passed_on_submesh[globalDog->v_to_submesh_idx(curves[0][e_idx].edge.v2)] = true;
	std::cout << "Passed on sub mesh " << cur_submesh << std::endl;
	auto adjacency_list = globalDog->get_submesh_adjacency();

	// 2) Go through all the neighbouring submeshes
	//queue<int> Q; for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) Q.push(adjacency_list[cur_submesh][i]);
	queue<int> Q; Q.push(submesh1); Q.push(submesh2);
	while (!Q.empty()) {
		cur_submesh = Q.front(); Q.pop();
		//std::cout << "cur submesh = " << cur_submesh << std::endl;
		// Go through all of the stitched curves around the submesh
		for (auto c_i: subm_to_curves[cur_submesh]) {
			// 3) Handle the stitched curves connected to the submesh if they were not passed on, then mark the relevant mesh as "passed", and same for the curves
			if (!passed_on_curves[c_i]) {
				int e_idx = curves[c_i].size()/2;
				//std::cout << "Curve " << c_i << " with submesh " << cur_submesh << std::endl;
				refFoldingConstrainsBuilder.add_fold(*globalDog, c_i, e_idx, passed_on_submesh);
				passed_on_curves[c_i] = true;
			}
		}

		passed_on_submesh[cur_submesh] = true;
		// Add all submeshes that were not processed
		for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) {
			int nb_submesh = adjacency_list[cur_submesh][i];
			if (!passed_on_submesh[nb_submesh]) Q.push(nb_submesh);
		}
	}

	Eigen::VectorXi b_MV; Eigen::VectorXd bc_MV; std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
	mvFoldingConstraintsBuilder.get_folds_constraint_indices(*globalDog, b_MV, edgePoints);
	mvFoldingConstraintsBuilder.get_folds_constraint_coords(fold_dihedral_angle, *globalDog, bc_MV, edgeCoords);

	Eigen::VectorXi b_ref; Eigen::VectorXd bc_ref;
	refFoldingConstrainsBuilder.get_folds_constraint_indices(*globalDog, b_ref);
	refFoldingConstrainsBuilder.get_folds_constraint_coords(*globalDog, bc_ref);

	Eigen::VectorXi b(b_MV.rows()+b_ref.rows()); Eigen::VectorXd bc(b_MV.rows()+b_ref.rows());
	int cnt = 0;
	for (int coord_i = 0; coord_i < 3; coord_i++) {
		for (int i = 0; i < b_MV.rows()/3; i++) {
			b(cnt) = b_MV(coord_i*b_MV.rows()/3+i);
			bc(cnt) = bc_MV(coord_i*b_MV.rows()/3+i);
			cnt++;
		}
		for (int i = 0; i < b_ref.rows()/3; i++) {
			b(cnt) = b_ref(coord_i*b_ref.rows()/3+i);
			bc(cnt) = bc_ref(coord_i*b_ref.rows()/3+i);
			cnt++;
		}
	}
	dogEditor.add_edge_point_constraints(edgePoints, edgeCoords);
	dogEditor.add_positional_constraints(b, bc);
}

void DeformationController::get_curve_fold_bias_obj() {
	bool dbg = true;
	CurvedFoldingBiasObjective tmpCurveFoldBiasObj(dbg);
	double alpha = 0.1; CurvedFoldingBiasSignObjective tmpCurveSignBiasSignObj(alpha,dbg);
	CurvedFoldBias curvedFoldBias;
	auto eS = globalDog->getEdgeStitching();
	
	int e_idx;	
	for (int fold_curve_idx = 0; fold_curve_idx < eS.stitched_curves.size(); fold_curve_idx++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
		auto edge_pt = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx); // currently not used

		e_idx = foldingCurve.size()/2;
		edge_pt = foldingCurve[e_idx];
		curvedFoldBias.ep_b = foldingCurve[e_idx-1]; curvedFoldBias.ep_f = foldingCurve[e_idx+1];
		curvedFoldBias.edge_t = edge_pt.t;
		globalDog->get_2_submeshes_vertices_from_edge(edge_pt.edge, curvedFoldBias.v1,curvedFoldBias.v2,curvedFoldBias.w1,curvedFoldBias.w2);
		tmpCurveFoldBiasObj.add_fold_bias(curvedFoldBias);
		tmpCurveSignBiasSignObj.add_fold_bias(curvedFoldBias);
	}
	double curve_fold_bias_obj = tmpCurveFoldBiasObj.obj(globalDog->getV_vector());
	std::cout << "Curve fold bias obj = " << curve_fold_bias_obj << std::endl;
	double curve_fold_bias_sign_obj = tmpCurveSignBiasSignObj.obj(globalDog->getV_vector());
	std::cout << "Curve fold bias sign obj approx = " << curve_fold_bias_sign_obj << std::endl;
}

void DeformationController::setup_fold_bias() {
	/*
	CurvedFoldBias curvedFoldBias;
	auto eS = globalDog->getEdgeStitching();
	
	int e_idx;	
	for (int fold_curve_idx = 0; fold_curve_idx < eS.stitched_curves.size(); fold_curve_idx++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
		auto edge_pt = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx); // currently not used

		e_idx = foldingCurve.size()/2;
		edge_pt = foldingCurve[e_idx];
		curvedFoldBias.ep_b = foldingCurve[e_idx-1]; curvedFoldBias.ep_f = foldingCurve[e_idx+1];
		curvedFoldBias.edge_t = edge_pt.t;
		globalDog->get_2_submeshes_vertices_from_edge(edge_pt.edge, curvedFoldBias.v1,curvedFoldBias.v2,curvedFoldBias.w1,curvedFoldBias.w2);
		curvedFoldingBiasObjective.add_fold_bias(curvedFoldBias);
	}
	Eigen::VectorXi b; Eigen::VectorXd bc;
	dogEditor.add_positional_constraints(b, bc);
	*/
	
	int vnum = globalDog->getV().rows();
	auto eS = globalDog->getEdgeStitching(); auto curves = eS.stitched_curves;
	int e_idx;
	int submesh_n = globalDog->get_submesh_n(); int curves_n = curves.size();
	vector<bool> passed_on_submesh(submesh_n, false); vector<bool> passed_on_curves(curves_n, false);
	
	
	// 0) Implement a method called submesh to curves. Goes through every curve points, then get some edge points and get all submeshes related to it - Done
	auto subm_to_curves = submeshes_to_curves(*globalDog);
	// 1) Mark one of the 2 submeshes as "passed"
	int sub1_v,sub2_v; globalDog->get_2_inner_vertices_from_edge(curves[0][e_idx].edge,sub1_v,sub2_v);
	int submesh1 =  globalDog->v_to_submesh_idx(sub1_v); int submesh2 =  globalDog->v_to_submesh_idx(sub2_v); 
	passed_on_submesh[submesh1] = true; passed_on_submesh[submesh2] = false; 
	int cur_submesh = submesh1;
	//passed_on_submesh[globalDog->v_to_submesh_idx(curves[0][e_idx].edge.v2)] = true;
	std::cout << "Passed on sub mesh " << cur_submesh << std::endl;
	auto adjacency_list = globalDog->get_submesh_adjacency();

	// 2) Go through all the neighbouring submeshes
	//queue<int> Q; for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) Q.push(adjacency_list[cur_submesh][i]);
	queue<int> Q; Q.push(submesh2);
	while (!Q.empty()) {
		cur_submesh = Q.front(); Q.pop();
		//std::cout << "cur submesh = " << cur_submesh << std::endl;
		// Go through all of the stitched curves around the submesh
		for (auto c_i: subm_to_curves[cur_submesh]) {
			// 3) Handle the stitched curves connected to the submesh if they were not passed on, then mark the relevant mesh as "passed", and same for the curves
			if (!passed_on_curves[c_i]) {
				int e_idx = curves[c_i].size()/2;
				//std::cout << "Curve " << c_i << " with submesh " << cur_submesh << std::endl;
				refFoldingConstrainsBuilder.add_fold(*globalDog, c_i, e_idx, passed_on_submesh);
				passed_on_curves[c_i] = true;
			}
			passed_on_submesh[cur_submesh] = true;
		}

		
		// Add all submeshes that were not processed
		for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) {
			int nb_submesh = adjacency_list[cur_submesh][i];
			if (!passed_on_submesh[nb_submesh]) Q.push(nb_submesh);
		}
	}

	Eigen::VectorXi b_ref; Eigen::VectorXd bc_ref;
	refFoldingConstrainsBuilder.get_folds_constraint_indices(*globalDog, b_ref);
	refFoldingConstrainsBuilder.get_folds_constraint_coords(*globalDog, bc_ref);

	dogEditor.add_positional_constraints(b_ref, bc_ref);


}

void DeformationController::setup_fold_constraints() {
	cout << "Setting up fold constraints!" << endl;
	int vnum = globalDog->getV().rows();
	int fold_curve_idx = 0; int e_idx; 
	EdgePoint edgePoint = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx);
	e_idx=globalDog->getEdgeStitching().stitched_curves[0].size()/2;

	bool is_mountain = true; bool keep_rigid_motion = true;
	mvFoldingConstraintsBuilder.add_fold(*globalDog, fold_curve_idx, e_idx, is_mountain, keep_rigid_motion);


	// Add another fold constraint
	fold_curve_idx = 1;
	CurvedFoldBias curvedFoldBias;
	auto eS = globalDog->getEdgeStitching(); is_mountain = true;
	if (eS.stitched_curves.size() > 1 ){
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
		for (fold_curve_idx = 1; fold_curve_idx < eS.stitched_curves.size(); fold_curve_idx++) {
			auto edge_pt = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx); // currently not used

			e_idx = foldingCurve.size()/2;
			edge_pt = foldingCurve[e_idx];
			curvedFoldBias.ep_b = foldingCurve[e_idx-1]; curvedFoldBias.ep_f = foldingCurve[e_idx+1];
			curvedFoldBias.edge_t = edge_pt.t;
			globalDog->get_2_submeshes_vertices_from_edge(edge_pt.edge, curvedFoldBias.v1,curvedFoldBias.v2,curvedFoldBias.w1,curvedFoldBias.w2);
			curvedFoldingBiasObjective.add_fold_bias(curvedFoldBias);

			//edgePoint = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx); keep_rigid_motion = false;
			//mvFoldingConstraintsBuilder.add_fold(*globalDog, fold_curve_idx, e_idx, is_mountain, keep_rigid_motion);

			//is_mountain = !is_mountain;
		}
	}
	
	//mvFoldingConstraintsBuilder
	/*
	
	*/
	
	Eigen::VectorXi b; Eigen::VectorXd bc; std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
	mvFoldingConstraintsBuilder.get_folds_constraint_indices(*globalDog, b, edgePoints);
	mvFoldingConstraintsBuilder.get_folds_constraint_coords(fold_dihedral_angle, *globalDog, bc, edgeCoords);

	dogEditor.add_edge_point_constraints(edgePoints, edgeCoords);
	dogEditor.add_positional_constraints(b, bc);
}

void DeformationController::setup_fold_constraints_old() {
	cout << "Setting up fold constraints!" << endl;
	int vnum = globalDog->getV().rows();

	// For now only handle the case of 1-crease

	// Find an edge on a polyline (choose the one with the most equal distance from both sides)
	int fold_curve_idx = 0; int e_idx; 
	EdgePoint edgePoint = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx, e_idx);
	cout << "Found edge with t = " << edgePoint.t << endl;
	int v1,v2; globalDog->get_2_inner_vertices_from_edge(edgePoint.edge,v1,v2);
	std::cout << "inner v1 = " << v1 << " inner_v2 = " << v2 << std::endl;
	// Add a pair-closeness-constraint between 2 inner vertices of different submeshes edge points
	for (int i = 0; i < 3; i++) dogEditor.add_pair_vertices_constraint(i*vnum+v1,i*vnum+v2);
	// Add an edge point constraint to an higher Z value
	auto z_offset = 0.1;
	Eigen::RowVector3d new_edge_pt_pos = edgePoint.getPositionInMesh(globalDog->getV()) + Eigen::RowVector3d(0,0,z_offset);
	dogEditor.add_edge_point_constraint(edgePoint, new_edge_pt_pos);

	Eigen::VectorXi pos_const_i(2); Eigen::VectorXd pos_const(2);
	pos_const_i << v1+2*vnum,v2+2*vnum; pos_const << 0,0;
	dogEditor.add_positional_constraints(pos_const_i, pos_const);

	/*
	auto eS = globalDog->getEdgeStitching();
	for (; fold_curve_idx < eS.stitched_curves.size(); fold_curve_idx++) {
		edgePoint = find_most_equally_spaced_edge_on_fold_curve(fold_curve_idx);
		globalDog->get_2_inner_vertices_from_edge(edgePoint.edge,v1,v2);
		for (int i = 0; i < 3; i++) dogEditor.add_pair_vertices_constraint(i*vnum+v1,i*vnum+v2);
	}
	*/
}

// t = 0.5 in the edge constraint means it is equally spaced
EdgePoint DeformationController::find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int &min_edge) {
	auto eS = globalDog->getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
	int curve_v_n = foldingCurve.size();
	min_edge = floor(curve_v_n/7); double min_dist_from_equal = abs(0.5-foldingCurve[floor(curve_v_n/7)].t);
	for (int ei = floor(curve_v_n/7)+1; ei < foldingCurve.size()-floor(curve_v_n/7); ei++) {
		double dist_from_equal = abs(0.5-foldingCurve[ei].t);
		if ( dist_from_equal < min_dist_from_equal) {
			min_edge = ei;
			min_dist_from_equal = dist_from_equal;
		}
	}
	return foldingCurve[min_edge];
}

void DeformationController::update_edited_mesh(int newEditedSubmeshI) {
	if (globalDog->get_submesh_n() == 1) return; // Only 1 submesh (no folds)
	if (newEditedSubmeshI == editedSubmeshI) return; // no change
	if ( (newEditedSubmeshI >= globalDog->get_submesh_n()) || (newEditedSubmeshI < -1 ) ) return;

	if (newEditedSubmeshI != editedSubmeshI) {
		if (editedSubmesh && (editedSubmesh!= globalDog)) delete editedSubmesh;
		editedSubmeshI = newEditedSubmeshI;
		if (newEditedSubmeshI == -1) {
			editedSubmesh = globalDog;
		} else {
			editedSubmesh = globalDog->get_submesh(editedSubmeshI);	
		}
		
		dogEditor.init_from_new_dog(*editedSubmesh);
	}
}

void DeformationController::propagate_submesh_constraints() {
	int submesh_n = globalDog->get_submesh_n();
	if ( (editedSubmeshI < 0) || (editedSubmeshI >= submesh_n) ) return;

	auto adjacency_list = globalDog->get_submesh_adjacency();
	auto eS = globalDog->getEdgeStitching();
	vector<bool> edge_constraint_set(eS.edge_const_1.size(),false); vector<Eigen::RowVector3d> const_value(eS.edge_const_1.size());

	vector<bool> passed_on_submesh(submesh_n, false);
	passed_on_submesh[editedSubmeshI] = true; globalDog->update_submesh_V(editedSubmeshI, editedSubmesh->getV());
	update_edge_constraints_from_submesh(editedSubmeshI, eS, edge_constraint_set, const_value);
	queue<int> Q; for (int i = 0; i < adjacency_list[editedSubmeshI].size(); i++) Q.push(adjacency_list[editedSubmeshI][i]);
	while (!Q.empty()) {
		int cur_submesh = Q.front(); Q.pop();
		//std::cout << "processing submesh " << cur_submesh << std::endl;
		deform_submesh_based_on_previous_submeshes(cur_submesh, eS, edge_constraint_set, const_value);
		update_edge_constraints_from_submesh(cur_submesh, eS, edge_constraint_set, const_value);

		passed_on_submesh[cur_submesh] = true;
		// Add all submeshes that were not processed
		for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) {
			int nb_submesh = adjacency_list[cur_submesh][i];
			if (!passed_on_submesh[nb_submesh]) Q.push(nb_submesh);
		}
	}
	editedSubmeshI = -1;
	editedSubmesh = globalDog;
	dogEditor.init_from_new_dog(*editedSubmesh);
}

void DeformationController::update_edge_constraints_from_submesh(int submesh_i, const DogEdgeStitching& eS, 
									std::vector<bool>& edge_constraint_set, std::vector<Eigen::RowVector3d>& const_value) {
	for (int const_i = 0; const_i < eS.edge_const_1.size(); const_i++) {
		if (!edge_constraint_set[const_i]) {
			double t = eS.edge_coordinates[const_i]; Edge edge1(eS.edge_const_1[const_i]),edge2(eS.edge_const_2[const_i]);
			// Make sure that the constraint submesh edge, if exists, is at edge1
			if (globalDog->v_to_submesh_idx(edge2.v1) == submesh_i) {std::swap(edge1,edge2);}
			// Check if this constraint involves this submesh
			if ( (globalDog->v_to_submesh_idx(edge1.v1) == submesh_i) ) {
				// This constraint involves this submesh. Get all duplicated edges
				int mult_edge_index = eS.edge_to_duplicates.at(edge1);
				int mult_edge_const_start = eS.multiplied_edges_start[mult_edge_index]; 
				int mult_edges_num = eS.multiplied_edges_num[mult_edge_index];
				Eigen::RowVector3d const_pos = EdgePoint(edge1, t).getPositionInMesh(globalDog->getV());

				// Go through all duplicated edges
				for (int set_const_i = mult_edge_const_start; set_const_i < mult_edge_const_start+mult_edges_num; set_const_i++) {
					edge_constraint_set[set_const_i] = true;
					const_value[set_const_i] = const_pos;
				}
				// set the loop index to the end of the duplicated edges, so that it will process new edges constraints
				const_i = mult_edge_const_start + mult_edges_num -1;
			}
		}
	}
}

void DeformationController::deform_submesh_based_on_previous_submeshes(int submesh_i, const DogEdgeStitching& eS,
					 const std::vector<bool>& edge_constraint_set, const std::vector<Eigen::RowVector3d>& const_value) {
	auto submeshDog = globalDog->get_submesh(submesh_i);

	// process submesh_i
	int submesh_v_min_i, submesh_v_max_i;
	globalDog->get_submesh_min_max_i(submesh_i, submesh_v_min_i, submesh_v_max_i, true);
	
	std::vector<EdgePoint> edgePoints; std::vector<Eigen::RowVector3d> edgePointsCoordsList;

	// take all edge constraints that involve this submesh, and are already "set" by other submeshes
	for (int edge_const_i = 0; edge_const_i < eS.edge_const_1.size(); edge_const_i++) {
			double t = eS.edge_coordinates[edge_const_i]; Edge edge_src(eS.edge_const_1[edge_const_i]),edge_target(eS.edge_const_2[edge_const_i]);

			// This makes sure that if this constraint involves the current submesh then the source edge concerns the current submesh
			if (globalDog->v_to_submesh_idx(edge_target.v1) == submesh_i) {std::swap(edge_src,edge_target);}

			// Now we want to make sure this constraint does involve this submesh, and that it was already set by other submeshes
			int src_submesh = globalDog->v_to_submesh_idx(edge_src.v1);
			if ((src_submesh == submesh_i) && (edge_constraint_set[edge_const_i]) ) {
				// Add this as edge point constraints
				// needs to change from global V index to local v index
				Edge submeshEdge(edge_src); submeshEdge.v1 -= submesh_v_min_i; submeshEdge.v2 -= submesh_v_min_i;
				edgePoints.push_back(EdgePoint(submeshEdge,t));
				edgePointsCoordsList.push_back(const_value[edge_const_i]);

				// At this point we need to skip the rest of the duplicated constraints that might involve edge_const_i
				int mult_edge_index = eS.edge_to_duplicates.at(edge_src);
				int mult_edge_const_start = eS.multiplied_edges_start[mult_edge_index]; 
				int mult_edges_num = eS.multiplied_edges_num[mult_edge_index];
				edge_const_i = mult_edge_const_start + mult_edges_num -1;

			}
	}
	// Convert it to a matrix
	Eigen::MatrixXd edgePointCoords(edgePointsCoordsList.size(),3); 
	for (int edgePtRow = 0; edgePtRow < edgePointCoords.rows(); edgePtRow++) edgePointCoords.row(edgePtRow) = edgePointsCoordsList[edgePtRow];
	EdgePointConstraints submeshEdgePtConst(edgePoints, edgePointCoords);
	
	// Use newton penalty optimization while gradually pushing both positional constraints and dog constraints
	auto init_x0 = submeshDog->getV_vector();
	SimplifiedBendingObjective bending(submeshDog->getQuadTopology(), init_x0); IsometryObjective isoObj(submeshDog->getQuadTopology(), init_x0);
	DogConstraints dogConst(submeshDog->getQuadTopology());
	QuadraticConstraintsSumObjective dogConstSoft(dogConst, init_x0);
	QuadraticConstraintsSumObjective edgePtSoft(submeshEdgePtConst, init_x0);
	
	Eigen::VectorXd x0(init_x0), x = x0;;

	double infeasability_epsilon = 0.001, infeasability_filter = 0.1; int max_newton_iters = 1; double merit_p = 1;
	double penalty = 1;
	for (int i = 0; i < 10 ; i++) {
		double dogConstWeight = penalty, posConstWeight = penalty;
		CompositeObjective compObj({&bending,&isoObj,&dogConstSoft, &edgePtSoft},
								  {dogEditor.p.bending_weight, dogEditor.p.isometry_weight, dogConstWeight, posConstWeight});
		
		NewtonKKT newtonSolver(infeasability_epsilon, infeasability_filter, max_newton_iters, merit_p);
		EdgePointConstraints emptyConstraints;
		for (int iter = 0; iter < 50; iter++) newtonSolver.solve_constrained(x, compObj,emptyConstraints, x);
		
		x0 = x;
		penalty*=2;
	}
	// TODO clip positional edge point constraints somehow...
	submeshDog->update_V_vector(x);
	globalDog->update_submesh_V(submesh_i, submeshDog->getV());
	delete submeshDog;
}

std::vector< std::vector<int> > DeformationController::submeshes_to_curves(const Dog& dog) {
	auto eS = dog.getEdgeStitching(); auto curves = eS.stitched_curves; int curves_n = curves.size();
	int submesh_n = dog.get_submesh_n();
	std::vector< std::vector<int> > submeshes_to_curves(submesh_n);

	for (int ci = 0; ci < curves_n; ci++) {
		for (auto ep: curves[ci]) {
			// find all edge points there
			auto edge = ep.edge;
			int mult_edge_index = eS.edge_to_duplicates.at(edge);
			int mult_edge_const_start = eS.multiplied_edges_start[mult_edge_index]; 
			int mult_edges_num = eS.multiplied_edges_num[mult_edge_index];

			// Go through all duplicated edges
			for (int set_const_i = mult_edge_const_start; set_const_i < mult_edge_const_start+mult_edges_num; set_const_i++) {
				auto edge1 = eS.edge_const_1[set_const_i], edge2 = eS.edge_const_2[set_const_i];

				submeshes_to_curves[dog.v_to_submesh_idx(edge1.v1)].push_back(ci);
				submeshes_to_curves[dog.v_to_submesh_idx(edge1.v2)].push_back(ci);
				submeshes_to_curves[dog.v_to_submesh_idx(edge2.v1)].push_back(ci);
				submeshes_to_curves[dog.v_to_submesh_idx(edge2.v2)].push_back(ci);
			}
		}
	}
	// Make the lists unique
	for (int subi = 0; subi < submesh_n; subi++) {
		set<int> s(submeshes_to_curves[subi].begin(), submeshes_to_curves[subi].end() );
		submeshes_to_curves[subi].clear();
		submeshes_to_curves[subi].assign( s.begin(), s.end() );
	}

	return submeshes_to_curves;
}