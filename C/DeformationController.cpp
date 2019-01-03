#include "DeformationController.h"

#include <queue>

void DeformationController::init_from_new_dog(Dog& dog) {
	if (globalDog) delete globalDog;
	if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = -1; // Editing the global dog
	dogEditor.init_from_new_dog(dog);
}

void DeformationController::update_edited_mesh(int newEditedSubmeshI) {
	if (globalDog->get_submesh_n() == 1) return; // Only 1 submesh (no folds)
	if (newEditedSubmeshI == editedSubmeshI) return; // no change
	// -1 means the global mesh
	if (newEditedSubmeshI == -1) {
		editedSubmeshI = -1;
		editedSubmesh = globalDog;
		return;
	}
	if ( (newEditedSubmeshI >= globalDog->get_submesh_n()) || (newEditedSubmeshI < -1 ) ) return;

	if (newEditedSubmeshI != editedSubmeshI) {
		if (editedSubmesh && (editedSubmesh!= globalDog)) delete editedSubmesh;
		editedSubmeshI = newEditedSubmeshI;
		editedSubmesh = globalDog->get_submesh(editedSubmeshI);
		dogEditor.init_from_new_dog(*editedSubmesh);
	}
}

void DeformationController::propagate_submesh_constraints() {
	int submesh_n = globalDog->get_submesh_n();
	if ( (editedSubmeshI < 0) || (editedSubmeshI >= submesh_n) ) return;

	auto adjacency_list = globalDog->get_submesh_adjacency();
	auto eS = globalDog->getEdgeStitching();
	std::vector<bool> edge_constraint_set(eS.edge_const_1.size(),false); std::vector<Eigen::RowVector3d> const_value(eS.edge_const_1.size());

	std::vector<bool> passed_on_submesh(submesh_n, false);
	passed_on_submesh[editedSubmeshI] = true; globalDog->update_submesh_V(editedSubmeshI, editedSubmesh->getV());
	update_edge_constraints_from_submesh(editedSubmeshI, eS, edge_constraint_set, const_value);
	std::queue<int> Q; for (int i = 0; i < adjacency_list[editedSubmeshI].size(); i++) Q.push(adjacency_list[editedSubmeshI][i]);
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
	SimplifiedBendingObjective bending(submeshDog->getQuadTopology()); IsometryObjective isoObj(submeshDog->getQuadTopology(), init_x0);
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