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

	std::vector<bool> passed_on_submesh(submesh_n, false);
	passed_on_submesh[editedSubmeshI] = true; globalDog->update_submesh_V(editedSubmeshI, editedSubmesh->getV());
	std::queue<int> Q; for (int i = 0; i < adjacency_list[editedSubmeshI].size(); i++) Q.push(adjacency_list[editedSubmeshI][i]);
	while (!Q.empty()) {
		int cur_submesh = Q.front(); Q.pop();

		// process cur_submesh
		int submesh_v_min_i, submesh_v_max_i;
		globalDog->get_submesh_min_max_i(cur_submesh, submesh_v_min_i, submesh_v_max_i, true);
		auto submeshDog = globalDog->get_submesh(cur_submesh);
		std::vector<EdgePoint> edgePoints; std::vector<Eigen::RowVector3d> edgePointsCoordsList;

		// take all edge constraints that involve this submesh, and are already "set" by other submeshes
		for (int edge_const_i = 0; edge_const_i < eS.edge_const_1.size(); edge_const_i++) {
				double t = eS.edge_coordinates[edge_const_i]; Edge edge_src(eS.edge_const_1[edge_const_i]),edge_target(eS.edge_const_2[edge_const_i]);

				// This makes sure that if this constraint involves the current submesh then the source edge concerns the current submesh
				if (globalDog->v_to_submesh_idx(edge_target.v1) == cur_submesh) {std::swap(edge_src,edge_target);}

				// Now we want to make sure this constraint does involve the submesh, and was already set by other submeshes
				int src_submesh = globalDog->v_to_submesh_idx(edge_src.v1);
				int target_submesh = globalDog->v_to_submesh_idx(edge_target.v1);
				if ((src_submesh == cur_submesh) && (passed_on_submesh[target_submesh]) ) {
					// Add this as edge point constraints
					// needs to change from global V index to local v index
					Edge submeshEdge(edge_src); submeshEdge.v1 -= submesh_v_min_i; submeshEdge.v2 -= submesh_v_min_i;
					edgePoints.push_back(EdgePoint(submeshEdge,t));
					EdgePoint targetEdgePt(edge_target,t);
					edgePointsCoordsList.push_back(targetEdgePt.getPositionInMesh(globalDog->getV()));
				}
		}
		// Convert it to a matrix
		Eigen::MatrixXd edgePointCoords(edgePointsCoordsList.size(),3); 
		for (int edgePtRow = 0; edgePtRow < edgePointCoords.rows(); edgePtRow++) edgePointCoords.row(edgePtRow) = edgePointsCoordsList[edgePtRow];
		EdgePointConstraints submeshEdgePtConst(edgePoints, edgePointCoords);

		globalDog->update_submesh_V(cur_submesh, submeshDog->getV());
		delete submeshDog;

		passed_on_submesh[cur_submesh] = true;
		// Add all submeshes that where not processed
		for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) {
			if (!passed_on_submesh[cur_submesh]) Q.push(adjacency_list[cur_submesh][i]);
		}
	}
}