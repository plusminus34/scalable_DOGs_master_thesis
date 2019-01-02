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
	auto adjacency_list = globalDog->get_submesh_adjacency();
	auto eS = globalDog->getEdgeStitching();

	std::vector<Edge> edge_const_1, edge_const_2;

	if ( (editedSubmeshI < 0) || (editedSubmeshI >= submesh_n) ) return;
	std::vector<bool> passed_on_submesh(submesh_n, false);
	passed_on_submesh[editedSubmeshI] = true; globalDog->update_submesh_V(editedSubmeshI, editedSubmesh->getV());
	std::queue<int> Q; for (int i = 0; i < adjacency_list[editedSubmeshI].size(); i++) Q.push(adjacency_list[editedSubmeshI][i]);
	while (!Q.empty()) {
		int cur_submesh = Q.front(); Q.pop();

		// process cur_submesh
		auto submeshDog = globalDog->get_submesh(cur_submesh);
		// take all edge constraints that involve this submesh, and are already "set" by other submeshes

		globalDog->update_submesh_V(cur_submesh, submeshDog->getV());
		delete submeshDog;

		passed_on_submesh[cur_submesh] = true;
		// Add all submeshes that where not processed
		for (int i = 0; i < adjacency_list[cur_submesh].size(); i++) {
			if (!passed_on_submesh[cur_submesh]) Q.push(adjacency_list[cur_submesh][i]);
		}
	}
}