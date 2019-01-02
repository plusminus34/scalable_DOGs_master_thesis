#include "DeformationController.h"

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