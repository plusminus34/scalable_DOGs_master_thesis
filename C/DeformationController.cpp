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
	if (newEditedSubmeshI == editedSubmeshI) return; // no change
	// -1 means the global mesh
	if (newEditedSubmeshI == -1) {
		editedSubmesh = globalDog;
		return;
	}
	//if (globalDog->get_submesh_n() == 1) {return;} // There are no submeshes
	// check submesh range
	//std::cout << "newEditedSubmeshI = " << newEditedSubmeshI << " and globalDog->get_submesh_n() = " << globalDog->get_submesh_n() << std::endl;
	if ( (newEditedSubmeshI >= globalDog->get_submesh_n()) || (newEditedSubmeshI < -1 ) ) return;

	if (newEditedSubmeshI != editedSubmeshI) {
		if (editedSubmesh && (editedSubmesh!= globalDog)) delete editedSubmesh;
		editedSubmeshI = newEditedSubmeshI;
		editedSubmesh = globalDog->get_submesh(editedSubmeshI);
		dogEditor.init_from_new_dog(*editedSubmesh);
	}
}