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
	std::cout << "newEditedSubmeshI = " << newEditedSubmeshI << " globalDog->get_submesh_n() = " << globalDog->get_submesh_n() << std::endl;
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
		std::cout << "editedSubmesh->getV().rows() = " << editedSubmesh->getV().rows() << std::endl;
		std::cout << "editedSubmesh->getF().rows() = " << editedSubmesh->getF().rows() << std::endl;
		dogEditor.init_from_new_dog(*editedSubmesh);
		std::cout << "alive?" << std::endl;
	}
}