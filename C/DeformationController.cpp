#include "DeformationController.h"

void DeformationController::init_from_new_dog(Dog& dog) {
	if (globalDog) delete globalDog;
	if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = 0; // Editing the global dog
	dogEditor.init_from_new_dog(dog);
}

void DeformationController::update_edited_mesh(int newEditedSumeshI) {
	if (newEditedSumeshI != editedSubmeshI) {
		if (editedSubmesh) delete editedSubmesh;
		editedSubmeshI = newEditedSumeshI;
		//editedSubmesh = globalDog->getSubmesh(editedSubmeshI);
		dogEditor.init_from_new_dog(*editedSubmesh);
	}
}