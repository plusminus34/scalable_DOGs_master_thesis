#pragma once

#include "Gui/DogEditor.h"

// The DOG structure is be constant (we have a global DOG that consists of multiple submeshses connected by creases/folds)
// The deformation controller holds pointer to the global Dog and to an "edited submesh DOG"
class DeformationController {
public:
	DeformationController() {globalDog = NULL; editedSubmesh = NULL;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i; dogEditor.init_viewer(viewer_i);}

	// pass it on to the editor
	void single_optimization() {return dogEditor.single_optimization();}
	void propagate_submesh_constraints();
	// if needed: change activeDog and update the editor accordingly
	void update_edited_mesh(int newEditedSubmeshI);
	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}

	DogEditor dogEditor;

private:
	void process_submesh(int submesh_i, const DogEdgeStitching& eS, const std::vector<bool>& passed_on_submesh);

	igl::opengl::glfw::Viewer* viewer;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh	
};