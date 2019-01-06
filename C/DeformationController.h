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
	void setup_fold_constraints();
	void propagate_submesh_constraints();
	// if needed: change activeDog and update the editor accordingly
	void update_edited_mesh(int newEditedSubmeshI);
	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}

	DogEditor dogEditor;

private:
	void update_edge_constraints_from_submesh(int submesh_i, const DogEdgeStitching& eS, 
									std::vector<bool>& edge_constraint_set, std::vector<Eigen::RowVector3d>& const_value);
	void deform_submesh_based_on_previous_submeshes(int submesh_i, const DogEdgeStitching& eS, 
		 const std::vector<bool>& edge_constraint_set, const std::vector<Eigen::RowVector3d>& const_value);
	Edge find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx);

	igl::opengl::glfw::Viewer* viewer;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh	
};