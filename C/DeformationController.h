#pragma once

#include "Gui/DogEditor.h"

#include "Folding/MVFolds.h"
#include "Folding/ReflectionFolds.h"

// The DOG structure is be constant (we have a global DOG that consists of multiple submeshses connected by creases/folds)
// The deformation controller holds pointer to the global Dog and to an "edited submesh DOG"
class DeformationController {
public:
	DeformationController() : dogEditor() {globalDog = NULL; editedSubmesh = NULL; curveConstraintsBuilder = NULL;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i; dogEditor.init_viewer(viewer_i);}

	// pass it on to the editor
	void single_optimization();

	
	void setup_curve_constraints();
	void update_edge_curve_constraints();

	// if needed: change activeDog and update the editor accordingly
	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}

	void reset_constraints() {dogEditor.reset_constraints(); is_curve_constraint = false;}
	bool is_folded();

	bool is_curve_constraint = false;
	DogEditor dogEditor;
	double curve_timestep = 0;
private:
	
	EdgePoint find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int& edge_index);

	igl::opengl::glfw::Viewer* viewer;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh	

	CurveInterpolationConstraintsBuilder* curveConstraintsBuilder;
};