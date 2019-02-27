#pragma once

#include "Gui/DogEditor.h"

// The DOG structure is be constant (we have a global DOG that consists of multiple submeshses connected by creases/folds)
// The deformation controller holds pointer to the global Dog and to an "edited submesh DOG"
class DeformationController {
public:
	DeformationController() {globalDog = NULL; editedSubmesh = NULL; dogSolver = NULL; curveConstraintsBuilder = NULL;}
	~DeformationController() {if (dogSolver) delete dogSolver;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i; dogEditor.init_viewer(viewer_i);}

	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}

	void single_optimization();

	void setup_curve_constraints();
	void update_edge_curve_constraints();

	void reset_constraints();
	bool is_folded();

	bool is_curve_constraint = false;
	DogEditor dogEditor;
	DogSolver::Params p;
	double curve_timestep = 0;
	double constraints_deviation;
	double objective;
private:
	void reset_dog_solver();
	EdgePoint find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int& edge_index);

	void update_edge_coords(Eigen::MatrixXd& edgeCoords_i);
	void update_point_coords(Eigen::VectorXd& bc_i);

	igl::opengl::glfw::Viewer* viewer;
	Eigen::VectorXd init_x0;
	// This needs to resest sometimes. 
	// For instance when a new soft constraint is added (but not when the constraint value change), or when a entirely new DOG is loaded
	//	Since this amounts to a different objective/hessian sparsity pattern
	DogSolver* dogSolver;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh	

	CurveInterpolationConstraintsBuilder* curveConstraintsBuilder;
};