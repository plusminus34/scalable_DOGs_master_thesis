#include "DogEditor.h"

std::vector<int> get_second_dog_row(Dog& dog);

void DogEditor::single_optimization() {
	if (dogSolver) dogSolver->single_iteration(constraints_deviation, objective);
}

void DogEditor::init_from_new_dog(igl::opengl::glfw::Viewer& viewer_i, Dog& dog, const QuadTopology& quadTop) {
	viewer = &viewer_i;
	auto init_x0 = dog.getV_vector();

	if (geoConstraintsBuilder) delete geoConstraintsBuilder;
	geoConstraintsBuilder = new CurveInterpolationConstraintsBuilder(dog.getV(), 
															get_second_dog_row(dog), curve_timestep);
	bool update_solver = false; //update_positional_constraints(update_solver);
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,quadTop,init_x0, p, b, bc, edgePoints, edgeCoords);
	if (editor) delete editor;
	editor = new Editor(*viewer,dog.getV(), dog.getFrendering(), b, bc, mouse_mode, select_mode);
}

void DogEditor::reset_dog_solver() {
	Dog& dog = dogSolver->getDog();
	auto init_x0 = dog.getV_vector();
	const QuadTopology& quadTop = dogSolver->getQuadTop();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,quadTop,init_x0, p, b, bc, edgePoints, edgeCoords);
}

void DogEditor::update_positional_constraints(bool update_solver) {
	//if (p.deformationType == DIHEDRAL_FOLDING) {
	//	state->angleConstraintsBuilder.get_positional_constraints(b,bc);	
	//} else if 
	//(p.deformationType == CURVE_DEFORMATION) {
	// update curve constrained folds
	SurfaceCurve surfaceCurve;
	geoConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	/*
	if (state->dog.has_creases()) {
		curveConstraintsBuilder.get_curve_constraints(surfaceCurve, edgeCoords);
	} else {
		
	}
	*/
	edgePoints = surfaceCurve.edgePoints;
	//}	

	if (update_solver) dogSolver->update_edge_coords(edgeCoords);
}

std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}

bool DogEditor::callback_mouse_down() {
	auto ret = editor->callback_mouse_down();
	if (editor->new_constraints) {
		reset_dog_solver();
		editor->new_constraints = false;
	}
	return ret;
}
bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	auto ret = editor->callback_mouse_move(mouse_x, mouse_y);
	if (bc.rows()) dogSolver->update_point_coords(bc);
	return ret;
}