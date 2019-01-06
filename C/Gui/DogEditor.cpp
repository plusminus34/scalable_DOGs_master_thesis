#include "DogEditor.h"

std::vector<int> get_second_dog_row(Dog& dog);

DogEditor::~DogEditor() {
	if (dogSolver) delete dogSolver;
	if (editor) delete editor;
}

void DogEditor::single_optimization() {
	if (dogSolver) dogSolver->single_iteration(constraints_deviation, objective);
}

void DogEditor::init_from_new_dog(Dog& dog) {
	b.resize(0); bc.resize(0);
	paired_vertices.clear();

	init_x0 = dog.getV_vector();

	if (geoConstraintsBuilder) delete geoConstraintsBuilder;
	geoConstraintsBuilder = new CurveInterpolationConstraintsBuilder(dog.getV(), 
															get_second_dog_row(dog), curve_timestep);
	bool update_solver = false; //update_positional_constraints(update_solver);
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, paired_vertices);
	if (editor) delete editor;
	FTriangular = dog.getFTriangular();
	editor = new Editor(*viewer,dog.getV(), FTriangular, b, bc, paired_vertices, mouse_mode, select_mode);
}

void DogEditor::reset_dog_solver() {
	Dog& dog = dogSolver->getDog();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, paired_vertices);
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

void DogEditor::add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc) {
	auto old_b = b; auto old_bc = bc;
	b.resize(old_b.rows()+new_b.rows()); bc.resize(b.rows());
	b << old_b,new_b; bc << old_bc,new_bc;
	reset_dog_solver();
}

void DogEditor::add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords) {
	edgePoints.insert(edgePoints.end(), new_edgePoints.begin(), new_edgePoints.end());
	auto old_edgeCoords = edgeCoords; edgeCoords.resize(old_edgeCoords.rows()+new_edgeCoords.rows(), old_edgeCoords.cols());
	edgeCoords << old_edgeCoords, new_edgeCoords;

	reset_dog_solver();
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