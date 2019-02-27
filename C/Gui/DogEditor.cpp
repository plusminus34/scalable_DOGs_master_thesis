#include "DogEditor.h"

std::vector<int> get_second_dog_row(Dog& dog);

DogEditor::DogEditor(bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : has_new_constraints(has_new_constraints),b(b),bc(bc),
		paired_vertices(paired_vertices),edgePoints(edgePoints), edgeCoords(edgeCoords) {
		// empty on purpose
}

DogEditor::~DogEditor() {
	if (editor) delete editor;
}

void DogEditor::init_from_new_dog(Dog& dog) {
	b.resize(0); bc.resize(0);
	paired_vertices.clear();

	if (editor) delete editor;
	FTriangular = dog.getFTriangular();
	editor = new Editor(*viewer,dog.getV(), FTriangular, has_new_constraints, b, bc, paired_vertices, mouse_mode, select_mode);
}

/*
std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}*/

bool DogEditor::callback_mouse_down() {
	return editor->callback_mouse_down();
}
bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	return editor->callback_mouse_move(mouse_x, mouse_y);
}