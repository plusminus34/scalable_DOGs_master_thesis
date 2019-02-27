#include "DogEditor.h"

std::vector<int> get_second_dog_row(Dog& dog);

DogEditor::DogEditor(Eigen::VectorXi& b, Eigen::VectorXd& bc, std::vector<std::pair<int,int>>& paired_vertices,
		std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : b(b),bc(bc),paired_vertices(paired_vertices),
		edgePoints(edgePoints), edgeCoords(edgeCoords) {
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
	editor = new Editor(*viewer,dog.getV(), FTriangular, b, bc, paired_vertices, mouse_mode, select_mode);
}

void DogEditor::add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords) {
	edgePoints.insert(edgePoints.end(), new_edgePoints.begin(), new_edgePoints.end());
	Eigen::MatrixXd old_edgeCoords = edgeCoords; edgeCoords.resize(old_edgeCoords.rows()+new_edgeCoords.rows(), old_edgeCoords.cols());
	if (old_edgeCoords.rows()) edgeCoords << old_edgeCoords, new_edgeCoords; else edgeCoords = new_edgeCoords; // Eigen's concatenate crashes if one of them is empty

	editor->new_constraints = true;
}
/*
void DogEditor::add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc) {
	Eigen::VectorXi old_b = b; Eigen::VectorXd old_bc = bc;
	b.resize(old_b.rows()+new_b.rows()); bc.resize(b.rows());
	if (old_b.rows()) {
		b << old_b,new_b; bc << old_bc,new_bc;
	} else {
		b = new_b; bc = new_bc; // Eigen's concatenate crashes if one of them is empty
	}
	
	reset_dog_solver();
}

void DogEditor::add_edge_point_constraint(const EdgePoint& new_edgePoint, const Eigen::RowVector3d& new_edgeCoords) {
	std::vector<EdgePoint> new_edgePoints = {new_edgePoint}; Eigen::MatrixXd newCoords(1,3); newCoords.row(0) = new_edgeCoords;
	add_edge_point_constraints(new_edgePoints, newCoords);
}

void DogEditor::add_pair_vertices_constraints(const std::vector<std::pair<int,int>>& new_pair_vertices) {
	paired_vertices.insert(paired_vertices.end(), new_pair_vertices.begin(), new_pair_vertices.end());
	reset_dog_solver();
}

void DogEditor::add_pair_vertices_constraint(int v1, int v2) {
	std::vector<std::pair<int,int>> new_pair_vertices{std::pair<int,int>(v1,v2)}; 
	add_pair_vertices_constraints(new_pair_vertices);
}

std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}*/

bool DogEditor::callback_mouse_down() {
	return editor->callback_mouse_down();
	/*
	auto ret = editor->callback_mouse_down();
	if (editor->new_constraints) {
		reset_dog_solver();
		editor->new_constraints = false;
	}
	return ret;
	*/
}
bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	new_point_constraints = (bc.rows() > 0);
	return editor->callback_mouse_move(mouse_x, mouse_y);
}