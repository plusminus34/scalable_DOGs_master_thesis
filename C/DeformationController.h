#pragma once

#include "Dog/DogSolver.h"
#include "Dog/Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Dog/Objectives/FoldingDihedralAngleConstraintsBuilder.h"
#include "Gui/DogEditor.h"

class DeformationController {
public:
	DeformationController();
	~DeformationController() {if (dogSolver) delete dogSolver;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}
	void apply_new_editor_constraint();
	void reset_new_editor_constraint() {dogEditor->reset_new_constraint();}

	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};

	void single_optimization();

	void setup_curve_constraints();
	void update_edge_curve_constraints();

	void reset_constraints();
	bool is_folded();

	DogEditor::EditMode edit_mode = DogEditor::NONE;
	DogEditor::SelectMode select_mode = DogEditor::VertexPicker;

	bool has_new_constraints = false;
	bool is_curve_constraint = false;
	DogEditor* dogEditor;
	DogSolver::Params p;
	double deformation_timestep = 0;
	double constraints_deviation;
	double objective;

private:
	void reset_dog_solver();
	EdgePoint find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int& edge_index);

	void add_edge_point_constraints(const std::vector<EdgePoint>& new_edgePoints, const Eigen::MatrixXd& new_edgeCoords);
	void update_edge_coords(Eigen::MatrixXd& edgeCoords_i);
	void update_point_coords(Eigen::VectorXd& bc_i);
	void add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc);
	void add_edge_point_constraint(const EdgePoint& new_edgePoints, const Eigen::RowVector3d& new_edgeCoords);
	void add_pair_vertices_constraints(const std::vector<std::pair<int,int>>& new_pair_vertices);
	void add_pair_vertices_constraint(int v1, int v2);

	igl::opengl::glfw::Viewer* viewer;

	// Points positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	// Edge positional constraints
	std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;
	// Point pair constraints
	std::vector<std::pair<int,int>> paired_vertices;
	// Edge angle constraints
	std::vector<std::pair<Edge,Edge>> edge_angle_pairs; std::vector<double> edge_cos_angles;

	Eigen::VectorXd init_x0;
	// This needs to resest sometimes. 
	// For instance when a new soft constraint is added (but not when the constraint value change), or when a entirely new DOG is loaded
	//	Since this amounts to a different objective/hessian sparsity pattern
	DogSolver* dogSolver;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh	

	CurveInterpolationConstraintsBuilder* curveConstraintsBuilder;
	FoldingDihedralAngleConstraintsBuilder* foldingDihedralAngleConstraintsBuilder;
};