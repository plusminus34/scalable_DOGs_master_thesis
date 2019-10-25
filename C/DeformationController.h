#pragma once

#include "igl/serialize.h"

#include "Dog/DogSolver.h"
#include "Dog/Objectives/CurveInterpolationConstraintsBuilder.h"
#include "Dog/Objectives/FoldingDihedralAngleConstraintsBuilder.h"
#include "Dog/Objectives/MVFoldingDihedralAngleConstraintsBuilder.h"
#include "Gui/DogEditor.h"


enum WallapaperType {
	XY = 0,
	UXY = 1,
	UXUY = 2,
	XUY = 3
};

class DeformationController : public igl::Serializable {
public:
	DeformationController();
	~DeformationController() {if (dogSolver) delete dogSolver;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}
	void apply_new_editor_constraint();
	void reset_new_editor_constraint() {dogEditor->reset_new_constraint();}
	void setup_optimization_measurements(std::string log_file_name);

	bool has_constraints();

	const Dog* getEditedSubmesh() const {return editedSubmesh;}
	int getEditedSubmeshI() const {return editedSubmeshI;}
	void change_submesh(int submesh_i = -2);

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void get_edge_point_constraints(std::vector<EdgePoint>& edgePoints_out, Eigen::MatrixXd& edgeCoords_out) const {edgePoints_out = edgePoints; edgeCoords_out = edgeCoords;};

	void single_optimization();

	void setup_curve_constraints();
	void update_time_deformations() {update_edge_curve_constraints();update_dihedral_constraints();}
	void update_edge_curve_constraints();
	void update_dihedral_constraints();
	void set_cylindrical_boundary_constraints();
	void set_wallpaper_constraints();

	void reset_constraints();
	bool is_folded();

	void InitSerialization() {
	      Add(edit_mode,std::string("edit_mode"));
	      Add(select_mode,std::string("select_mode"));
	      Add(deformed_curve_idx,std::string("deformed_curve_idx"));

	      Add(curve_k_translation,std::string("curve_k_translation")); Add(curve_k_mult,std::string("curve_k_mult")); Add(curve_t_addition,std::string("curve_t_addition"));
	      Add(max_curve_points,std::string("max_curve_points"));
	      Add(z_only_editing,std::string("z_only_editing"));

	      Add(has_new_constraints,std::string("has_new_constraints"));
	      Add(is_curve_constraint,std::string("is_curve_constraint"));
	      Add(is_time_dependent_deformation,std::string("is_time_dependent_deformation"));

	      Add(p,std::string("deformation_params"));
	      Add(deformation_timestep,std::string("deformation_timestep"));
	      Add(deformation_timestep_diff,std::string("deformation_timestep_diff"));
	      Add(paired_boundary_bending_weight_mult,std::string("paired_boundary_bending_weight_mult"));
	      Add(dst_dihedral_angle,std::string("dst_dihedral_angle"));
	      Add(constraints_deviation,std::string("constraints_deviation"));
	      Add(objective,std::string("objective"));

	      Add(matching_curve_pts_x,std::string("matching_curve_pts_x"));
	      Add(matching_curve_pts_y,std::string("matching_curve_pts_y"));

	      Add(dihedral_constrained,std::string("dihedral_constrained"));

	      Add(b,std::string("b"));
	      Add(bc,std::string("bc"));
	      Add(edgePoints,std::string("edgePoints"));
	      Add(edgeCoords,std::string("edgeCoords"));
	      Add(paired_vertices,std::string("paired_vertices"));
	      Add(edge_angle_pairs,std::string("edge_angle_pairs")); Add(edge_cos_angles,std::string("edge_cos_angles"));
	      Add(bnd_vertices_pairs,std::string("bnd_vertices_pairs"));
	      Add(init_x0,std::string("init_x0"));
    }

	DogEditor::EditMode edit_mode = DogEditor::NONE;
	DogEditor::SelectMode select_mode = DogEditor::VertexPicker;
	DogSolver::SolverMode solver_mode = DogSolver::mode_experimental;

	int deformed_curve_idx = 0;
	double curve_k_translation = 0; double curve_k_mult = 2; double curve_t_addition = 0;
	int max_curve_points = 1000; // maximum constrained curve points
	bool z_only_editing = false;

	bool has_new_constraints = false;
	bool is_curve_constraint = false;
	bool is_time_dependent_deformation = false;

	DogEditor* dogEditor;
	DogSolver::Params p;
	double deformation_timestep = 0;
	double deformation_timestep_diff = 0.004;
	double paired_boundary_bending_weight_mult = 1;
	double dst_dihedral_angle = 0; // Later possibly have a list for every edge_angle_pair, so we could have different ones
	double constraints_deviation;
	double objective;

	int stored_iterations = 0 ;
	void store_data(int num_iterations);

	std::pair<std::vector<int>,std::vector<int>> matching_curve_pts_x;
	std::pair<std::vector<int>,std::vector<int>> matching_curve_pts_y;

	std::vector<EdgePoint> dihedral_constrained; // used for also plotting the dihedral constraints
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
	// Edge angle constraint and dihedral constraints
	std::vector<std::pair<Edge,Edge>> edge_angle_pairs; std::vector<double> edge_cos_angles;
	// MV dihedral angle folds data (typically used for straight creases but sometimes for dihedral angle with an M/V assignment)
	std::vector<MVTangentCreaseFold> mvTangentCreaseAngleParams; std::vector<double> mv_cos_angles;

	// Pairs matching boundary curvature
	std::vector<std::pair<int,int>> bnd_vertices_pairs;


	Eigen::VectorXd init_x0;
	// This needs to reset sometimes.
	// For instance when a new soft constraint is added (but not when the constraint value change), or when a entirely new DOG is loaded
	//	Since this amounts to a different objective/hessian sparsity pattern
	DogSolver* dogSolver;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh

	CurveInterpolationConstraintsBuilder* curveConstraintsBuilder;
	FoldingDihedralAngleConstraintsBuilder* foldingDihedralAngleConstraintsBuilder;
	MVFoldingDihedralAngleConstraintsBuilder* mvFoldingDihedralAngleConstraintsBuilder;

	bool optimization_measurements;
	std::ofstream* opt_measurements_log;

	// Objective data stored
	int current_iteration = 0;
	std::vector<double> obj_data;
};
