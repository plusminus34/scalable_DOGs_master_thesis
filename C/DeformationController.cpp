#include "DeformationController.h"

#include <queue>
using namespace std;

DeformationController::DeformationController() : dogEditor(b,bc,paired_vertices,edgePoints,edgeCoords), globalDog(NULL),
			 editedSubmesh(NULL), dogSolver(NULL), curveConstraintsBuilder(NULL) {
	// empty on purpose
}

void DeformationController::init_from_new_dog(Dog& dog) {
	if (globalDog) delete globalDog;
	if (editedSubmesh) delete editedSubmesh;
	globalDog = &dog;
	editedSubmesh = globalDog;
	editedSubmeshI = -1; // Editing the global dog

	dogEditor.init_from_new_dog(dog);

	init_x0 = dog.getV_vector();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, paired_vertices);
}

void DeformationController::single_optimization() {
	if (dogEditor.has_new_constraints()) {
		reset_dog_solver();
		dogEditor.signal_handled_new_constraints();
	}
	if (is_curve_constraint) update_edge_curve_constraints();
	if (dogEditor.has_new_point_constraints()) dogSolver->update_point_coords(bc);
	dogSolver->single_iteration(constraints_deviation, objective);
}

void DeformationController::setup_curve_constraints() {
	if (curveConstraintsBuilder) delete curveConstraintsBuilder;
	curveConstraintsBuilder = new CurveInterpolationConstraintsBuilder(globalDog->getV(), 
															globalDog->getEdgeStitching(), curve_timestep);
	SurfaceCurve surfaceCurve; Eigen::MatrixXd edgeCoords;
	curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	is_curve_constraint = true;
	dogEditor.add_edge_point_constraints(surfaceCurve.edgePoints,edgeCoords);
}

void DeformationController::update_edge_curve_constraints() {
	SurfaceCurve surfaceCurve; Eigen::MatrixXd edgeCoords;
	curveConstraintsBuilder->get_curve_constraints(surfaceCurve, edgeCoords);
	update_edge_coords(edgeCoords);
}

void DeformationController::update_edge_coords(Eigen::MatrixXd& edgeCoords_i) {
	edgeCoords = edgeCoords_i; dogSolver->update_edge_coords(edgeCoords);
}

void DeformationController::update_point_coords(Eigen::VectorXd& bc_i) {
	bc = bc_i; dogSolver->update_point_coords(bc);
}

void DeformationController::reset_constraints() {
	b.resize(0);
	bc.resize(0); 
	paired_vertices.clear(); 
	edgePoints.clear(); 
	edgeCoords.resize(0,3); 
	dogEditor.editor->clearHandles(); 
	reset_dog_solver(); 
	is_curve_constraint = false;
}

// t = 0.5 in the edge constraint means it is equally spaced
EdgePoint DeformationController::find_most_equally_spaced_edge_on_fold_curve(int fold_curve_idx, int &min_edge) {
	auto eS = globalDog->getEdgeStitching(); const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];
	int curve_v_n = foldingCurve.size();
	min_edge = 1; double min_dist_from_equal = abs(0.5-foldingCurve[1].t);
	for (int ei = 1; ei < foldingCurve.size()-1; ei++) {
		double dist_from_equal = abs(0.5-foldingCurve[ei].t);
		if ( dist_from_equal < min_dist_from_equal) {
			min_edge = ei;
			min_dist_from_equal = dist_from_equal;
		}
	}
	return foldingCurve[min_edge];
}

void DeformationController::reset_dog_solver() {
	Dog& dog = dogSolver->getDog();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc, edgePoints, edgeCoords, paired_vertices);
}
