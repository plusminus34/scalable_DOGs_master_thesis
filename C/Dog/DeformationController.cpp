#include "DeformationController.h"

std::vector<int> get_second_dog_row(Dog& dog);

void DeformationController::init_from_new_dog(Dog& dog, const QuadTopology& quadTop) {
	auto init_x0 = dog.getV_vector();

	if (geoConstraintsBuilder) delete geoConstraintsBuilder;
	geoConstraintsBuilder = new CurveInterpolationConstraintsBuilder(dog.getV(), 
															get_second_dog_row(dog), p.curve_timestep);
	update_positional_constraints();

	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,quadTop,init_x0, p, b, bc, edgePoints, edgeCoords);
}

void DeformationController::update_positional_constraints() {
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

	dogSolver.update_edge_coords(edgeCoords);
}

std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}