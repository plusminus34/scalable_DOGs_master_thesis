#include "DogSolver.h"

#include "Objectives/StitchingConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"
#include "Objectives/LaplacianSimilarity.h"

#include "../Optimization/CompositeObjective.h"
#include "../Optimization/CompositeConstraints.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/EdgePointConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"
#include "../Optimization/Solvers/LBFGS.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"
#include "Objectives/StitchingConstraints.h"

using namespace std;

std::vector<int> get_second_dog_row(Dog& dog);
DogSolver::State::State(Dog& dog, const QuadTopology& quadTop, const DogSolver::Params& p) 
					: dog(dog), quadTop(quadTop), p(p),
					flowProject(dog, 1., 1,p.max_lbfgs_routines, p.penalty_repetitions),
          dogGuess(dog, p.align_procrustes, p.arap_guess),
					angleConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.folding_angle),
					curveConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.curve_timestep),
          geoConstraintsBuilder(dog.getV(), get_second_dog_row(dog), p.curve_timestep) {
	// empty on purpose
}

void DogSolver::update_positional_constraints() {
	// TODO support curve constraints as well
	if (p.deformationType == DIHEDRAL_FOLDING) {
		state->angleConstraintsBuilder.get_positional_constraints(b,bc);	
	} else if (p.deformationType == CURVE_DEFORMATION) {
		// update curve constrained folds
		SurfaceCurve surfaceCurve;
    if (state->dog.has_creases()) {
      state->curveConstraintsBuilder.get_curve_constraints(surfaceCurve, edgeCoords);
    } else {
      state->geoConstraintsBuilder.get_curve_constraints(surfaceCurve, edgeCoords);
    }
		edgePoints = surfaceCurve.edgePoints;
	}
	
}

std::vector<int> get_second_dog_row(Dog& dog) {
  std::vector<int> curve_i; int v_n = dog.getV().rows();
  for (int i = sqrt(v_n); i < 2*sqrt(v_n); i++) {curve_i.push_back(i);}
  return curve_i;
}

void DogSolver::single_optimization() {
	if (!state) return; // No optimizer

  cout << "guessing!" << endl;
  if (state->dog.has_creases()) {
    PositionalConstraints posConst(b,bc);
    StitchingConstraints stitchingConstraints(state->quadTop,state->dog.getEdgeStitching());
    EdgePointConstraints edgePtConst(edgePoints, edgeCoords);

    state->dogGuess.guess(state->dog, posConst, stitchingConstraints, edgePtConst);
  }

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(state->dog.getV_vector()), x(x0);

  	// Constraints
  	DogConstraints dogConst(state->quadTop);

  CompositeConstraints compConst({&dogConst});

  if (state->dog.has_creases()) {
    StitchingConstraints stitchingConstraints(state->quadTop,state->dog.getEdgeStitching());
    compConst.add_constraints(&stitchingConstraints);

    // Check for any positional constraints (for now these will only be folding constraints)
    /*
    if (state->b.rows()) {
      PositionalConstraints posConst(state->b,state->bc);
      compConst.add_constraints(&posConst);

      FoldingAngleConstraints
    }
    */
  }

  // Objectives
  SimplifiedBendingObjective bending(state->quadTop);
  IsometryObjective isoObj(state->quadTop,x0);
  QuadraticConstraintsSumObjective constObjBesidesPos(compConst);
  LaplacianSimilarity laplacianSimilarity(state->dog,x0);


  //CompositeObjective compObj({&bending, &isoObj,&constObjBesidesPos}, {bending_weight,isometry_weight,const_obj_penalty});
  CompositeObjective compObj({&bending, &isoObj, &laplacianSimilarity}, {p.bending_weight,p.isometry_weight,p.laplacian_similarity_weight});
  if (b.rows()) {
    PositionalConstraints posConst(b,bc);
    
    /*
    QuadraticConstraintsSumObjective softPosConst(posConst);
    compObj.add_objective(&softPosConst,const_obj_penalty,true);
    */
    compConst.add_constraints(&posConst);
  }
  if (edgeCoords.rows()) {
    EdgePointConstraints edgePtConst(edgePoints, edgeCoords);
    compConst.add_constraints(&edgePtConst);
  }
  
  std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;

  switch (p.solverType) {
    case SOLVE_FLOW_PROJECT: {
      state->flowProject.solve_single_iter(x0, compObj, compConst, x);
      //state->flowProject.solve_constrained(x0, compObj, compConst, x);
      state->flowProject.resetSmoother();
      break;
    }
    case SOLVE_PENALTY: {
      LBFGSWithPenalty lbfsgSolver(p.max_lbfgs_routines, p.penalty_repetitions);
      lbfsgSolver.solve_constrained(x0, compObj, compConst, x);
    }
    case SOLVE_IPOPT: {
      //TODO implement
      break;
    }
    case SOLVE_NONE: {
      // Empty on purpose
      break;
    }
  }
  state->dog.update_V_vector(x);
  
  constraints_deviation = compConst.Vals(x).squaredNorm();
  objective = compObj.obj(x);
}

void DogSolver::init_from_new_dog(Dog& dog, const QuadTopology& quadTop) {
	init_solver_state(dog,quadTop);
}

void DogSolver::init_solver_state(Dog& dog, const QuadTopology& quadTop) {
	if (state) delete state;
	state = new DogSolver::State(dog,quadTop,p);
}