#include "DogSolver.h"

#include "Objectives/StitchingConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/SimplifiedBendingObjective.h"

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


DogSolver::State::State(Dog& dog, const QuadTopology& quadTop, const DogSolver::Params& p) 
					: dog(dog), quadTop(quadTop), p(p),
					flowProject(dog, 1., 1,p.max_lbfgs_routines, p.penalty_repetitions),
          dogGuess(dog, p.align_procrustes, p.arap_guess),
					angleConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.folding_angle),
					curveConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.curve_timestep) {
	// empty on purpose
}

void DogSolver::update_positional_constraints() {
	// TODO support curve constraints as well
	if (p.deformationType == DIHEDRAL_FOLDING) {
		state->angleConstraintsBuilder.get_positional_constraints(b,bc);	
	} else if (p.deformationType == CURVE_DEFORMATION) {
		// update curve constrained folds
		SurfaceCurve surfaceCurve;
		state->curveConstraintsBuilder.get_curve_constraints(surfaceCurve, edgeCoords);
		edgePoints = surfaceCurve.edgePoints;
	}
	
}

void DogSolver::single_optimization() {
	if (!state) return; // No optimizer

  cout << "guessing!" << endl;
  if (state->dog.has_creases()) {
    PositionalConstraints posConst(b,bc);
    StitchingConstraints stitchingConstraints(state->quadTop,state->dog.getEdgeStitching());
    state->dogGuess.guess(state->dog, posConst, stitchingConstraints);
  }

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(state->dog.getV_vector()),x;

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


  //CompositeObjective compObj({&bending, &isoObj,&constObjBesidesPos}, {bending_weight,isometry_weight,const_obj_penalty});
  CompositeObjective compObj({&bending, &isoObj}, {p.bending_weight,p.isometry_weight});
  if (b.rows()) {
    PositionalConstraints posConst(b,bc);
    
    /*
    QuadraticConstraintsSumObjective softPosConst(posConst);
    compObj.add_objective(&softPosConst,const_obj_penalty,true);
    */
    compConst.add_constraints(&posConst);
  }
  

  //state->flowProject.solve_single_iter(x0, compObj, compConst, x);
  //state->flowProject.solve_constrained(x0, compObj, compConst, x);
  LBFGSWithPenalty lbfsgSolver(p.max_lbfgs_routines, p.penalty_repetitions);
  lbfsgSolver.solve_constrained(x0, compObj, compConst, x);

  //state->flowProject.resetSmoother();
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