#include "DogSolver.h"

#include "Objectives/EqualDiagObjective.h"

#include "../Optimization/CompositeObjective.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/EdgePointConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"

#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"

using namespace std;

std::vector<int> get_second_dog_row(Dog& dog);
DogSolver::State::State(Dog& dog, const QuadTopology& quadTop, const DogSolver::Params& p, const Eigen::VectorXd& init_x0) 
					: dog(dog), quadTop(quadTop), init_x0(init_x0), p(p), obj(dog, quadTop, init_x0), constraints(dog, quadTop),
					flowProject(dog, 1., 1,p.max_lbfgs_routines, p.penalty_repetitions),
          lbfsgSolver(p.max_lbfgs_routines),
          newtonKKT(p.merit_p),
          dogGuess(dog, p.align_procrustes),
					angleConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.folding_angle),
					curveConstraintsBuilder(dog.getV(), dog.getEdgeStitching(), p.curve_timestep),
          geoConstraintsBuilder(dog.getV(), get_second_dog_row(dog), p.curve_timestep) {

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
  //if (state->dog.has_creases()) {
    PositionalConstraints posConst(b,bc);
    StitchingConstraints stitchingConstraints(state->quadTop,state->dog.getEdgeStitching());
    EdgePointConstraints edgePtConst(edgePoints, edgeCoords);

    state->dogGuess.guess(state->dog, posConst, stitchingConstraints, edgePtConst);
  //}

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(state->dog.getV_vector()), x(x0);

  	// Constraints

  //compConst.add_constraints_permanent(&state->constraints.dogConst);
/*
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
    
  }*/

  //CompositeObjective compObj({&bending, &isoObj,&constObjBesidesPos}, {bending_weight,isometry_weight,const_obj_penalty});
  CompositeObjective compObj({&state->obj.bending, &state->obj.isoObj, &state->obj.laplacianSimilarity}, {p.bending_weight,p.isometry_weight,p.laplacian_similarity_weight});
  
  if (edgeCoords.rows()) {
    EdgePointConstraints edgePtConst(edgePoints, edgeCoords);
    /*
    if (p.solverType != SOLVE_FLOW_PROJECT) {compConst.add_constraints(&edgePtConst);};
    if (p.solverType == SOLVE_FLOW_PROJECT) {
      QuadraticConstraintsSumObjective edgePosConst(edgePtConst);
      compObj.add_objective(&edgePosConst,1,true);
    }
    */
    QuadraticConstraintsSumObjective edgePosConst(edgePtConst);
    compObj.add_objective(&edgePosConst,p.const_obj_penalty,true);
  }
  
  std::vector<EdgePoint> edgePoints; Eigen::MatrixXd edgeCoords;

  switch (p.solverType) {
    case SOLVE_NEWTON_PENALTY: {
      CompositeObjective compObj2;
      compObj2.add_objective(&state->obj.bending,p.bending_weight,true);
      compObj2.add_objective(&state->obj.isoObj,p.isometry_weight,true);
      QuadraticConstraintsSumObjective edgePosConst(edgePtConst);
      compObj2.add_objective(&edgePosConst,p.const_obj_penalty,true);

      double penalty = 1;
      for (int i = 0; i < p.penalty_repetitions  ; i++) {
        CompositeObjective compObj3(compObj2);
        QuadraticConstraintsSumObjective dogConstSoft(state->constraints.dogConst);
        compObj3.add_objective(&dogConstSoft,penalty*p.const_obj_penalty,true);
        state->newton.solve(x0, compObj3, x); 
        penalty*=2;
      }
      break;
    }
    case SOLVE_NEWTON_FLOW: {
      //EqualDiagObjective eqDiag(state->quadTop);
      //compObj.add_objective(&eqDiag,p.diag_length_weight);
      CompositeObjective compObj2;
      compObj2.add_objective_permanent(state->obj.bending,p.bending_weight,true);
      compObj2.add_objective_permanent(state->obj.isoObj,p.isometry_weight,true);
      QuadraticConstraintsSumObjective edgePosConst(edgePtConst);
      compObj2.add_objective_permanent(edgePosConst,p.const_obj_penalty,true);

      state->newtonKKT.solve_constrained(x0, compObj2, state->constraints.compConst, x);
      break;
    }
    case SOLVE_NONE: {
      // Empty on purpose
      break;
    }
  }
  state->dog.update_V_vector(x);
  
  constraints_deviation = state->constraints.compConst.Vals(x).squaredNorm();
  objective = compObj.obj(x);
}

void DogSolver::init_from_new_dog(Dog& dog, const QuadTopology& quadTop) {
  auto init_x0 = dog.getV_vector();
	init_solver_state(dog,quadTop, init_x0);
}

void DogSolver::init_solver_state(Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0) {
	if (state) delete state;
	state = new DogSolver::State(dog,quadTop,p, init_x0);
}