#include "DogSolver.h"

#include "Objectives/EqualDiagObjective.h"

#include "../Optimization/CompositeObjective.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/EdgePointConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"

#include "Objectives/FoldingAnglePositionalConstraintsBuilder.h"

using namespace std;


DogSolver::DogSolver(Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0, 
        const DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : 

          dog(dog), quadTop(quadTop), init_x0(init_x0), p(p),
          constraints(dog, quadTop, b, bc, edgePoints, edgeCoords), 
          obj(dog, quadTop, init_x0, constraints.posConst, constraints.edgePtConst),
          newtonKKT(p.merit_p),
          dogGuess(dog, p.align_procrustes) {
    // Empty on purpose
}

DogSolver::Constraints::Constraints(const Dog& dog, const QuadTopology& quadTop,
      Eigen::VectorXi& b, Eigen::VectorXd& bc, 
      std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords) : 
                    dogConst(quadTop),
                    stitchingConstraints(quadTop, dog.getEdgeStitching()),
                    posConst(b,bc),
                    edgePtConst(edgePoints, edgeCoords),
                    compConst({&dogConst, &stitchingConstraints}) /*the positional/edge constraints are soft and go into the objective*/ {
    // Empty on purpose
}

DogSolver::Objectives::Objectives(const Dog& dog, const QuadTopology& quadTop, const Eigen::VectorXd& init_x0,
          PositionalConstraints& posConst,
          EdgePointConstraints& edgePtConst) : 
        bending(quadTop), isoObj(quadTop, init_x0), /*laplacianSimilarity(dog,init_x0),*/
        pointsPosSoftConstraints(posConst),
        edgePosSoftConstraints(edgePtConst),
        compObj(
          {&obj.bending, &obj.isoObj, &obj.pointsPosSoftConstraints, &obj.edgePosSoftConstraints},
          {p.bending_weight,p.isometry_weight. p.soft_pos_weight, p.soft_pos_weight})
        //compObj({&state->obj.bending, &state->obj.isoObj, &state->obj.laplacianSimilarity}, {p.bending_weight,p.isometry_weight,p.laplacian_similarity_weight})
          {
    // Empty on purpose
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
	if (!state) return; // No optimizer

  cout << "guessing!" << endl;
  dogGuess.guess(dog, constraints.posConst, constraints.stitchingConstraints, constraints.edgePtConst);

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(dog.getV_vector()), x(x0);
  obj.compObj.update_weights({p.bending_weight,p.isometry_weight. p.soft_pos_weight, p.soft_pos_weight});

  switch (p.solverType) {
    case SOLVE_NEWTON_PENALTY: {
      /*
      CompositeObjective compObj2;
      compObj2.add_objective(&state->obj.bending,p.bending_weight,true);
      compObj2.add_objective(&state->obj.isoObj,p.isometry_weight,true);
      QuadraticConstraintsSumObjective edgePosConst(edgePtConst);
      compObj2.add_objective(&edgePosConst,p.const_obj_penalty,true);

      double penalty = 1;
      for (int i = 0; i < p.penalty_repetitions  ; i++) {
        CompositeObjective compObj3(compObj2);
        QuadraticConstraintsSumObjective dogConstSoft(state->constraints.dogConst);
        compObj3.add_objective(&dogConstSoft,penalty,true);
        state->newton.solve(x0, compObj3, x); 
        penalty*=2;
      }
      */
      break;
    }
    case SOLVE_NEWTON_FLOW: {
      newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x);
      break;
    }
    case SOLVE_NONE: {
      // Empty on purpose
      break;
    }
  }
  dog.update_V_vector(x);
  
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  objective = obj.compObj.obj(x);
}