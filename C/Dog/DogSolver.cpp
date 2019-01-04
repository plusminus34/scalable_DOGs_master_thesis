#include "DogSolver.h"

using namespace std;


DogSolver::DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, 
        const DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
        std::vector<std::pair<int,int>>& pairs) : 

          dog(dog), init_x0(init_x0), p(p),
          constraints(dog, b, bc, edgePoints, edgeCoords, pairs), 
          obj(dog, init_x0, constraints.posConst, constraints.edgePtConst,constraints.ptPairConst,p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          dogGuess(dog, p.align_procrustes) {
    
    is_constrained = (b.rows() + edgePoints.size())>0;
}

DogSolver::Constraints::Constraints(const Dog& dog,
      Eigen::VectorXi& b, Eigen::VectorXd& bc, 
      std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
      std::vector<std::pair<int,int>>& pairs) : 
                    dogConst(dog.getQuadTopology()),
                    stitchingConstraints(dog.getQuadTopology(), dog.getEdgeStitching()),
                    posConst(b,bc),
                    edgePtConst(edgePoints, edgeCoords),
                    ptPairConst(pairs),
                    compConst({&dogConst, &stitchingConstraints}) /*the positional/edge constraints are soft and go into the objective*/ {
    // Empty on purpose
}

DogSolver::Objectives::Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
          PositionalConstraints& posConst,
          EdgePointConstraints& edgePtConst,
          PointPairConstraints& ptPairConst,
          const DogSolver::Params& p) : 
        bending(dog.getQuadTopology(), init_x0), isoObj(dog.getQuadTopology(), init_x0), laplacianSimilarity(dog,init_x0),
        pointsPosSoftConstraints(posConst, init_x0),
        edgePosSoftConstraints(edgePtConst, init_x0),
        ptPairSoftConst(ptPairConst, init_x0),
        
        compObj(
          {&bending, &isoObj, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst},
          {p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight})
          {
    // Empty on purpose
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
  if (is_constrained) {
    cout << "guessing!" << endl;
    dogGuess.guess(dog, constraints.posConst, constraints.stitchingConstraints, constraints.edgePtConst);
  }

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(dog.getV_vector()), x(x0);
  obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight});

  /*
  std::cout << "current objective = " << obj.compObj.obj(x0) << std::endl; int wait; std::cin >> wait;
  std::cout <<"obj.bending.obj(x0) = " << obj.bending.obj(x0) << std::endl; std::cin >> wait;
  std::cout <<"obj.isoObj.obj(x0) = " << obj.isoObj.obj(x0) << std::endl; std::cin >> wait;
  */

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