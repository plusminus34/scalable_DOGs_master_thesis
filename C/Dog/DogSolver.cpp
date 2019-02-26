#include "DogSolver.h"

#include "../Folding/CurvedFoldingBiasObjective.h"

using namespace std;


DogSolver::DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, 
        const DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
        std::vector<std::pair<int,int>>& pairs
        /*CurvedFoldingBiasObjective& curvedFoldingBiasObj*/) :

          dog(dog),
          foldingBinormalBiasConstraints(dog),
          init_x0(init_x0), p(p),
          constraints(dog, b, bc, edgePoints, edgeCoords, pairs), 
          obj(dog, init_x0, constraints.posConst, constraints.edgePtConst,constraints.ptPairConst,
                foldingBinormalBiasConstraints, /*curvedFoldingBiasObj,*/ p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          interiorPt(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
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
          /*CurvedFoldingBiasObjective& curvedFoldingBiasObj,*/
          FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
          const DogSolver::Params& p) : 
        bending(dog.getQuadTopology(), init_x0), isoObj(dog.getQuadTopology(), init_x0), laplacianSimilarity(dog,init_x0),
        pointsPosSoftConstraints(posConst, init_x0),
        edgePosSoftConstraints(edgePtConst, init_x0),
        ptPairSoftConst(ptPairConst, init_x0),
        foldingBinormalBiasObj(foldingBinormalBiasConstraints, init_x0),
        /*curvedFoldingBiasObj(curvedFoldingBiasObj),*/
        
        compObj(
          {&bending, &isoObj, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst, &foldingBinormalBiasObj},
          {p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.fold_bias_weight})
          {
    // Empty on purpose
}

bool DogSolver::is_folded() {
  bool is_folded = true;
  auto eS = dog.getEdgeStitching();
  for (int fold_curve_idx = 0; fold_curve_idx < eS.stitched_curves.size(); fold_curve_idx++) {
    const vector<EdgePoint>& foldingCurve = eS.stitched_curves[fold_curve_idx];

    for (int e_idx = 1; e_idx < foldingCurve.size()-1; e_idx++) {
      auto ep_b = foldingCurve[e_idx-1]; auto ep_f = foldingCurve[e_idx+1];
      auto edge_pt = foldingCurve[e_idx]; double edge_t = edge_pt.t;
      if (eS.get_vertex_edge_point_deg(edge_pt.edge) != 1) continue;
      int v1,v2,w1,w2; dog.get_2_submeshes_vertices_from_edge(edge_pt.edge, v1,v2,w1,w2);
      Eigen::VectorXd V1_p(dog.getV().row(v1)), V2_p(dog.getV().row(v2));
      Eigen::VectorXd W1_p(dog.getV().row(w1)), W2_p(dog.getV().row(w2));
      auto ep_p = edge_pt.getPositionInMesh(dog.getV()); auto ep_b_p = ep_b.getPositionInMesh(dog.getV()); auto ep_f_p = ep_f.getPositionInMesh(dog.getV());
      Eigen::VectorXd B = ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).normalized();

      Eigen::VectorXd e1 = V1_p-V2_p, e2 = W1_p-W2_p;
      //std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
      double sign1 = B.dot(e1), sign2 = B.dot(e2);
      if (sign1*sign2 > 0) {
        is_folded = false;
        //cout << "Change!" << endl;
        //std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
        cout << "Curve = " << fold_curve_idx << ", e_idx = " << e_idx << ": sign1 = " << sign1 << " sign2 = " << sign2 << endl;
        //cout << "The entire curve's length is " << foldingCurve.size() << endl;
        //break;
      }
    }
  }
  std::cout << "is_folded = " << is_folded << std::endl;
  return is_folded;
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
  if (is_constrained) {
    cout << "guessing!" << endl;
    dogGuess.guess(dog, constraints.posConst, constraints.stitchingConstraints, constraints.edgePtConst);
  }

	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(dog.getV_vector()), x(x0);
  //obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.fold_bias_weight/*,p.fold_bias_weight*/});
  obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.fold_bias_weight});

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
      //std::cout << "before foldingBinormalBiasConstraints.Vals(x).norm() = " << foldingBinormalBiasConstraints.Vals(x).norm() << std::endl;
      std::pair<Edge,Edge> edge_pairs(Edge(210,220),Edge(220,230));
      EdgesAngleConstraints edgeAngleConstraints(x,{edge_pairs},{0.8});
      QuadraticConstraintsSumObjective edgeAngleConstraintsObj(edgeAngleConstraints,x0);

      
      CompositeObjective compObjNew({&obj.compObj,&edgeAngleConstraintsObj},{1.,10.});
      newtonKKT.solve_constrained(x0, compObjNew, constraints.compConst, x);

      //newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x);
      is_folded();
      //std::cout << "after foldingBinormalBiasConstraints.Vals(x).norm() = " << foldingBinormalBiasConstraints.Vals(x).norm() << std::endl;
      /*
      int wait;// cin >> wait;
      for (int i = 0; i < 100; i++) {
        newtonKKT.solve_constrained(x, obj.compObj, constraints.compConst, x);
        //std::cout << "iter = " << i << " foldingBinormalBiasConstraints.Vals(x).norm() = " << foldingBinormalBiasConstraints.Vals(x).norm() << std::endl;
      }
      */
      //cin >> wait;
      //interiorPt.solve_constrained(x0, obj.compObj, constraints.compConst, mvFoldingConstraints, x);
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