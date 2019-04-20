#include "DogSolver.h"

using namespace std;


DogSolver::DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, 
        DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
        std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
        std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
        std::vector<std::pair<int,int>>& pairs,
        std::pair<vector<int>,vector<int>>& matching_curve_pts_y,
        std::pair<vector<int>,vector<int>>& matching_curve_pts_x,
        std::ofstream* time_measurements_log) :

          dog(dog),
          foldingBinormalBiasConstraints(dog),
          init_x0(init_x0), p(p),
          constraints(dog, init_x0, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles, mvTangentCreaseAngleParams, 
                      mv_cos_angles, pairs),
          obj(dog, init_x0, constraints, foldingBinormalBiasConstraints, matching_curve_pts_y, matching_curve_pts_x, p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          //interiorPt(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          time_measurements_log(time_measurements_log)
           {
    
    is_constrained = (b.rows() + edgePoints.size())>0;
    if (time_measurements_log) {
      p.max_newton_iters = 1;
    }
}

DogSolver::Constraints::Constraints(const Dog& dog, const Eigen::VectorXd& init_x0,
      Eigen::VectorXi& b, Eigen::VectorXd& bc, 
      std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
      std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
      std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
      std::vector<std::pair<int,int>>& pairs) : 
                    dogConst(dog.getQuadTopology()),
                    stitchingConstraints(dog.getQuadTopology(), dog.getEdgeStitching()),
                    posConst(b,bc),
                    edgePtConst(edgePoints, edgeCoords),
                    edgeAngleConst(init_x0, edge_angle_pairs, edge_cos_angles),
                    mvTangentCreaseAngleConst(init_x0, mvTangentCreaseAngleParams, mv_cos_angles),
                    ptPairConst(pairs),
                    compConst({&dogConst, &stitchingConstraints}) /*the positional/edge constraints are soft and go into the objective*/ {
    // Empty on purpose
}

DogSolver::Objectives::Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
        Constraints& constraints,/*
          PositionalConstraints& posConst,
          EdgePointConstraints& edgePtConst,
          EdgesAngleConstraints& edgeAnglesConst,
          PointPairConstraints& ptPairConst,*/
          FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
          std::pair<vector<int>,vector<int>>& matching_curve_pts_y,
          std::pair<vector<int>,vector<int>>& matching_curve_pts_x,
          const DogSolver::Params& p) : 
        bending(dog.getQuadTopology(), init_x0), isoObj(dog.getQuadTopology(), init_x0),
        pointsRigidAlignmentY(matching_curve_pts_y.first, matching_curve_pts_y.second),
        pointsRigidAlignmentX(matching_curve_pts_x.first, matching_curve_pts_x.second),
        pointsPosSoftConstraints(constraints.posConst, init_x0),
        edgePosSoftConstraints(constraints.edgePtConst, init_x0),
        edgeAnglesSoftConstraints(constraints.edgeAngleConst, init_x0),
        mvTangentCreaseSoftConstraints(constraints.mvTangentCreaseAngleConst, init_x0),
        ptPairSoftConst(constraints.ptPairConst, init_x0),
        foldingBinormalBiasObj(foldingBinormalBiasConstraints, init_x0),
        //allConstQuadraticObj(constraints, init_x0),
        /*curvedFoldingBiasObj(curvedFoldingBiasObj),*/
        
        compObj(
          {&bending, &isoObj, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst, &edgeAnglesSoftConstraints, &mvTangentCreaseSoftConstraints, &foldingBinormalBiasObj, &pointsRigidAlignmentY, &pointsRigidAlignmentX},
          {p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.wallpaper_curve_weight,p.wallpaper_curve_weight})
          /*
        compObj(
          {&bending, &isoObj, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst, &edgeAnglesSoftConstraints, &foldingBinormalBiasObj, &allConstQuadraticObj},
          {p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.fold_bias_weight,1})
          */
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
      if ((eS.get_vertex_edge_point_deg(edge_pt.edge) == 1) && dog.is_crease_vertex_flat(fold_curve_idx,e_idx) ) continue;
      int v1,v2,w1,w2; dog.get_2_submeshes_vertices_from_edge(edge_pt.edge, v1,v2,w1,w2);
      Eigen::VectorXd V1_p(dog.getV().row(v1)), V2_p(dog.getV().row(v2));
      Eigen::VectorXd W1_p(dog.getV().row(w1)), W2_p(dog.getV().row(w2));
      auto ep_p = edge_pt.getPositionInMesh(dog.getV()); auto ep_b_p = ep_b.getPositionInMesh(dog.getV()); auto ep_f_p = ep_f.getPositionInMesh(dog.getV());
      Eigen::VectorXd B = ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).normalized();

      Eigen::VectorXd e1 = V1_p-V2_p, e2 = W1_p-W2_p;
      //std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
      double sign1 = B.dot(e1), sign2 = B.dot(e2); double flat_tolerance = 1e-12; // ignore flat points..
      //if ( (sign1*sign2 > 0) && (((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() > flat_tolerance) ) {
      if ( sign1*sign2 > 0) {
        is_folded = false;
        //cout << "Change!" << endl;
        std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
        cout << "Curve = " << fold_curve_idx << ", e_idx = " << e_idx << ": sign1 = " << sign1 << " sign2 = " << sign2 << " sign1*sign2 = " << sign1*sign2 << endl;
        //cout << "The entire curve's length is " << foldingCurve.size() << endl;
        break;
      }
    }
  }
  //std::cout << "is_folded = " << is_folded << std::endl;
  return is_folded;
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
  if (p.folding_mode) return single_iteration_fold(constraints_deviation, objective);
  else return single_iteration_normal(constraints_deviation, objective);
}

void DogSolver::single_iteration_fold(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(dog.getV_vector()), x(x0);
  if (!is_folded()) {
    cout << "Error: Not folded" << endl; exit(1);
  }
  //obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.fold_bias_weight/*,p.fold_bias_weight*/});
  obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.wallpaper_curve_weight,p.wallpaper_curve_weight});
  //obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.fold_bias_weight,1});
  newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
  dog.update_V_vector(x);
  
  if (is_folded()) {
    p.fold_bias_weight = 1;
  } else {
    while (!is_folded() && (p.fold_bias_weight < 1e14)) {
      cout << "Not folded, fold bias = " << p.fold_bias_weight << " reverting back and making the bias stronger" << endl;
      x = x0;
      cout << "x.norm() = " << x.norm() << endl;
      dog.update_V_vector(x);
      cout << "Rolled back and is_folded = " << is_folded() << endl;
      p.fold_bias_weight *= 10;
      obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.dihedral_weight,p.fold_bias_weight,p.wallpaper_curve_weight,p.wallpaper_curve_weight});
      newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
      dog.update_V_vector(x);
    }
  }
  cout << "Finished: fold bias = " << p.fold_bias_weight << " and is_folded = " << is_folded() << endl;

  dog.update_V_vector(x);
  
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  objective = obj.compObj.obj(x);
}

void DogSolver::single_iteration_normal(double& constraints_deviation, double& objective) {
  cout << "running a single optimization routine" << endl;
  Eigen::VectorXd x0(dog.getV_vector()), x(x0);

  obj.compObj.update_weights({p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, 0.1*p.soft_pos_weight, p.dihedral_weight, p.dihedral_weight,p.fold_bias_weight, p.wallpaper_curve_weight,p.wallpaper_curve_weight});

  newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
  is_folded();

  dog.update_V_vector(x);
  
  objective = obj.compObj.obj(x);
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  if (time_measurements_log) {
    *time_measurements_log << objective << "," << constraints_deviation << endl;
  }
}