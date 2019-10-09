#include "DogSolver.h"

#include <igl/polar_svd.h>

using namespace std;

DogSolver::DogSolver(Dog& dog, const Eigen::VectorXd& init_mesh_vars,
        DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
        std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
        std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
        std::vector<std::pair<int,int>>& pairs,
        std::vector<pair<int,int>>& bnd_vertices_pairs,
        std::ofstream* time_measurements_log) :

          dog(dog), x(init_variables(init_mesh_vars)),
          foldingBinormalBiasConstraints(dog),
          foldingMVBiasConstraints(dog, p.flip_sign), p(p),
          constraints(dog, x, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles, mvTangentCreaseAngleParams,
                      mv_cos_angles, pairs),
          obj(dog, x, constraints, foldingBinormalBiasConstraints, foldingMVBiasConstraints, bnd_vertices_pairs, p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          time_measurements_log(time_measurements_log)
           {
    is_constrained = (b.rows() + edgePoints.size())>0;
    if (time_measurements_log) {
      p.max_newton_iters = 1;
    }

    /*
    // output to find out how stuff works
    int v_num = dog.get_v_num();
    cout << "Dog has " << v_num << " vertices [pos]\n";
    cout << "dogsolver has x of size " << x.size() << " (and that should be 3x as much as the dog.v_num) [pos]\n";
    cout << "position constraints of the new dogSolver:\n";
    if(b.size()>0){
      for(int i=0;i<b.size();++i)cout << "pos\t" << b[i] << "\t" << bc[i] << "\tin submesh " << dog.v_to_submesh_idx(b[i]%v_num) << endl;
    }
    cout << "edge constraints of the new dogSolver: "<<edgePoints.size()<<"\n";
    cout << "paired vertices of the new dogSolver:\n";
    if(pairs.size()>0){
      for(int i=0;i<pairs.size();++i)cout << "pai\t" << pairs[i].first << "\t" << pairs[i].second << endl;
    }
    */

    int num_submeshes = dog.get_submesh_n();
    cout << "DogSolver constructor: There are " << num_submeshes << " submeshes\n";
    for(int i=0; i<num_submeshes; ++i){
      int mini,maxi;
      dog.get_submesh_min_max_i(i, mini, maxi, true);
      cout << "\tsubmesh " << i << " has (mini, maxi) (" << mini << ", " << maxi << ")\n";
    }

    if(num_submeshes>1){
      sub_dog.resize(num_submeshes);
      sub_dogsolver.resize(num_submeshes);

      int num_v = dog.get_v_num();

      // Absolutely no constraints for now
      empty_xi.resize(0);
      empty_xd.resize(0);
      empty_ep.clear();
      empty_mat.resize(0,3);
      empty_egg.clear();
      empty_d.clear();
      empty_thing.clear();
      empty_pair.clear();

      const DogEdgeStitching& edgeStitching = dog.getEdgeStitching();
      for(int j=0;j<edgeStitching.stitched_curves.size(); ++j){
        for(int k=0;k<edgeStitching.stitched_curves[j].size(); ++k){
          cout << "edgeStitching.stitched_curves["<<j<<"]["<<k<<"] is\t" << edgeStitching.stitched_curves[j][k].edge.v1 <<"-"<<edgeStitching.stitched_curves[j][k].edge.v2 << " at " << edgeStitching.stitched_curves[j][k].t<<endl;
        }
      }
      for(int i=0; i<num_v; ++i){
        cout << "v_to_submesh_idx("<<i<<") = " << dog.v_to_submesh_idx(i)<<endl;// gives index of submesh
        cout << "v_in_submesh("<<i<<") = "<<dog.v_in_submesh(i)<<endl;//gives index of v in its submesh
      }
      vector< vector<int> > stored_stitch(num_submeshes);
      constrained_edge_points.resize(num_submeshes);
      for(int i=0; i<num_submeshes; ++i) {constrained_edge_points[i].clear(); stored_stitch[i].clear();}
      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v11 = edgeStitching.edge_const_1[j].v1;
        int v12 = edgeStitching.edge_const_1[j].v2;
        int v21 = edgeStitching.edge_const_2[j].v1;
        int v22 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];

        cout << "edgeStitching " << j << ":\t" << v11 << "-" << v12 << " corresponding to " << v21 << "-" << v22
            << " at coord " << t <<endl;

        int submesh_1 = dog.v_to_submesh_idx(v11);
        int submesh_2 = dog.v_to_submesh_idx(v21);
        constrained_edge_points[submesh_1].push_back(EdgePoint(Edge(dog.v_in_submesh(v11), dog.v_in_submesh(v12)),t));
        constrained_edge_points[submesh_2].push_back(EdgePoint(Edge(dog.v_in_submesh(v21), dog.v_in_submesh(v22)),t));//1-t?

        stored_stitch[submesh_1].push_back(j);
        stored_stitch[submesh_2].push_back(j);//.... :/
      }
      for(int i=0;i<constrained_edge_points.size();++i){
        for(int j=0;j<constrained_edge_points[i].size();++j){
          cout << "ij["<<i<<"]["<<j<<"] is\t" << constrained_edge_points[i][j].edge.v1<<endl;
        }
      }


      for(int i=0; i<num_submeshes; ++i){
        sub_dog[i] = dog.get_submesh(i);
        Eigen::VectorXd sub_x = sub_dog[i]->getV_vector();

        /*
        // PositionalConstraints for submeshes
        vector<int> sb(0);
        vector<double> sbc(0);
        int mini,maxi;
        dog.get_submesh_min_max_i(i, mini, maxi, true);
        for(int j=0; j<b.size(); ++j){
          if(dog.v_to_submesh_idx(b[j]%num_v) == i){
            sb.push_back( b[j]%num_v - mini + (b[j]/num_v)*sub_dog[i]->get_v_num() );//TODO right is?
            sbc.push_back(bc[j]);
            cout << "b[j]:\t" << b[j] << "\t" << bc[j] <<endl;
            cout << "sb:\t" << sb[sb.size()-1] << "\t" << sbc[sbc.size()-1]<<endl;
          }
        }
        Eigen::VectorXi sub_b(sb.size());
        Eigen::VectorXd sub_bc(sb.size());
        for(int j; j<sub_b.size(); ++j){
          sub_b[j] = sb[j];
          sub_bc[j] = sbc[j];
        }
        */

        Eigen::MatrixXd sub_m = EdgePoint::getPositionInMesh(constrained_edge_points[i], sub_dog[i]->getV());
        /*
        for(int j=0; j<sub_m.rows(); ++j){
          sub_m.rows(j) << constrained_edge_points[i][j].
        }
        */

        cout << "constructing subsolver "<<i<<endl;
        sub_dogsolver[i] = new DogSolver(*sub_dog[i], sub_x, p,
              //sub_b,sub_bc,
              empty_xi, empty_xd,
              constrained_edge_points[i], sub_m,
              //empty_ep, empty_mat,
              empty_egg, empty_d,
              empty_thing, empty_d,
              empty_pair,
              empty_pair,
              NULL);
      }
      cout << "Sub_dogsolvers constructed\n";
    } else {
      sub_dog.clear();
      sub_dogsolver.clear();
    }
}

DogSolver::~DogSolver(){
  if(sub_dog.size() > 0){
    for(int i=0; i<sub_dog.size(); ++i) delete sub_dog[i];
  }
  if(sub_dogsolver.size() > 0){
    for(int i=0; i<sub_dogsolver.size(); ++i) delete sub_dogsolver[i];
  }
}

Eigen::VectorXd DogSolver::init_variables(const Eigen::VectorXd& init_mesh_vars) {
  int dog_vars_num = init_mesh_vars.rows();
  int var_num = dog_vars_num;
  // Here we can increase var_num and add more optimization variables at the end of x, if needed for other stuff (like symmetry for instance)
  Eigen::VectorXd x(var_num); x.setZero();
  x.head(dog_vars_num) = init_mesh_vars;
  return x;
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
                    edgePtConst(dog.getQuadTopology(),edgePoints, edgeCoords),
                    edgeAngleConst(dog.getQuadTopology(), init_x0, edge_angle_pairs, edge_cos_angles),
                    mvTangentCreaseAngleConst(dog.getQuadTopology(), init_x0, mvTangentCreaseAngleParams, mv_cos_angles),
                    ptPairConst(pairs),
                    //compConst({&dogConst, &stitchingConstraints}) /*the positional/edge constraints are soft and go into the objective*/ {
                    compConst({&dogConst}) /*the positional/edge constraints are soft and go into the objective*/ {cout<<"Constraints done\n";
    // Empty on purpose
}

DogSolver::Objectives::Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
        Constraints& constraints,/*
          PositionalConstraints& posConst,
          EdgePointConstraints& edgePtConst,
          EdgesAngleConstraints& edgeAnglesConst,
          PointPairConstraints& ptPairConst,*/
          FoldingBinormalBiasConstraints& foldingBinormalBiasConstraints,
          FoldingMVBiasConstraints& foldingMVBiasConstraints,
          std::vector<std::pair<int,int>>& bnd_vertices_pairs,
          const DogSolver::Params& p) :
        bending(dog.getQuadTopology(), init_x0), isoObj(dog.getQuadTopology(), init_x0),
        pointsPosSoftConstraints(constraints.posConst, init_x0),
        edgePosSoftConstraints(constraints.edgePtConst, init_x0),
        edgeAnglesSoftConstraints(constraints.edgeAngleConst, init_x0),
        mvTangentCreaseSoftConstraints(constraints.mvTangentCreaseAngleConst, init_x0),
        ptPairSoftConst(constraints.ptPairConst, init_x0),
        foldingBinormalBiasObj(foldingBinormalBiasConstraints, init_x0),
        foldingMVBiasObj(foldingMVBiasConstraints, init_x0),
        stitchingConstraintsPenalty(constraints.stitchingConstraints, init_x0),
        pairedBndVertBendingObj(dog.getQuadTopology(), bnd_vertices_pairs, init_x0, Eigen::Vector3d(1,0,0)),
        //allConstQuadraticObj(constraints, init_x0),
        /*curvedFoldingBiasObj(curvedFoldingBiasObj),*/

        compObj(
          {&bending, &isoObj, &stitchingConstraintsPenalty, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst, &edgeAnglesSoftConstraints, &mvTangentCreaseSoftConstraints, &foldingBinormalBiasObj, &foldingMVBiasObj,&pairedBndVertBendingObj},
          {p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.stitching_weight,p.soft_pos_weight, p.soft_pos_weight, p.pair_weight, p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight, p.paired_boundary_bending_weight})
          /*
        compObj(
          {&bending, &isoObj, &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst, &edgeAnglesSoftConstraints, &foldingBinormalBiasObj, &allConstQuadraticObj},
          {p.bending_weight,p.isometry_weight, p.soft_pos_weight, p.soft_pos_weight, pair_weight, p.dihedral_weight, p.fold_bias_weight,1})
          */
          {cout << "objectives done\n";
    // Empty on purpose
  //std::cout << "p.isometry_weight/dog.getQuadTopology().E.rows() = " << p.isometry_weight/dog.getQuadTopology().E.rows() << std::endl; exit(1);
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
      if ( (abs((ep_p-ep_b_p).normalized().dot(e1.normalized())) > 0.9995) || (abs((ep_p-ep_b_p).normalized().dot(e1.normalized()))>0.9995) ) {std::cout << "e_idx = " << e_idx << std::endl; continue;}
      if ( (edge_t < 0.05 ) || ( edge_t > 0.95) ) continue;
      //if ( ((ep_p-ep_b_p).normalized().dot(e1.normalized()) < 0.1) || ((ep_p-ep_b_p).normalized().dot(e1.normalized())<0.1) ) {std::cout << "e_idx = " << e_idx << std::endl; continue;}
      if ( sign1*sign2 > 0) {
        is_folded = false;
        //cout << "Change!" << endl;

        cout << "edge_t = " << edge_t << std::endl;
        cout << "(ep_p-ep_b_p).norm() = " << (ep_p-ep_b_p).norm()  << std::endl;
        cout << "(ep_f_p-ep_p).norm() = " << (ep_f_p-ep_p).norm()  << std::endl;
        cout << "(ep_p-ep_b_p).normalized().dot(e1.normalized()) " << (ep_p-ep_b_p).normalized().dot(e1.normalized()) << std::endl;
        cout << "(ep_f_p-ep_p).normalized().dot(e1.normalized()) " << (ep_f_p-ep_p).normalized().dot(e1.normalized()) << std::endl;
        cout << "(ep_p-ep_b_p).normalized().dot(e2.normalized()) " << (ep_p-ep_b_p).normalized().dot(e2.normalized()) << std::endl;
        cout << "(ep_f_p-ep_p).normalized().dot(e2.normalized()) " << (ep_f_p-ep_p).normalized().dot(e2.normalized()) << std::endl;
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

bool DogSolver::is_mountain_valley_correct(const Eigen::VectorXd& x) {
  return foldingMVBiasConstraints.is_mv_assignment_correct(x);
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
  if(mode_subsolvers) return single_iteration_subsolvers(constraints_deviation, objective);
  if (!p.folding_mode) p.fold_bias_weight = 0;
  if (p.folding_mode) return single_iteration_fold(constraints_deviation, objective);
  else return single_iteration_normal(constraints_deviation, objective);
}

void DogSolver::single_iteration_fold(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(x);
  if (!is_folded()) {
    cout << "Error: Not folded" << endl; exit(1);
  }

  obj.compObj.update_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight, p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,p.paired_boundary_bending_weight});
  newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
  dog.update_V_vector(x.head(3*dog.get_v_num()));

  if (is_folded()) {
    p.fold_bias_weight = 1;
  } else {
    while (!is_folded() && (p.fold_bias_weight < 1e14)) {
      cout << "Not folded, fold bias = " << p.fold_bias_weight << " reverting back and making the bias stronger" << endl;
      x = x0;
      cout << "x.norm() = " << x.norm() << endl;
      dog.update_V_vector(x.head(3*dog.get_v_num()));
      cout << "Rolled back and is_folded = " << is_folded() << endl;
      p.fold_bias_weight *= 10;
      obj.compObj.update_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight, p.dihedral_weight, p.dihedral_weight,p.fold_bias_weight,p.mv_bias_weight,p.paired_boundary_bending_weight});
      newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
      dog.update_V_vector(x.head(3*dog.get_v_num()));
    }
  }
  cout << "Finished: fold bias = " << p.fold_bias_weight << " and is_folded = " << is_folded() << endl;

  dog.update_V_vector(x.head(3*dog.get_v_num()));

  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  objective = obj.compObj.obj(x);

  if (p.mv_bias_weight) {
    if (!is_mountain_valley_correct(x)){
      std::cout << std::endl << std::endl << "M/V NOT CORRECT!" << std::endl;
    } else {
      std::cout << std::endl << std::endl << "M/V IS OK!" << std::endl;
    }
  }
}

void DogSolver::single_iteration_subsolvers(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine" << endl;
	Eigen::VectorXd x0(x);
  //TODO update constrained edge points ...?
  bool use_subsolvers = (sub_dog.size()>1);
  if(!use_subsolvers){
    cout << " no subsolvers\n";
    obj.compObj.update_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight, p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,p.paired_boundary_bending_weight});
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";
      Eigen::MatrixXd sub_m = EdgePoint::getPositionInMesh(constrained_edge_points[i], sub_dog[i]->getV());
      sub_dogsolver[i]->update_edge_coords(sub_m);

      sub_dogsolver[i]->set_opt_vars(sub_dog[i]->getV_vector());
      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_normal(sub_cd, sub_obj);
      cout << " subsolver " << i << " done\n";
    }

    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
  	   dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    cout << "all subsolvers done\n";
  }
}

void DogSolver::single_iteration_normal(double& constraints_deviation, double& objective) {
  //p.fold_bias_weight = 0;
  cout << "running a single optimization routine" << endl;
  Eigen::VectorXd x0(x);

  obj.compObj.update_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight, p.dihedral_weight, p.dihedral_weight,p.fold_bias_weight, p.mv_bias_weight,p.paired_boundary_bending_weight});
  newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
  is_folded();

  dog.update_V_vector(x.head(3*dog.get_v_num()));

  objective = obj.compObj.obj(x);
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  if (time_measurements_log) {
    *time_measurements_log << objective << "," << constraints_deviation << endl;
  }
}

void DogSolver::get_x_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T) {
  int vn = 3*dog.get_v_num();
  T << x(vn),x(vn+1),x(vn+2);
  R.row(0) << x(vn+3),x(vn+4),x(vn+5);
  R.row(1) << x(vn+6),x(vn+7),x(vn+8);
  R.row(2) << x(vn+9),x(vn+10),x(vn+11);
}

void DogSolver::set_x_rotation(Eigen::Matrix3d& R) {
  int vn = 3*dog.get_v_num();
  x(vn+3) = R(0,0); x(vn+4) = R(0,1); x(vn+5) = R(0,2);
  x(vn+6) = R(1,0); x(vn+7) = R(1,1); x(vn+8) = R(1,2);
  x(vn+9) = R(2,0); x(vn+10) = R(2,1); x(vn+11) = R(2,2);
}

void DogSolver::set_y_rotation(Eigen::Matrix3d& R) {
  int vn = 3*dog.get_v_num()+12;
  x(vn+3) = R(0,0); x(vn+4) = R(0,1); x(vn+5) = R(0,2);
  x(vn+6) = R(1,0); x(vn+7) = R(1,1); x(vn+8) = R(1,2);
  x(vn+9) = R(2,0); x(vn+10) = R(2,1); x(vn+11) = R(2,2);
}

void DogSolver::get_y_rigid_motion(Eigen::Matrix3d& R, Eigen::RowVector3d& T) {
  int vn = 3*dog.get_v_num()+12;
  T << x(vn),x(vn+1),x(vn+2);
  R.row(0) << x(vn+3),x(vn+4),x(vn+5);
  R.row(1) << x(vn+6),x(vn+7),x(vn+8);
  R.row(2) << x(vn+9),x(vn+10),x(vn+11);
}
