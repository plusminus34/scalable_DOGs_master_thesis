#include "DogSolver.h"

#include <igl/polar_svd.h>
#include <igl/procrustes.h>

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

          dog(dog), x(dog.getV_vector()),
          foldingBinormalBiasConstraints(dog),
          foldingMVBiasConstraints(dog, p.flip_sign), p(p),
          constraints(dog, init_mesh_vars, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles, mvTangentCreaseAngleParams,
                      mv_cos_angles, pairs),
          obj(dog, init_mesh_vars, constraints, foldingBinormalBiasConstraints, foldingMVBiasConstraints, bnd_vertices_pairs, p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          time_measurements_log(time_measurements_log)
           {
    is_constrained = (b.rows() + edgePoints.size())>0;
    if (time_measurements_log) {
      p.max_newton_iters = 1;
    }

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

      int v_num = dog.get_v_num();

      for(int i=0; i<num_submeshes; ++i) sub_dog[i] = dog.get_submesh(i);

      // Mostly no constraints for now
      empty_xi.resize(0);
      empty_xd.resize(0);
      empty_ep.clear();
      empty_mat.resize(0,3);
      empty_egg.clear();
      empty_d.clear();
      empty_thing.clear();
      empty_pair.clear();

      //Subsolver positional constraints from global positional constraints
      {
        sub_b.resize(num_submeshes);
        sub_bc.resize(num_submeshes);
        sub_ij_to_bc.resize(num_submeshes);
        for(int i=0; i<num_submeshes; ++i) {
          sub_b[i].resize(0);
        }
        for(int i=0; i<b.size()/3; ++i){
          int submesh_i = dog.v_to_submesh_idx(b[i]);
          sub_b[submesh_i].resize(sub_b[submesh_i].size()+3);
        }
        for(int i=0; i<num_submeshes; ++i) {
          sub_bc[i].resize(sub_b[i].size());
          sub_ij_to_bc[i].resize(sub_b[i].size());
        }
        vector<int> current_j(num_submeshes,0);
        for(int i=0; i<b.size()/3; ++i){
          int submesh_i = dog.v_to_submesh_idx(b[i]);
          int size_i = sub_b[submesh_i].size() / 3;
          int local_idx = dog.v_in_submesh(b[i]);
          //cout << "writing element " << i << " into [" << submesh_i << "][" << current_j[submesh_i] << ", "<< (current_j[submesh_i] + size_i)<<", "<<(current_j[submesh_i] + 2*size_i)<<"] where sub_b has size "<<sub_b[submesh_i].size()<<"\n";
          sub_b[submesh_i][current_j[submesh_i] ] = local_idx;
          sub_b[submesh_i][ current_j[submesh_i] + size_i ] = local_idx + sub_dog[submesh_i]->get_v_num();
          sub_b[submesh_i][current_j[submesh_i] + 2*size_i ] = local_idx + 2*sub_dog[submesh_i]->get_v_num();
          sub_bc[submesh_i][current_j[submesh_i] ] = bc[i];
          sub_bc[submesh_i][current_j[submesh_i] + size_i ] = bc[i + b.size()/3] ;
          sub_bc[submesh_i][current_j[submesh_i] + 2*size_i ] = bc[i + 2*b.size()/3];
          sub_ij_to_bc[submesh_i][ current_j[submesh_i] ] = i;
          sub_ij_to_bc[submesh_i][ current_j[submesh_i] + size_i ] = i + b.size()/3;
          sub_ij_to_bc[submesh_i][ current_j[submesh_i] + 2*size_i ] = i + 2*b.size()/3;
          ++current_j[submesh_i];
        }
      }

      //Subsolver edge point constraints from stitching
      const DogEdgeStitching& edgeStitching = dog.getEdgeStitching();
      constrained_edge_points.resize(num_submeshes);
      sub_edgeCoords.resize(num_submeshes);
      corresponding_edge_points.resize(num_submeshes);

      for(int i=0; i<num_submeshes; ++i) {
        constrained_edge_points[i].clear();
        corresponding_edge_points[i].clear();
      }

      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v11 = edgeStitching.edge_const_1[j].v1;
        int v12 = edgeStitching.edge_const_1[j].v2;
        int v21 = edgeStitching.edge_const_2[j].v1;
        int v22 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];

        int submesh_1 = dog.v_to_submesh_idx(v11);
        int submesh_2 = dog.v_to_submesh_idx(v21);

        corresponding_edge_points[submesh_2].push_back(pair<int,int>(submesh_1, constrained_edge_points[submesh_1].size()));
        corresponding_edge_points[submesh_1].push_back(pair<int,int>(submesh_2, constrained_edge_points[submesh_2].size()));

        constrained_edge_points[submesh_1].push_back(EdgePoint(Edge(dog.v_in_submesh(v11), dog.v_in_submesh(v12)),t));
        constrained_edge_points[submesh_2].push_back(EdgePoint(Edge(dog.v_in_submesh(v21), dog.v_in_submesh(v22)),t));
      }

      // initialize sub_edgeCoords
      for(int i=0; i<num_submeshes; ++i) {
        sub_edgeCoords[i] = Eigen::MatrixXd(constrained_edge_points[i].size(), 3);
      }
      vector<int> current_row(num_submeshes,0);
      for(int i=0; i<edgeStitching.edge_coordinates.size(); ++i){
        int v11 = edgeStitching.edge_const_1[i].v1;
        int v12 = edgeStitching.edge_const_1[i].v2;
        int v21 = edgeStitching.edge_const_2[i].v1;
        int v22 = edgeStitching.edge_const_2[i].v2;
        double t = edgeStitching.edge_coordinates[i];
        int submesh_1 = dog.v_to_submesh_idx(v11);
        int submesh_2 = dog.v_to_submesh_idx(v21);
        sub_edgeCoords[submesh_1].row(current_row[submesh_1]) =
          EdgePoint(Edge(v21, v22),t).getPositionInMesh(dog.getV());
        sub_edgeCoords[submesh_2].row(current_row[submesh_2]) =
          EdgePoint(Edge(v11, v12),t).getPositionInMesh(dog.getV());
        current_row[submesh_1]++;
        current_row[submesh_2]++;
      }

      is_constrained = (is_constrained || edgeStitching.edge_coordinates.size()>0);


      std::vector< Eigen::SparseMatrix<double> > sub_admm_A(num_submeshes);

      for(int i=0; i<num_submeshes; ++i){
        sub_admm_A[i] = Eigen::SparseMatrix<double>(edgeStitching.edge_coordinates.size(), 3*sub_dog[i]->get_v_num());
        sub_admm_A[i].reserve(6*sub_admm_A[i].rows());
      }

      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v11 = edgeStitching.edge_const_1[j].v1;
        int v12 = edgeStitching.edge_const_1[j].v2;
        int v21 = edgeStitching.edge_const_2[j].v1;
        int v22 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];
        int submesh_1 = dog.v_to_submesh_idx(v11);
        int submesh_2 = dog.v_to_submesh_idx(v21);

        int subsize = dog.get_submesh_i_size(submesh_1);
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v11)) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v11) + subsize) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v11) + 2*subsize) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v12)) = 1.0-t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v12) + subsize) = 1.0-t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v12) + 2*subsize) = 1.0-t;

        subsize = dog.get_submesh_i_size(submesh_2);
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v21)) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v21) + subsize) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v21) + 2*subsize) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v22)) = t-1.0;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v22) + subsize) = t-1.0;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(v22) + 2*subsize) = t-1.0;
      }

      //edge pairs for the subsolvers
      vector< vector< pair<Edge,Edge> > > sub_egg(num_submeshes);
      vector< vector<double> > sub_angles(num_submeshes);
      vector< Eigen::MatrixXd > sub_outsidepts(num_submeshes);
      {
        sub_edge_angle_ww.resize(num_submeshes);
        sub_idx_to_angle_idx.resize(num_submeshes);
        for(int i=0; i<num_submeshes; ++i){
          sub_egg[i].clear();
          sub_angles[i].clear();
          sub_outsidepts[i].resize(0,6);
          sub_idx_to_angle_idx[i].resize(0);
        }
        for(int i=0; i<edge_angle_pairs.size(); ++i){
          int v1 = edge_angle_pairs[i].first.v1;
          int v2 = edge_angle_pairs[i].first.v2;
          int w1 = edge_angle_pairs[i].second.v1;
          int w2 = edge_angle_pairs[i].second.v2;
          int submesh_i = dog.v_to_submesh_idx(v1);
          int sub_v1 = dog.v_in_submesh(v1);
          int sub_v2 = dog.v_in_submesh(v2);
          sub_egg[submesh_i].push_back(pair<Edge,Edge>(Edge(sub_v1, sub_v2),Edge(sub_v1, sub_v2)));
          sub_angles[submesh_i].push_back(edge_cos_angles[i]);


          Eigen::MatrixXi rep1(sub_edge_angle_ww[submesh_i].rows()+1, 2);
          rep1 << sub_edge_angle_ww[submesh_i], w1, w2;
          sub_edge_angle_ww[submesh_i] = rep1;
          sub_idx_to_angle_idx[submesh_i].push_back(i);

          submesh_i = dog.v_to_submesh_idx(w1);
          sub_v1 = dog.v_in_submesh(w1);
          sub_v2 = dog.v_in_submesh(w2);
          sub_egg[submesh_i].push_back(pair<Edge,Edge>(Edge(sub_v1, sub_v2),Edge(sub_v1, sub_v2)));
          sub_angles[submesh_i].push_back(edge_cos_angles[i]);

          Eigen::MatrixXi rep2(sub_edge_angle_ww[submesh_i].rows()+1, 2);
          rep2 << sub_edge_angle_ww[submesh_i], v1, v2;
          sub_edge_angle_ww[submesh_i] = rep2;
          sub_idx_to_angle_idx[submesh_i].push_back(i);
        }
      }

      //2-patch Procrustes stuff
      if(num_submeshes==2){
        proc_T.resize(2);
        int num_edgepts = edgeStitching.edge_coordinates.size();
        proc_T[0].resize(num_edgepts, sub_dog[0]->get_v_num());
        proc_T[1].resize(num_edgepts, sub_dog[1]->get_v_num());
        proc_T[0].reserve(2*num_edgepts);
        proc_T[1].reserve(2*num_edgepts);
        for(int i=0; i<num_edgepts; ++i){
          int v11 = edgeStitching.edge_const_1[i].v1;
          int v12 = edgeStitching.edge_const_1[i].v2;
          int v21 = edgeStitching.edge_const_2[i].v1;
          int v22 = edgeStitching.edge_const_2[i].v2;
          double t = edgeStitching.edge_coordinates[i];
          proc_T[0].insert(i, dog.v_in_submesh(v11)) = t;
          proc_T[0].insert(i, dog.v_in_submesh(v12)) = 1.0-t;
          proc_T[1].insert(i, dog.v_in_submesh(v21)) = t;
          proc_T[1].insert(i, dog.v_in_submesh(v22)) = 1.0-t;
        }
      }

      for(int i=0; i<num_submeshes; ++i){
        int mini, maxi, sub_v_num;
        sub_v_num = sub_dog[i]->get_v_num();
        dog.get_submesh_min_max_i(i, mini, maxi);

        Eigen::VectorXd sub_x0(3*sub_v_num);
        for(int j=0; j<sub_v_num; ++j){
          sub_x0[j] = init_mesh_vars[mini + j];
          sub_x0[j + sub_v_num] = init_mesh_vars[mini + v_num + j];
          sub_x0[j + 2*sub_v_num] = init_mesh_vars[mini + 2*v_num + j];
        }

        cout << "constructing subsolver "<<i<<endl;
        sub_dogsolver[i] = new DogSolver(*sub_dog[i], sub_x0, p,
              sub_b[i], sub_bc[i],
              //empty_xi, empty_xd,
              constrained_edge_points[i], sub_edgeCoords[i],
              //empty_ep, empty_mat,
              sub_egg[i], sub_angles[i],
              //empty_egg, empty_d,
              empty_thing, empty_d,
              empty_pair,
              empty_pair,
              NULL);

        sub_dogsolver[i]->set_lambda(Eigen::VectorXd::Zero(sub_admm_A[i].rows()));
        sub_dogsolver[i]->set_z(Eigen::VectorXd::Zero(sub_admm_A[i].rows()));
        sub_dogsolver[i]->build_VSADMMObjective(sub_admm_A[i]);

        Eigen::SparseMatrix<double> I_i(sub_x0.size(), sub_x0.size());
        I_i.setIdentity(); I_i *= 0.1*(num_submeshes-1) *0.01;
        sub_dogsolver[i]->build_ProximalObjective(I_i);
        sub_dogsolver[i]->remake_compobj();
        cout << " constructed subsolver "<<i<<endl;
      }
      cout << "Sub dogsolvers constructed\n";
    } else {
      sub_dog.clear();
      sub_dogsolver.clear();
    }

}

DogSolver::~DogSolver(){
  if(vsadmm_obj) delete vsadmm_obj;
  if(sub_dog.size() > 0){
    for(int i=0; i<sub_dog.size(); ++i) delete sub_dog[i];
  }
  if(sub_dogsolver.size() > 0){
    for(int i=0; i<sub_dogsolver.size(); ++i) delete sub_dogsolver[i];
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
                    edgePtConst(dog.getQuadTopology(),edgePoints, edgeCoords),
                    edgeAngleConst(dog.getQuadTopology(), init_x0, edge_angle_pairs, edge_cos_angles),
                    subEdgesAngleConst(dog.get_v_num(), init_x0, edge_angle_pairs, edge_cos_angles),
                    mvTangentCreaseAngleConst(dog.getQuadTopology(), init_x0, mvTangentCreaseAngleParams, mv_cos_angles),
                    ptPairConst(pairs),
                    //compConst({&dogConst, &stitchingConstraints}) /*the positional/edge constraints are soft and go into the objective*/ {
                    compConst({&dogConst}) /*the positional/edge constraints are soft and go into the objective*/ {
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
          {
            if (dog.get_submesh_n()==1 ){//subsolvers use different edge angle constraints
              edgeAnglesSoftConstraints = QuadraticConstraintsSumObjective(constraints.subEdgesAngleConst, init_x0);
            }
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
  if(mode == mode_subsolvers) {
    single_iteration_subsolvers(constraints_deviation, objective);
  } else if (mode == mode_vsadmm || mode == mode_jadmm || mode == mode_proxjadmm) {
    single_iteration_ADMM(constraints_deviation, objective);
  } else if (mode == mode_serial) {
    single_iteration_serial(constraints_deviation, objective);
  } else if (mode == mode_procrustes) {
    single_iteration_procrustes(constraints_deviation, objective);
  } else if (mode == mode_experimental) {
    single_iteration_experimental(constraints_deviation, objective);
  } else {
    // mode == mode_standard
    if (!p.folding_mode) p.fold_bias_weight = 0;
    if (p.folding_mode) single_iteration_fold(constraints_deviation, objective);
    else single_iteration_normal(constraints_deviation, objective);
  }
}

void DogSolver::single_iteration_fold(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (fold)" << endl;
	x0 = x;
  if (!is_folded()) {
    cout << "Error: Not folded" << endl;
    exit(1);
  }
  update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
    p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
    p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
    p.paired_boundary_bending_weight, 0, 0, 0});
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
      update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
        p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
        p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
        p.paired_boundary_bending_weight, 0, 0, 0});
      newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
      dog.update_V_vector(x.head(3*dog.get_v_num()));
    }
  }
  cout << "Finished: fold bias = " << p.fold_bias_weight << " and is_folded = " << is_folded() << endl;

  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  objective = obj.compObj.obj(x);

  if (p.mv_bias_weight) {
    if (!is_mountain_valley_correct(x)){
      std::cout << std::endl << std::endl << "M/V NOT CORRECT!" << std::endl;
    } else {
      std::cout << std::endl << std::endl << "M/V IS OK!" << std::endl;
    }
  }
  ++iter_i;
}

void DogSolver::single_iteration_subsolvers(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (subsolvers)" << endl;
	x0 = x;
  if(is_subsolver()){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();

    constraints_deviation = 0.0;
    objective = 0.0;

    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0, 0, 0});

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";

      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);

      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_subsolvers(sub_cd, sub_obj);
      constraints_deviation += sub_cd;
      objective += sub_obj;

      cout << " subsolver " << i << " done\n";
    }

    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
      dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    x = dog.getV_vector();

    update_sub_edgeCoords();

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_ADMM(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (ADMM)" << endl;
	x0 = x;
  if(is_subsolver()){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();

    double rho2 = (mode==mode_proxjadmm ? p.admm_rho : 0);
    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0.5*p.admm_rho, rho2, 0});

    constraints_deviation = 0.0;
    objective = 0.0;

    vector<Eigen::VectorXd> sub_Ax(num_submeshes);
    Eigen::VectorXd sum_Ax, sum_lambda;
    vector<Eigen::VectorXd> sub_z(num_submeshes);
    for(int i=0; i<num_submeshes; ++i){
      sub_Ax[i] = sub_dogsolver[i]->get_Ax();
      if (i==0) sum_Ax = sub_Ax[i];
      else sum_Ax += sub_Ax[i];
      if(mode == mode_vsadmm){
        if (i==0) sum_lambda = sub_dogsolver[i]->get_lambda();
        else sum_lambda += sub_dogsolver[i]->get_lambda();
      }
    }
    for(int i=0; i<num_submeshes; ++i){
      if (mode == mode_vsadmm) {
        sub_z[i] = sub_Ax[i] - 1.0/num_submeshes * (sum_Ax - sum_lambda/p.admm_rho);
      } else if (mode == mode_jadmm || mode == mode_proxjadmm) {
        sub_z[i] = sub_Ax[i] - sum_Ax + sub_dogsolver[0]->get_lambda()/p.admm_rho;
      }
      sub_dogsolver[i]->set_z(sub_z[i]);
    }

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";

      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);

      double sub_cd = constraints_deviation;
      double sub_obj = objective;

      sub_dogsolver[i]->single_iteration_ADMM(sub_cd, sub_obj);

      constraints_deviation += sub_cd;
      objective += sub_obj;

      if(mode == mode_vsadmm){
        Eigen::VectorXd sub_lambda = sub_dogsolver[i]->get_lambda();
        sub_lambda -= p.admm_rho*(sub_dogsolver[i]->get_Ax() - sub_z[i] + sub_lambda/p.admm_rho);
        sub_dogsolver[i]->set_lambda(sub_lambda);
      }

      cout << " subsolver " << i << " done\n";
    }

    if(mode == mode_jadmm || mode == mode_proxjadmm){
      for(int i=0; i<num_submeshes; ++i){
        sub_Ax[i] = sub_dogsolver[i]->get_Ax();
        if(i==0) sum_Ax = sub_Ax[i];
        else sum_Ax += sub_Ax[i];
      }
      Eigen::VectorXd lambda = sub_dogsolver[0]->get_lambda() - p.admm_gamma * p.admm_rho * sum_Ax;
      for(int i=0; i<num_submeshes; ++i) sub_dogsolver[i]->set_lambda(lambda);
    }

    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
  	   dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    x = dog.getV_vector();

    update_sub_edgeCoords();

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_serial(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (serial)" << endl;
	x0 = x;
  if(is_subsolver()){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();

    update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0.5*p.admm_rho, 0.0, p.admm_gamma});

    Eigen::VectorXd sum_Ax = sub_dogsolver[0]->get_Ax();
    for(int i=1; i<num_submeshes; ++i){
      sum_Ax += sub_dogsolver[i]->get_Ax();
    }

    constraints_deviation = 0.0;
    objective = 0.0;

    //Assumption: patches are already ordered
    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";

      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);
      sum_Ax -= sub_dogsolver[i]->get_Ax();
      sub_dogsolver[i]->set_z(-sum_Ax);

      double sub_cd = constraints_deviation;
      double sub_obj = objective;

      sub_dogsolver[i]->single_iteration_serial(sub_cd, sub_obj);

      constraints_deviation += sub_cd;
      objective += sub_obj;

      dog.update_submesh_V(i, sub_dog[i]->getV());
      x = dog.getV_vector();
      sum_Ax += sub_dogsolver[i]->get_Ax();

      update_sub_edgeCoords();

      cout << " subsolver " << i << " done\n";
    }

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_procrustes(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (procrustes)" << endl;
	x0 = x;
  if(is_subsolver()){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();
    if(num_submeshes != 2) cout <<"WARNING: You shouldn't use this solver mode for more than 2 patches\n";

    update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0.5*p.admm_rho, 0, p.admm_rho});
      //0,p.admm_rho,0});
      //LinearConstraints, ProximalObjective, SomeSerialObjective});

    constraints_deviation = 0.0;
    objective = 0.0;

    cout << "solving submesh 0\n";
    sub_dogsolver[0]->set_z(-sub_dogsolver[1]->get_Ax());
    sub_dogsolver[0]->single_iteration_procrustes(constraints_deviation, objective);

    dog.update_submesh_V(0, sub_dog[0]->getV());
    x = dog.getV_vector();
    update_sub_edgeCoords();
    sub_dogsolver[1]->update_edge_coords(sub_edgeCoords[1]);

    Eigen::MatrixXd V0 = sub_dog[0]->getV();
    Eigen::MatrixXd V1 = sub_dog[1]->getV();
    Eigen::MatrixXd target_Vt = (proc_T[0]*V0);
    Eigen::MatrixXd source_Vt = (proc_T[1]*V1);

    Eigen::MatrixXd R;
    Eigen::VectorXd t;
    double scale;
                                      //scale,reflect
    igl::procrustes(source_Vt, target_Vt,true,false, scale,R,t);
    double alpha=p.admm_gamma;
    V1 = ((V1 * scale * R).rowwise() + t.transpose())*alpha + (1-alpha)*V1;
    sub_dog[1]->update_V(V1);
    sub_dogsolver[1]->set_opt_vars(sub_dog[1]->getV_vector());

    double sub_cd = 0.0;
    double sub_obj = 0.0;
    cout << "Solving submesh 1\n";
    sub_dogsolver[1]->set_z(-sub_dogsolver[0]->get_Ax());
    sub_dogsolver[1]->single_iteration_procrustes(sub_cd, sub_obj);
    constraints_deviation += sub_cd;
    objective += sub_obj;

    dog.update_submesh_V(1, sub_dog[1]->getV());
    x = dog.getV_vector();

    update_sub_edgeCoords();
    sub_dogsolver[0]->update_edge_coords(sub_edgeCoords[0]);

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_experimental(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (experimental)" << endl;
	x0 = x;
  if(is_subsolver()){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();

    constraints_deviation = 0.0;
    objective = 0.0;

    //Do iteration on global mesh as "guess"
    double d0=0,d1=0;
    if(p.admm_gamma>0.5){
      single_iteration_normal(d0,d1);
    }else{
      single_iteration_fold(d0,d1);
    }
    --iter_i;

    //update_sub_edgeCoords();//this one uses sub_dog
    {//this one uses dog
      const DogEdgeStitching& edgeStitching = dog.getEdgeStitching();
      vector<int> k(num_submeshes,0);
      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v11 = edgeStitching.edge_const_1[j].v1;
        int v12 = edgeStitching.edge_const_1[j].v2;
        int v21 = edgeStitching.edge_const_2[j].v1;
        int v22 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];

        int submesh_1 = dog.v_to_submesh_idx(v11);
        int submesh_2 = dog.v_to_submesh_idx(v21);
        EdgePoint ep1(Edge(v11, v12), t);
        EdgePoint ep2(Edge(v21, v22), t);
        sub_edgeCoords[submesh_1].row(k[submesh_1]) = ep1.getPositionInMesh(dog.getV());
        sub_edgeCoords[submesh_2].row(k[submesh_2]) = ep1.getPositionInMesh(dog.getV());//currently using the same point for both submeshes
        k[submesh_1]++;
        k[submesh_2]++;
      }
      //for angles
      if(p.admm_gamma>0.5){
        for(int i=0;i<sub_dog.size();++i){
          Eigen::MatrixXd w_coords(sub_edge_angle_ww[i].rows(), 6);
          int v_num = dog.get_v_num();
          for(int j=0; j<sub_edge_angle_ww[i].rows(); ++j){
            int w1 = sub_edge_angle_ww[i](j,0);
            int w2 = sub_edge_angle_ww[i](j,1);
            w_coords.row(j) << x(w1), x(w1+v_num), x(w1+2*v_num), x(w2), x(w2+v_num), x(w2+2*v_num);
          }
          sub_dogsolver[i]->update_w_coords(w_coords);
        }
      }

    }

    update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0, 0, 0});

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";
      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);
      //sub_dogsolver[i]->set_opt_vars(sub_dog[i]->getV_vector());
      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_normal(sub_cd, sub_obj);
      constraints_deviation += sub_cd;
      objective += sub_obj;

      cout << " subsolver " << i << " done\n";
    }
    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
  	   dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    x = dog.getV_vector();

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_normal(double& constraints_deviation, double& objective) {
  cout << "running a single optimization routine (normal)" << endl;
  x0 = x;

  update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
    p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
    p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
    p.paired_boundary_bending_weight, 0, 0, 0});
  newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
  dog.update_V_vector(x.head(3*dog.get_v_num()));

  objective = obj.compObj.obj(x);
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
  if (time_measurements_log) {
    *time_measurements_log << objective << "," << constraints_deviation << endl;
  }
  ++iter_i;
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

void DogSolver::update_point_coords(Eigen::VectorXd& bc){
  constraints.posConst.update_coords(bc);
  if(sub_dogsolver.size() > 0){
    for(int i=0; i<sub_dogsolver.size(); ++i){
      for(int j=0; j<sub_ij_to_bc[i].size(); ++j){
        sub_bc[i][j] = bc[sub_ij_to_bc[i][j]];
      }
      sub_dogsolver[i]->update_point_coords(sub_bc[i]);
    }
  }
}

void DogSolver::update_sub_edgeCoords(){
  vector<Eigen::MatrixXd> old_coords = sub_edgeCoords;
  for(int i=0; i<corresponding_edge_points.size(); ++i){
    for(int j=0; j<corresponding_edge_points[i].size(); ++j){
      int other_submesh = corresponding_edge_points[i][j].first;
      int other_index = corresponding_edge_points[i][j].second;
      sub_edgeCoords[i].row(j) = constrained_edge_points[other_submesh][other_index].getPositionInMesh(sub_dog[other_submesh]->getV());
    }
  }
}

void DogSolver::set_solver_mode(SolverMode mode_new) {
  mode = mode_new;
  if(sub_dogsolver.size() > 0){
    for(int i=0; i<sub_dogsolver.size(); ++i){
      sub_dogsolver[i]->set_solver_mode(mode_new);
    }
  }
}

void DogSolver::build_VSADMMObjective(const Eigen::SparseMatrix<double>& A){
  admm_A = A;
  vsadmmConst.set_pointers(admm_A, admm_z);
  if(vsadmm_obj) delete vsadmm_obj;
  vsadmm_obj = new QuadraticConstraintsSumObjective(vsadmmConst, x);
  //and also the serial one
  if(serial_obj) delete serial_obj;
  serial_obj = new SomeSerialObjective(admm_A, admm_lambda);
}

void DogSolver::build_ProximalObjective(const Eigen::SparseMatrix<double>& P){
  admm_P = P;
  if(pjadmm_obj) delete pjadmm_obj;
  pjadmm_obj = new ProximalObjective(admm_P, x0);
}

void DogSolver::remake_compobj(){
  obj.compObj = CompositeObjective(
    {&obj.bending, &obj.isoObj, &obj.stitchingConstraintsPenalty,
     &obj.pointsPosSoftConstraints, &obj.edgePosSoftConstraints,
     &obj.ptPairSoftConst, &obj.edgeAnglesSoftConstraints,
     &obj.mvTangentCreaseSoftConstraints, &obj.foldingBinormalBiasObj,
     &obj.foldingMVBiasObj,&obj.pairedBndVertBendingObj,
     vsadmm_obj, pjadmm_obj, serial_obj},
    {p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
     p.stitching_weight,p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
     p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
     p.paired_boundary_bending_weight, p.admm_rho, p.admm_rho, 1});
}

void DogSolver::update_obj_weights(const std::vector<double>& weights_i){
  obj.compObj.update_weights(weights_i);
  if(!is_subsolver()){
    for(int i=0; i<sub_dogsolver.size(); ++i){
      sub_dogsolver[i]->update_obj_weights(weights_i);
    }
  }
}

void DogSolver::update_edge_angles(const std::vector<double> cos_angles_i) {
  constraints.edgeAngleConst.set_angles(cos_angles_i);

  if(!is_subsolver()){
    int v_num = dog.get_v_num();

    for(int i=0; i<sub_dogsolver.size(); ++i){
      vector<double> angles(sub_idx_to_angle_idx[i].size());
      Eigen::MatrixXd w_coords(sub_edge_angle_ww[i].rows(), 6);

      for(int j=0; j<sub_idx_to_angle_idx[i].size(); ++j){
        angles[j] = cos_angles_i[ sub_idx_to_angle_idx[i][j] ];
      }

      for(int j=0; j<sub_edge_angle_ww[i].rows(); ++j){
        int w1 = sub_edge_angle_ww[i](j,0);
        int w2 = sub_edge_angle_ww[i](j,1);
        w_coords.row(j) << x(w1), x(w1+v_num), x(w1+2*v_num), x(w2), x(w2+v_num), x(w2+2*v_num);
      }
      sub_dogsolver[i]->update_w_coords(w_coords);

      sub_dogsolver[i]->update_edge_angles(angles);
    }
  } else {
    constraints.subEdgesAngleConst.set_angles(cos_angles_i);
  }
}

void DogSolver::update_w_coords(const Eigen::MatrixXd& W){
  //must be initialized in first iteration
  if(iter_i==0) constraints.subEdgesAngleConst.init_outside_points(W, x);
  else constraints.subEdgesAngleConst.update_outside_points(W);
}

Eigen::VectorXd DogSolver::get_obj_parts(){
  Eigen::VectorXd res(3);
  res(0) = obj.bending.obj(x);
  res(1) = obj.isoObj.obj(x);
  res(2) = obj.compObj.obj(x);
  return res;
}
