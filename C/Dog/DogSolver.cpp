#include "DogSolver.h"

#include <igl/polar_svd.h>
#include <igl/procrustes.h>

using namespace std;

DogSolver::DogSolver(Dog& dog, Dog& coarse_dog, FineCoarseConversion& conversion,
        const Eigen::VectorXd& init_mesh_vars, const Eigen::VectorXd& coarse_x0, DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc,
        std::vector<EdgePoint>& edgePoints, Eigen::MatrixXd& edgeCoords,
        std::vector<std::pair<Edge,Edge>>& edge_angle_pairs, std::vector<double>& edge_cos_angles,
        std::vector<MVTangentCreaseFold>& mvTangentCreaseAngleParams, std::vector<double>& mv_cos_angles,
        std::vector<std::pair<int,int>>& pairs,
        std::vector<pair<int,int>>& bnd_vertices_pairs,
        Eigen::VectorXi& link_b,
        std::ofstream* time_measurements_log, bool ismainsolver) :
          dog(dog), coarse_dog(coarse_dog), fine_coarse(conversion), x(dog.getV_vector()),
          foldingBinormalBiasConstraints(dog),
          foldingMVBiasConstraints(dog, p.flip_sign), p(p),
          link_bc(link_b.size()), linkVerticesConst(link_b, link_bc),
          constraints(dog, init_mesh_vars, b, bc, edgePoints, edgeCoords, edge_angle_pairs, edge_cos_angles, mvTangentCreaseAngleParams,
                      mv_cos_angles, pairs),
          obj(dog, init_mesh_vars, constraints, foldingBinormalBiasConstraints, foldingMVBiasConstraints, bnd_vertices_pairs, linkVerticesConst, p),
          newtonKKT(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p),
          time_measurements_log(time_measurements_log), is_main_solver(ismainsolver)
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

    if(is_main_solver){
      sub_dog.resize(num_submeshes);
      sub_dogsolver.resize(num_submeshes);

      int v_num = dog.get_v_num();

      for(int i=0; i<num_submeshes; ++i) sub_dog[i] = dog.get_submesh(i);
      const DogEdgeStitching& coarse_es = coarse_dog.getEdgeStitching();

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

      for(int i=0; i<edgeStitching.edge_coordinates.size(); ++i){
        int v1 = edgeStitching.edge_const_1[i].v1;
        int v2 = edgeStitching.edge_const_1[i].v2;
        int w1 = edgeStitching.edge_const_2[i].v1;
        int w2 = edgeStitching.edge_const_2[i].v2;
        double t = edgeStitching.edge_coordinates[i];

        int submesh_1 = dog.v_to_submesh_idx(v1);
        int submesh_2 = dog.v_to_submesh_idx(w1);

        corresponding_edge_points[submesh_2].push_back(pair<int,int>(submesh_1, constrained_edge_points[submesh_1].size()));
        corresponding_edge_points[submesh_1].push_back(pair<int,int>(submesh_2, constrained_edge_points[submesh_2].size()));

        constrained_edge_points[submesh_1].push_back(EdgePoint(Edge(dog.v_in_submesh(v1), dog.v_in_submesh(v2)),t));
        constrained_edge_points[submesh_2].push_back(EdgePoint(Edge(dog.v_in_submesh(w1), dog.v_in_submesh(w2)),t));
      }

      // initialize sub_edgeCoords and curve_ep_to_sub_edgeCoords
      for(int i=0; i<num_submeshes; ++i) {
        sub_edgeCoords[i].resize(constrained_edge_points[i].size(), 3);
        sub_edgeCoords[i].setConstant(-42);
      }
      curve_ep_to_sub_edgeCoords.resize(edgeStitching.stitched_curves.size());
      for(int i=0; i<edgeStitching.stitched_curves.size(); ++i){
        curve_ep_to_sub_edgeCoords[i].resize(edgeStitching.stitched_curves[i].size(), 4);
      }
      vector<int> current_row(num_submeshes,0);
      for(int i=0; i<edgeStitching.edge_coordinates.size(); ++i){
        int v1 = edgeStitching.edge_const_1[i].v1;
        int v2 = edgeStitching.edge_const_1[i].v2;
        int w1 = edgeStitching.edge_const_2[i].v1;
        int w2 = edgeStitching.edge_const_2[i].v2;
        double t = edgeStitching.edge_coordinates[i];
        int submesh_1 = dog.v_to_submesh_idx(v1);
        int submesh_2 = dog.v_to_submesh_idx(w1);

        for(int j=0; j<edgeStitching.stitched_curves.size(); ++j){
          for(int k=0; k<edgeStitching.stitched_curves[j].size(); ++k){
            int vc1 = edgeStitching.stitched_curves[j][k].edge.v1;
            int vc2 = edgeStitching.stitched_curves[j][k].edge.v2;
            if(vc1 == v1 && vc2 == v2){
              sub_edgeCoords[submesh_1].row(current_row[submesh_1]) =
                EdgePoint(Edge(w1, w2),t).getPositionInMesh(dog.getV());
              sub_edgeCoords[submesh_2].row(current_row[submesh_2]) =
                EdgePoint(Edge(v1, v2),t).getPositionInMesh(dog.getV());
              curve_ep_to_sub_edgeCoords[j].row(k) << submesh_1, current_row[submesh_1], submesh_2, current_row[submesh_2];
              ++current_row[submesh_1];
              ++current_row[submesh_2];
            }
          }
        }
      }
      is_constrained = (is_constrained || edgeStitching.edge_coordinates.size()>0);

      std::vector< Eigen::SparseMatrix<double> > sub_admm_A(num_submeshes);

      for(int i=0; i<num_submeshes; ++i){
        sub_admm_A[i] = Eigen::SparseMatrix<double>(edgeStitching.edge_coordinates.size(), 3*sub_dog[i]->get_v_num());
        sub_admm_A[i].reserve(6*sub_admm_A[i].rows());
      }

      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v1 = edgeStitching.edge_const_1[j].v1;
        int v2 = edgeStitching.edge_const_1[j].v2;
        int w1 = edgeStitching.edge_const_2[j].v1;
        int w2 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];
        int submesh_1 = dog.v_to_submesh_idx(v1);
        int submesh_2 = dog.v_to_submesh_idx(w1);

        int subsize = dog.get_submesh_i_size(submesh_1);
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v1)) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v1) + subsize) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v1) + 2*subsize) = t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v2)) = 1.0-t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v2) + subsize) = 1.0-t;
        sub_admm_A[submesh_1].insert(j, dog.v_in_submesh(v2) + 2*subsize) = 1.0-t;

        subsize = dog.get_submesh_i_size(submesh_2);
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w1)) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w1) + subsize) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w1) + 2*subsize) = -t;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w2)) = t-1.0;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w2) + subsize) = t-1.0;
        sub_admm_A[submesh_2].insert(j, dog.v_in_submesh(w2) + 2*subsize) = t-1.0;
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

      //link vertex "constraints" (they're more like a suggestion)
      Eigen::VectorXi fine_link_b = fine_coarse.get_fine_link_b();
      Eigen::VectorXi coarse_link_b = fine_coarse.get_coarse_link_b();
      vector<Eigen::VectorXi> sub_link_b(num_submeshes);
      vector<Eigen::VectorXd> sub_link_bc(num_submeshes);
      {
        sub_links_size.resize(num_submeshes, 0);
        int num_links = fine_link_b.size() / 3;
        for(int i=0; i<num_links; ++i){
          link_bc[i] = x[fine_link_b[i]];
          link_bc[i + num_links] = x[fine_link_b[i + num_links]];
          link_bc[i + 2 * num_links] = x[fine_link_b[i + 2 * num_links]];
          ++sub_links_size[dog.v_to_submesh_idx(fine_link_b[i])];
        }
        update_link_vertices_constraints(link_bc);
        for(int i=0; i<num_submeshes; ++i){
          sub_link_b[i].resize(3 * sub_links_size[i]);
          sub_link_bc[i].resize(3 * sub_links_size[i]);
        }
        vector<int> sub_i(num_submeshes, 0);
        for(int i=0; i<num_links; ++i){
          int u = fine_link_b[i];
          int patch = dog.v_to_submesh_idx(u);
          sub_link_b[patch][sub_i[patch]] = dog.v_in_submesh(u);
          sub_link_b[patch][sub_i[patch] + sub_links_size[patch]] = dog.v_in_submesh(u) + sub_dog[patch]->get_v_num();
          sub_link_b[patch][sub_i[patch] + 2 * sub_links_size[patch]] = dog.v_in_submesh(u) + 2 * sub_dog[patch]->get_v_num();
          sub_link_bc[patch][sub_i[patch]] = x[u];
          sub_link_bc[patch][sub_i[patch] + sub_links_size[patch]] = x[u + v_num];
          sub_link_bc[patch][sub_i[patch] + 2 * sub_links_size[patch]] = x[u + 2 * v_num];
          ++sub_i[patch];
        }
      }

      // Procrustes stuff
      proc_T.resize(sub_dog.size());
      vector<int> T_i(sub_dog.size(), 0);

      for(int i=0; i<edgeStitching.stitched_curves.size(); ++i){
        for(int j=0; j<edgeStitching.stitched_curves[i].size(); ++j){
          Edge e1 = edgeStitching.stitched_curves[i][j].edge;
          int c = edgeStitching.edge_to_duplicates.at(e1);
          int k = edgeStitching.multiplied_edges_start[c];
          //the number of duplicates is ignored because it's only >1 in cases that are buggy
          Edge e2 = edgeStitching.edge_const_2[k];
          int patch1 = dog.v_to_submesh_idx(e1.v1);
          int patch2 = dog.v_to_submesh_idx(e2.v1);
          ++T_i[patch1];
          ++T_i[patch2];
        }
      }
      for(int i=0; i<sub_dog.size(); ++i) {
        proc_T[i].resize(T_i[i], sub_dog[i]->get_v_num());
        T_i[i] = 0;
      }
      for(int i=0; i<edgeStitching.stitched_curves.size(); ++i){
        for(int j=0; j<edgeStitching.stitched_curves[i].size(); ++j){
          Edge e1 = edgeStitching.stitched_curves[i][j].edge;
          int c = edgeStitching.edge_to_duplicates.at(e1);
          int k = edgeStitching.multiplied_edges_start[c];
          //the number of duplicates is ignored because it's only >1 in cases that are buggy
          Edge e2 = edgeStitching.edge_const_2[k];
          double t = edgeStitching.edge_coordinates[k];
          int patch1 = dog.v_to_submesh_idx(e1.v1);
          int patch2 = dog.v_to_submesh_idx(e2.v1);
          proc_T[patch1].insert(T_i[patch1], dog.v_in_submesh(e1.v1)) = t;
          proc_T[patch1].insert(T_i[patch1], dog.v_in_submesh(e1.v2)) = 1.0-t;
          proc_T[patch2].insert(T_i[patch2], dog.v_in_submesh(e2.v1)) = t;
          proc_T[patch2].insert(T_i[patch2], dog.v_in_submesh(e2.v2)) = 1.0-t;
          ++T_i[patch1];
          ++T_i[patch2];
        }
      }

      //Construct sub dogsolvers
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
        sub_dogsolver[i] = new DogSolver(*sub_dog[i], empty_dog, empty_conversion,
              sub_x0, empty_xd, p,
              sub_b[i], sub_bc[i],
              //empty_xi, empty_xd,
              constrained_edge_points[i], sub_edgeCoords[i],
              //empty_ep, empty_mat,
              sub_egg[i], sub_angles[i],
              //empty_egg, empty_d,
              empty_thing, empty_d,
              empty_pair,
              empty_pair,
              sub_link_b[i],
              NULL, false);

        sub_dogsolver[i]->update_link_vertices_constraints(sub_link_bc[i]);

        sub_dogsolver[i]->set_lambda(Eigen::VectorXd::Zero(sub_admm_A[i].rows()));
        sub_dogsolver[i]->set_z(Eigen::VectorXd::Zero(sub_admm_A[i].rows()));
        sub_dogsolver[i]->build_VSADMMObjective(sub_admm_A[i]);

        Eigen::SparseMatrix<double> I_i(sub_x0.size(), sub_x0.size());
        I_i.setIdentity(); I_i *= 0.1*(num_submeshes-1) *0.01;
        sub_dogsolver[i]->build_ProximalObjective(I_i);
        //sub_dogsolver[i]->remake_compobj();//incompatible with link vertex constraints

        cout << " constructed subsolver "<<i<<endl;
      }

      cout << "Sub dogsolvers constructed\n";

      // Construct coarse solver
      int coarse_v_num = coarse_dog.get_v_num();
      //position constraints
      vector<int> which_b(0);
      int bthird = b.size()/3;
      coarse_b_to_bi.clear();
      for(int i=0; i<bthird; ++i){
        int coarse_v = fine_coarse.fine_to_coarse(b[i]);
        if(coarse_v > -1) {
          which_b.push_back(coarse_v);
          coarse_b_to_bi.push_back(i);
        }
        else cout << "Warning: Position constraint "<<i<<" is not in coarse: fine "<<b[i]<<"   coarse "<<coarse_v<<"\n";
      }
      coarse_b.resize(which_b.size()*3);
      coarse_bc.resize(coarse_b.size());
      for(int i=0; i<which_b.size(); ++i){
        coarse_b(i) = which_b[i];
        coarse_b(i + which_b.size()) = which_b[i] + coarse_v_num;
        coarse_b(i + 2*which_b.size()) = which_b[i] + 2*coarse_v_num;
        coarse_bc(i) = bc(coarse_b_to_bi[i])*0.5;//coarse scale
        coarse_bc(i + which_b.size()) = bc(coarse_b_to_bi[i] + bthird)*0.5;//coarse scale
        coarse_bc(i + 2*which_b.size()) = bc(coarse_b_to_bi[i] + 2*bthird)*0.5;//coarse scale
      }

      //coarse angle constraints
      const QuadTopology& coarse_qt = coarse_dog.getQuadTopology();
      vector<std::pair<Edge,Edge>> coarse_edge_angle_pairs(0);
      vector<double> coarse_edge_cos_angles(0);
      coarse_angle_idx.clear();
      coarse_angle_vertices.clear();
      coarse_angle_t.clear();
      for(int i=0; i<edge_angle_pairs.size(); ++i){
        int v1 = edge_angle_pairs[i].first.v1;
        int v2 = edge_angle_pairs[i].first.v2;
        int c_v1 = fine_coarse.fine_to_coarse(v1);
        int c_v2 = fine_coarse.fine_to_coarse(v2);
        if(c_v1 < 0 && c_v2 < 0) continue;//One vertex has to be a link point
        int w1 = edge_angle_pairs[i].second.v1;
        int w2 = edge_angle_pairs[i].second.v2;
        int edge_idx_v, edge_idx_w;
        vector<double> c_t(2,0.5);//coarse scale ?
        if(c_v2 < 0){
          edge_idx_v = fine_coarse.fine_to_coarse_edge(v2);
          edge_idx_w = fine_coarse.fine_to_coarse_edge(w2);
          c_t[0] = 1;
        } else {
          edge_idx_v = fine_coarse.fine_to_coarse_edge(v1);
          edge_idx_w = fine_coarse.fine_to_coarse_edge(w1);
          c_t[1] = 0;
        }
        vector<int> c_points(4);
        c_points[0] = coarse_qt.E(edge_idx_v, 0);
        c_points[1] = coarse_qt.E(edge_idx_v, 1);
        c_points[2] = coarse_qt.E(edge_idx_w, 0);
        c_points[3] = coarse_qt.E(edge_idx_w, 1);
        //cout << "Edges are "<<edge_idx_v<<" and "<<edge_idx_w<<"\n";
        Edge ev = Edge(c_points[0], c_points[1]);
        Edge ew = Edge(c_points[2], c_points[3]);
        //cout << "Edge vertices "<<ev.v1<<", "<<ev.v2<<" and "<<ew.v1<<", "<<ew.v2<<"\n";
        coarse_edge_angle_pairs.push_back(pair<Edge,Edge>(ev, ew));
        coarse_edge_cos_angles.push_back(edge_cos_angles[i]);
        coarse_angle_idx.push_back(i);

        coarse_angle_vertices.push_back(c_points);
        coarse_angle_t.push_back(c_t);
        Eigen::RowVector3d cv1 = 2* ( c_t[0]*coarse_dog.getV().row(c_points[0]) +(1-c_t[0])*coarse_dog.getV().row(c_points[1]) );
        Eigen::RowVector3d cv2 = 2* ( c_t[1]*coarse_dog.getV().row(c_points[0]) +(1-c_t[1])*coarse_dog.getV().row(c_points[1]) );
        Eigen::RowVector3d cw1 = 2* ( c_t[0]*coarse_dog.getV().row(c_points[2]) +(1-c_t[0])*coarse_dog.getV().row(c_points[3]) );
        Eigen::RowVector3d cw2 = 2* ( c_t[1]*coarse_dog.getV().row(c_points[2]) +(1-c_t[1])*coarse_dog.getV().row(c_points[3]) );

      }

      cout << "building coarse solver\n";
      coarse_solver = new DogSolver(coarse_dog, empty_dog, empty_conversion,
            coarse_x0, empty_xd, p,
            coarse_b, coarse_bc,
            empty_ep, empty_mat,
            coarse_edge_angle_pairs, coarse_edge_cos_angles,
            empty_thing, empty_d,
            empty_pair,
            empty_pair,
            coarse_link_b,
            NULL, false);
      cout << "Coarse solver constructed\n";
      for(int i=0; i<link_bc.size(); ++i){
        link_bc[i] = coarse_x0[coarse_link_b[i]];
      }
      coarse_solver->update_link_vertices_constraints(link_bc);
      coarse_curves.resize(coarse_es.stitched_curves.size());
      for(int i=0; i<coarse_curves.size(); ++i)
        coarse_curves[i].edgePoints = coarse_es.stitched_curves[i];


      int num_curves = edgeStitching.stitched_curves.size();
      for(int i=0; i<num_curves; ++i){
        int fine_begin = fine_coarse.coarse_to_fine_curve(i, 0);
        int coarse_end = coarse_es.stitched_curves[i].size()-1;
        int fine_end = fine_coarse.coarse_to_fine_curve(i, coarse_end);
        for(int j=0; j<edgeStitching.stitched_curves[i].size(); ++j){
          if(j<fine_begin || j > fine_end){
            int submesh_1 = curve_ep_to_sub_edgeCoords[i](j, 0);
            int k1 = curve_ep_to_sub_edgeCoords[i](j, 1);
            sub_dogsolver[submesh_1]->use_sub_edgepoint_constraint(k1, false);
            int submesh_2 = curve_ep_to_sub_edgeCoords[i](j, 2);
            int k2 = curve_ep_to_sub_edgeCoords[i](j, 3);
            sub_dogsolver[submesh_2]->use_sub_edgepoint_constraint(k2, false);
          }
        }
      }

      cout<<"main solver done\n";

    } else {
      sub_dog.clear();
      sub_dogsolver.clear();
      coarse_solver = nullptr;
    }

  // Anderson Acceleration
  aa_dF.resize(x.size(), aa_m);
  aa_dG.resize(x.size(), aa_m);
  aa_x = x;
  aa_M.resize(aa_m, aa_m); aa_M.setZero();
  aa_scale.resize(aa_m); aa_scale.setZero();
}

DogSolver::~DogSolver(){
  if(vsadmm_obj) delete vsadmm_obj;
  if(sub_dog.size() > 0){
    for(int i=0; i<sub_dog.size(); ++i) delete sub_dog[i];
  }
  if(sub_dogsolver.size() > 0){
    for(int i=0; i<sub_dogsolver.size(); ++i) delete sub_dogsolver[i];
  }
  if(coarse_solver) delete coarse_solver;
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
          PositionalConstraints& linkVerticesConstraints,
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
        linkVerticesObj(linkVerticesConstraints, init_x0),
        //allConstQuadraticObj(constraints, init_x0),
        /*curvedFoldingBiasObj(curvedFoldingBiasObj),*/

        compObj(
          {&bending, &isoObj, &stitchingConstraintsPenalty,
           &pointsPosSoftConstraints, &edgePosSoftConstraints, &ptPairSoftConst,
           &edgeAnglesSoftConstraints, &mvTangentCreaseSoftConstraints,
           &foldingBinormalBiasObj, &foldingMVBiasObj,&pairedBndVertBendingObj,
           &linkVerticesObj},
          {p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
           p.stitching_weight,p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
           p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight,
           p.mv_bias_weight, p.paired_boundary_bending_weight, 0.0})
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
      // Skip degenerate edges where ep is too close to one vertex
      if ( (edge_t < 0.05 ) || ( edge_t > 0.95) ) continue;

      if ((eS.get_vertex_edge_point_deg(edge_pt.edge) == 1) && dog.is_crease_vertex_flat(fold_curve_idx,e_idx) ) continue;
      int v1,v2,w1,w2; dog.get_2_submeshes_vertices_from_edge(edge_pt.edge, v1,v2,w1,w2);
      Eigen::VectorXd V1_p(dog.getV().row(v1)), V2_p(dog.getV().row(v2));
      Eigen::VectorXd W1_p(dog.getV().row(w1)), W2_p(dog.getV().row(w2));
      auto ep_p = edge_pt.getPositionInMesh(dog.getV()); auto ep_b_p = ep_b.getPositionInMesh(dog.getV()); auto ep_f_p = ep_f.getPositionInMesh(dog.getV());
      Eigen::VectorXd B = ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).normalized();

      Eigen::VectorXd e1 = V1_p-V2_p, e2 = W1_p-W2_p;
      //std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
      double sign1 = B.dot(e1), sign2 = B.dot(e2); double flat_tolerance = 1e-12; // ignore flat points..
      // Skip degenerate angles, where curve and edge are almost parallel
      if ( (abs((ep_p-ep_b_p).normalized().dot(e1.normalized())) > 0.9995) ||
           (abs((ep_f_p-ep_p).normalized().dot(e1.normalized())) > 0.9995) ) continue;
      auto bp_dir = (ep_b_p-ep_p).normalized(); auto pf_dir = (ep_p-ep_f_p).normalized();
      double curve_angle = acos( bp_dir.dot(pf_dir) );
      double curvature = sin(curve_angle) / (ep_f_p - ep_b_p).norm();
      // Skip degenerate osculating plane, where both curve segments have too similar direction
      if ( curvature < 1e-5 ) continue;
      //if ( ((ep_p-ep_b_p).normalized().dot(e1.normalized()) < 0.1) || ((ep_p-ep_b_p).normalized().dot(e1.normalized())<0.1) ) {std::cout << "e_idx = " << e_idx << std::endl; continue;}
      if ( sign1*sign2 > 0) {
        is_folded = false;
        //cout << "Change!" << endl;
        //cout << "Edge: "<<v1<<", "<<v2<<" on one submesh, "<<w1<<", "<<w2<<" on the other\n";
        /*
        cout << "(ep_p-ep_b_p).norm() = " << (ep_p-ep_b_p).norm()  << std::endl;
        cout << "(ep_f_p-ep_p).norm() = " << (ep_f_p-ep_p).norm()  << std::endl;
        cout << "(ep_p-ep_b_p).normalized().dot(e1.normalized()) " << (ep_p-ep_b_p).normalized().dot(e1.normalized()) << std::endl;
        cout << "(ep_f_p-ep_p).normalized().dot(e1.normalized()) " << (ep_f_p-ep_p).normalized().dot(e1.normalized()) << std::endl;
        cout << "(ep_p-ep_b_p).normalized().dot(e2.normalized()) " << (ep_p-ep_b_p).normalized().dot(e2.normalized()) << std::endl;
        cout << "(ep_f_p-ep_p).normalized().dot(e2.normalized()) " << (ep_f_p-ep_p).normalized().dot(e2.normalized()) << std::endl;
        std::cout << "((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() = " << ((ep_p-ep_b_p).cross(ep_f_p-ep_p)).norm() << endl;
        */
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
  igl::Timer timer;
  double init_time = timer.getElapsedTime();
  double t = init_time;
  if( dog.get_submesh_n() > 1 && mode != mode_standard) {
    if(mode == mode_subsolvers) {
      single_iteration_subsolvers(constraints_deviation, objective);
    } else if (mode == mode_vsadmm || mode == mode_jadmm || mode == mode_proxjadmm) {
      single_iteration_ADMM(constraints_deviation, objective);
    } else if (mode == mode_serial) {
      single_iteration_serial(constraints_deviation, objective);
    } else if (mode == mode_procrustes) {
      single_iteration_procrustes(constraints_deviation, objective);
    } else if (mode == mode_cheatguess) {
      single_iteration_cheat_guess(constraints_deviation, objective);
    } else if (mode == mode_coarseguess) {
      single_iteration_coarse_guess(constraints_deviation, objective);
    } else if (mode == mode_coarseprocrustes) {
      single_iteration_coarse_procrustes(constraints_deviation, objective);
    } else if (mode == mode_experimental) {
      single_iteration_experimental(constraints_deviation, objective);
    }
  } else {
    if (!p.folding_mode) p.fold_bias_weight = 0;
    if (p.folding_mode) single_iteration_fold(constraints_deviation, objective);
    else single_iteration_normal(constraints_deviation, objective);
  }
  t = timer.getElapsedTime();
  iteration_time = t - init_time;
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
      if (p.fold_bias_weight <= 0.0) {
        cout << "Error: Nonpositive fold bias weight" << endl;
        exit(1);
      }
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
  if(!is_main_solver){
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
  cout << "WARNING: This method is no longer supported (honestly, it wasn't good in the first place)\n";
	x0 = x;
  if(!is_main_solver){
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
  if(!is_main_solver){
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
  if(!is_main_solver){
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

void DogSolver::single_iteration_cheat_guess(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (global cheat iteration)" << endl;
	x0 = x;
  if(!is_main_solver){
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
    single_iteration_normal(d0,d1);
    --iter_i;

    //update_sub_edgeCoords();//this one uses sub_dog
    {//this one uses dog
      const DogEdgeStitching& edgeStitching = dog.getEdgeStitching();
      vector<int> k(num_submeshes,0);
      for(int j=0; j<edgeStitching.edge_coordinates.size(); ++j){
        int v1 = edgeStitching.edge_const_1[j].v1;
        int v2 = edgeStitching.edge_const_1[j].v2;
        int w1 = edgeStitching.edge_const_2[j].v1;
        int w2 = edgeStitching.edge_const_2[j].v2;
        double t = edgeStitching.edge_coordinates[j];

        int submesh_1 = dog.v_to_submesh_idx(v1);
        int submesh_2 = dog.v_to_submesh_idx(w1);
        EdgePoint ep1(Edge(v1, v2), t);
        EdgePoint ep2(Edge(w1, w2), t);
        Eigen::RowVector3d target = 0.5*(ep1.getPositionInMesh(dog.getV()) + ep2.getPositionInMesh(dog.getV()));
        sub_edgeCoords[submesh_1].row(k[submesh_1]) = target;
        sub_edgeCoords[submesh_2].row(k[submesh_2]) = target;//currently using the same point for both submeshes
        k[submesh_1]++;
        k[submesh_2]++;
      }
      //for angles
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
      sub_dogsolver[i]->single_iteration_cheat_guess(sub_cd, sub_obj);
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

void DogSolver::single_iteration_coarse_guess(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (coarse guess)" << endl;
	x0 = x;
  if(!is_main_solver){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();

    constraints_deviation = 0.0;
    objective = 0.0;

    cout << "coarse solver does global iteration\n";
    update_obj_weights({p.bending_weight, p.isometry_weight/coarse_dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, p.ftc_weight});
    double d0,d1;
    if(p.folding_mode) coarse_solver->single_iteration_fold(d0,d1);
    else coarse_solver->single_iteration_normal(d0,d1);

    // update stitching constraint coordinates
    coarse_to_fine_update();

    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      0, 0, p.fold_bias_weight, p.mv_bias_weight,// no angle constraints on subsolvers!
      p.paired_boundary_bending_weight, p.ftc_weight});

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";
      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);
      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_coarse_guess(sub_cd, sub_obj);
      constraints_deviation += sub_cd;
      objective += sub_obj;

      cout << " subsolver " << i << " done\n";
    }
    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
  	   dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    x = dog.getV_vector();

    //update coarse mesh?
    if(p.ftc_update) fine_to_coarse_update();

    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0});

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_coarse_procrustes(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (coarse + procrustes)" << endl;
	x0 = x;
  if(!is_main_solver){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));
    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();
    constraints_deviation = 0.0;
    objective = 0.0;

    update_obj_weights({p.bending_weight, p.isometry_weight/coarse_dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0, 0, 0});

    cout << "coarse solver does global iteration\n";
    double d0,d1;
    if(p.folding_mode) coarse_solver->single_iteration_fold(d0,d1);
    else coarse_solver->single_iteration_normal(d0,d1);

    // update stitching constraint coordinates
    vector<Eigen::MatrixXd> target_coords(num_submeshes);
    vector<int> row_i(num_submeshes, 0);
    for(int i=0; i<num_submeshes; ++i) target_coords[i].resize(proc_T[i].rows(), 3);
    const Eigen::MatrixXd& coarse_V = coarse_solver->getDog().getV();
    for(int i=0; i<coarse_curves.size(); ++i){
      Eigen::MatrixXd fine_coords = fine_coarse.getInterpolatedCurveCoords(dog, coarse_dog, i) /0.5;//coarse scale

      for(int j=0; j<fine_coords.rows(); ++j){
        int submesh_1 = curve_ep_to_sub_edgeCoords[i](j,0);
        int k1 = curve_ep_to_sub_edgeCoords[i](j,1);
        int submesh_2 = curve_ep_to_sub_edgeCoords[i](j,2);
        int k2 = curve_ep_to_sub_edgeCoords[i](j,3);

        target_coords[submesh_1].row(row_i[submesh_1]++) =
         sub_edgeCoords[submesh_1].row(k1) = fine_coords.row(j);
        target_coords[submesh_2].row(row_i[submesh_2]++) =
         sub_edgeCoords[submesh_2].row(k2) = fine_coords.row(j);
      }
    }

    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      0, 0, p.fold_bias_weight, p.mv_bias_weight,// no angle constraints on subsolvers!
      p.paired_boundary_bending_weight, 0, 0, 0});

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";

      Eigen::MatrixXd target_Vt = target_coords[i];

      Eigen::MatrixXd Vi = sub_dog[i]->getV();
      Eigen::MatrixXd source_Vt = proc_T[i] * Vi;

      Eigen::MatrixXd R;
      Eigen::VectorXd t;
      double scale;
      //                                   scale, reflect
      igl::procrustes(source_Vt, target_Vt, false,false, scale,R,t);
      //double alpha=p.admm_gamma;
      Vi = (Vi /* scale*/ * R).rowwise() + t.transpose();
      sub_dog[i]->update_V(Vi);
      sub_dogsolver[i]->set_opt_vars(sub_dog[i]->getV_vector());

      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);
      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_coarse_procrustes(sub_cd, sub_obj);
      constraints_deviation += sub_cd;
      objective += sub_obj;

      cout << " subsolver " << i << " done\n";
    }
    //propagate local solutions to global dog
    for(int i=0; i<num_submeshes; ++i){
      dog.update_submesh_V(i, sub_dog[i]->getV());
    }
    x = dog.getV_vector();

    if(p.ftc_update) fine_to_coarse_update();

    cout << "all subsolvers done\n";
  }
  ++iter_i;
}

void DogSolver::single_iteration_experimental(double& constraints_deviation, double& objective) {
	cout << "running a single optimization routine (experimental)" << endl;
	x0 = x;
  if(!is_main_solver){
    newtonKKT.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);
    dog.update_V_vector(x.head(3*dog.get_v_num()));

    constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
    objective = obj.compObj.obj(x);
  } else {
    cout << " with subsolvers\n";
    int num_submeshes = dog.get_submesh_n();
    constraints_deviation = 0.0;
    objective = 0.0;

    update_obj_weights({p.bending_weight, p.isometry_weight/coarse_dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0, 0, 0});

    cout << "coarse solver does global iteration\n";
    double d0,d1;
    if(p.folding_mode) coarse_solver->single_iteration_fold(d0,d1);
    else coarse_solver->single_iteration_normal(d0,d1);
    coarse_solver->anderson_acceleration();

    // update stitching constraint coordinates
    vector<Eigen::MatrixXd> target_coords(num_submeshes);
    vector<int> row_i(num_submeshes, 0);
    for(int i=0; i<num_submeshes; ++i) target_coords[i].resize(proc_T[i].rows(), 3);
    const Eigen::MatrixXd& coarse_V = coarse_solver->getDog().getV();
    for(int i=0; i<coarse_curves.size(); ++i){
      Eigen::MatrixXd fine_coords = fine_coarse.getInterpolatedCurveCoords(dog, coarse_dog, i) /0.5;//coarse scale

      for(int j=0; j<fine_coords.rows(); ++j){
        int submesh_1 = curve_ep_to_sub_edgeCoords[i](j,0);
        int k1 = curve_ep_to_sub_edgeCoords[i](j,1);
        int submesh_2 = curve_ep_to_sub_edgeCoords[i](j,2);
        int k2 = curve_ep_to_sub_edgeCoords[i](j,3);

        target_coords[submesh_1].row(row_i[submesh_1]++) =
         sub_edgeCoords[submesh_1].row(k1) = fine_coords.row(j);
        target_coords[submesh_2].row(row_i[submesh_2]++) =
         sub_edgeCoords[submesh_2].row(k2) = fine_coords.row(j);
      }
    }

    update_obj_weights({p.bending_weight, p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, 0.5*p.stitching_weight, p.pair_weight,
      0, 0, p.fold_bias_weight, p.mv_bias_weight,// no angle constraints on subsolvers!
      p.paired_boundary_bending_weight, 0, 0, 0});

    //solve on submeshes
    for(int i=0; i<num_submeshes; ++i){
      cout << "subsolver " << i << " does single iteration\n";

      Eigen::MatrixXd target_Vt = target_coords[i];

      Eigen::MatrixXd Vi = sub_dog[i]->getV();
      Eigen::MatrixXd source_Vt = proc_T[i] * Vi;

      Eigen::MatrixXd R;
      Eigen::VectorXd t;
      double scale;
      //                                   scale, reflect
      igl::procrustes(source_Vt, target_Vt, false,false, scale,R,t);
      //double alpha=p.admm_gamma;
      Vi = (Vi /* scale*/ * R).rowwise() + t.transpose();
      sub_dog[i]->update_V(Vi);
      sub_dogsolver[i]->set_opt_vars(sub_dog[i]->getV_vector());

      sub_dogsolver[i]->update_edge_coords(sub_edgeCoords[i]);
      double sub_cd = constraints_deviation;
      double sub_obj = objective;
      sub_dogsolver[i]->single_iteration_experimental(sub_cd, sub_obj);
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
  if(is_main_solver){
    update_obj_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
      p.stitching_weight, p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
      p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
      p.paired_boundary_bending_weight, 0, 0, 0});
  }
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
  if(coarse_solver){
    for(int i=0; i<coarse_bc.size()/3; ++i){
      //TODO probably wrong (see in constructor)
      coarse_bc(i) = bc(coarse_b_to_bi[i])*0.5;//coarse scale
      coarse_bc(i + coarse_bc.size()/3) = bc(coarse_b_to_bi[i] + bc.size()/3)*0.5;//coarse scale
      coarse_bc(i + 2*(coarse_bc.size()/3)) = bc(coarse_b_to_bi[i] + 2*(bc.size()/3))*0.5;//coarse scale
    }
    coarse_solver->update_point_coords(coarse_bc);
  }
}

void DogSolver::update_sub_edgeCoords(){
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
     &obj.linkVerticesObj, vsadmm_obj, pjadmm_obj, serial_obj},
    {p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(),
     p.stitching_weight,p.soft_pos_weight, p.soft_pos_weight, p.pair_weight,
     p.dihedral_weight, p.dihedral_weight, p.fold_bias_weight, p.mv_bias_weight,
     p.paired_boundary_bending_weight, 0, p.admm_rho, p.admm_rho, 1});
}

void DogSolver::update_obj_weights(const std::vector<double>& weights_i){
  obj.compObj.update_weights(weights_i);
  if(is_main_solver){
    for(int i=0; i<sub_dogsolver.size(); ++i){
      sub_dogsolver[i]->update_obj_weights(weights_i);
    }
    coarse_solver->update_obj_weights(weights_i);
  }
}

void DogSolver::update_edge_angles(const std::vector<double> cos_angles_i) {
  if(cos_angles_i.size() == 0) return;
  constraints.edgeAngleConst.set_angles(cos_angles_i);

  if(is_main_solver){
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
    vector<double> angles(coarse_angle_idx.size());
    for(int i=0; i<coarse_angle_idx.size(); ++i){
      angles[i] = cos_angles_i[ coarse_angle_idx[i] ];
    }
    coarse_solver->update_edge_angles(angles);
  } else {
    //TODO not for coarse_solver
    constraints.subEdgesAngleConst.set_angles(cos_angles_i);
  }
}

void DogSolver::update_w_coords(const Eigen::MatrixXd& W){
  //must be initialized in first iteration
  if(iter_i==0) constraints.subEdgesAngleConst.init_outside_points(W, x);
  else constraints.subEdgesAngleConst.update_outside_points(W);
}

void DogSolver::coarse_to_fine_update(){
  for(int i=0; i<coarse_curves.size(); ++i){
    Eigen::MatrixXd fine_coords = fine_coarse.getInterpolatedCurveCoords(dog, coarse_dog, i) /0.5;//coarse scale
    int coarse_start = fine_coarse.coarse_to_fine_curve(i,0);
    for(int j=0; j<fine_coords.rows(); ++j){
      int submesh_1 = curve_ep_to_sub_edgeCoords[i](j + coarse_start, 0);
      int k1 = curve_ep_to_sub_edgeCoords[i](j + coarse_start, 1);
      int submesh_2 = curve_ep_to_sub_edgeCoords[i](j + coarse_start, 2);
      int k2 = curve_ep_to_sub_edgeCoords[i](j + coarse_start, 3);

      sub_edgeCoords[submesh_1].row(k1) = fine_coords.row(j);
      sub_edgeCoords[submesh_2].row(k2) = fine_coords.row(j);
    }
  }
  if(p.ftc_weight > 0){
    Eigen::VectorXi fine_links = fine_coarse.get_fine_link_b();
    int num_links = fine_links.size() / 3;
    for(int i=0; i<num_links; ++i){
      auto row = coarse_dog.getV().row(fine_coarse.fine_to_coarse(fine_links[i])) *2;//coarse scale
      link_bc[i] = row[0];
      link_bc[i + num_links] = row[1];
      link_bc[i + 2 * num_links] = row[2];
    }
    int v_num = dog.get_v_num();
    int num_submeshes = sub_dog.size();
    vector<Eigen::VectorXd> sub_link_bc(num_submeshes);
    for(int i=0; i<num_submeshes; ++i){
      sub_link_bc[i].resize(3 * sub_links_size[i]);
    }
    vector<int> sub_i(num_submeshes, 0);
    for(int i=0; i<num_links; ++i){
      int u = fine_links[i];
      int patch = dog.v_to_submesh_idx(u);
      sub_link_bc[patch][sub_i[patch]] = link_bc[i];
      sub_link_bc[patch][sub_i[patch] + sub_links_size[patch]] = link_bc[i + num_links];
      sub_link_bc[patch][sub_i[patch] + 2 * sub_links_size[patch]] = link_bc[i + 2 * num_links];
      ++sub_i[patch];
    }
    update_link_vertices_constraints(link_bc);
    for(int i=0; i<num_submeshes; ++i){
      sub_dogsolver[i]->update_link_vertices_constraints(sub_link_bc[i]);
    }
  }
}

void DogSolver::fine_to_coarse_update(){
  auto links = fine_coarse.get_coarse_link_b();
  int num_links = links.size()/3;
  for(int i=0; i<num_links; ++i){
    auto row = dog.getV().row(fine_coarse.coarse_to_fine(links[i])) *0.5;//coarse scale
    link_bc[i] = row[0];
    link_bc[i + num_links] = row[1];
    link_bc[i + 2 * num_links] = row[2];
  }
  coarse_solver->update_link_vertices_constraints(link_bc);
}

double DogSolver::get_obj_val() const {
  return obj.compObj.obj(x);
}
double DogSolver::get_bending_obj_val() const {
  return obj.bending.obj(x);
}
double DogSolver::get_isometry_obj_val() const {
  return obj.isoObj.obj(x);
}
double DogSolver::get_pos_obj_val() const {
  return obj.pointsPosSoftConstraints.obj(x) + obj.edgePosSoftConstraints.obj(x);
}
double DogSolver::get_stitching_obj_val() const {
  return obj.stitchingConstraintsPenalty.obj(x);
}

double DogSolver::get_last_iteration_time() const {
  return iteration_time;
}

void DogSolver::anderson_acceleration(){
  //and now comes Anderson Acceleration
  Eigen::VectorXd aa_G = x;
  Eigen::VectorXd aa_F = aa_G - aa_x;
  if(iter_i == 1){
    // iteration 0 (or 1, because it's increased earlier apparently)
    aa_dF.col(0) = -aa_F;
    aa_dG.col(0) = -aa_G;
    aa_x = aa_G;
  } else {
    // later iterations
    int col_idx = (iter_i+ aa_m-2) % aa_m;
    cout << "aa iteration "<<iter_i <<" uses col "<<col_idx<<"\n";
    aa_dF.col(col_idx) += aa_F;
    aa_dG.col(col_idx) += aa_G;
    double eps = 1e-10;

    double scale = aa_dF.col(col_idx).norm();
    if(scale < eps){cout << "aa scale is "<<scale<<"\tyeah, that's gonna be a problem\n";scale=eps;}
    aa_dF.col(col_idx) /= scale;
    aa_scale(col_idx) = scale;

    int m = std::min(aa_m, iter_i-1);

    Eigen::VectorXd thetas(m); thetas.setZero();
    if (m == 1) {
      //iteration 1
      aa_M(0, 0) = aa_dF.col(0).squaredNorm();//=1 it seems
      //cout << "aa_M 00: "<<aa_M(0,0)<<endl;
      //if(aa_scale(col_idx)>eps) thetas(0) = aa_dF.col(col_idx).dot(aa_F);
      //cout << "aa dfcol 0: "<<aa_dF.col(0).norm()<<endl;
      //cout << "aa F: "<<aa_F.norm()<<endl;
      thetas(0) = aa_dF.col(0).dot(aa_F);
      //cout << "aa thetas 0: "<<thetas(0)<<"\n";
    } else {
      //normal case
      Eigen::VectorXd new_col = (aa_dF.col(col_idx).transpose() * aa_dF.block(0, 0, x.size(), m)).transpose();
      cout << " aa newcol "<<new_col.transpose()<<"\n";
      aa_M.block(col_idx, 0, 1, m) = new_col.transpose();
      aa_M.block(0, col_idx, m, 1) = new_col;

      //Eigen::MatrixXd A = aa_M.block(0, 0, m, m);
      //Eigen::MatrixXd b = aa_dF.block(0, 0, x.size(), m).transpose() * aa_F;

      // Solve normal equation
      //Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;//TODO why is this not working?
      //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> ls_solver;
      //ls_solver.compute(A);//.compute(M)?
      //thetas.head(m) = ls_solver.solve(b);
      //thetas = aa_M.block(0, 0, m, m).colPivHouseholderQr().solve(aa_dF.block(0, 0, x.size(), m).transpose() * aa_F);
      //cout << "aa_dF size "<<aa_dF.rows()<<" x "<<aa_dF.cols()<<"\n";
      //cout << "aa_F size "<<aa_F.size()<<" and btw m="<<m<<"\n";
      //cout << "aa A size "<<A.rows()<<" x "<<A.cols()<<"\n";
      //Eigen::MatrixXd A = aa_dF.block(0,0,aa_dF.rows(),m);

      //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> ls_solver;
      //ls_solver.compute(A.transpose()*A);
      //thetas = ls_solver.solve(A.transpose()*aa_F);
      //thetas = (A.transpose()*A).ldlt().solve(A.transpose()*aa_F);//best result yet
      thetas = (aa_M.block(0, 0, m, m)*aa_M.block(0, 0, m, m)).colPivHouseholderQr().solve(aa_M.block(0, 0, m, m)*aa_F);

      //thetas = (A.transpose()*A).ldlt().solve(b);
    }
    // and a result
    cout << "aa thetas "<<thetas.transpose()<<"\n";

    aa_x = aa_G - aa_dG.block(0, 0, x.size(), m)
                * ((thetas.head(m).array() / aa_scale.head(m).array())
                    .matrix());

    //aa_x = aa_G - thetas*aa_dG.block(0,0,aa_dG.rows(),m);

    int next_col_idx = (col_idx+1) % aa_m;
    aa_dF.col(next_col_idx) = -aa_F;
    aa_dG.col(next_col_idx) = -aa_G;
    //cout << "aa_G "<<aa_G.size()<<"   "<<aa_G.norm()<<"\n";
    //cout << "aa_dG "<<aa_dG.rows()<<" x "<<aa_dG.cols()<<"\n";
    //cout << "thetas "<<thetas<<"\n";
    //cout << "aa_scale "<<aa_scale.transpose()<<"\n";
    //cout << "aa_x "<<aa_x.head(7).transpose()<<"\n";
  }
  double glob_obj = obj.compObj.obj(x);
  double aa_obj = obj.compObj.obj(aa_x);
  cout << "Obj comparison "<<glob_obj<<" <-global  aa-> "<<aa_obj<<endl;
  if(aa_obj < glob_obj){
    cout << "  yaaay\n";
    x = aa_x;
    dog.update_V_vector(x);
  } else {
    aa_x = x;
  }
}

void DogSolver::update_sub_w_coords_from_coarse(){
  if(!is_main_solver) return;
  vector<Eigen::MatrixXd> w_coords(sub_dogsolver.size());
  int cnt = 0;
  for(int i=0; i<sub_dogsolver.size(); ++i){
    //resize assumes all angle constraints are coarse
    w_coords[i].resize(sub_edge_angle_ww[i].rows(), 6);
    cnt += sub_edge_angle_ww[i].rows();
  }
  if(cnt == 0) return;
  vector<int> sub_i(sub_dogsolver.size(), 0);
  for(int i=0; i<coarse_angle_vertices.size(); ++i){
    int v1 = coarse_angle_vertices[i][0];
    int v2 = coarse_angle_vertices[i][1];
    int w1 = coarse_angle_vertices[i][2];
    int w2 = coarse_angle_vertices[i][3];
    double t1 = coarse_angle_t[i][0];
    double t2 = coarse_angle_t[i][1];
    int submesh_v = coarse_dog.v_to_submesh_idx(v1);
    int submesh_w = coarse_dog.v_to_submesh_idx(w1);
    Eigen::RowVector3d v1_vec = coarse_dog.getV().row(v1);
    Eigen::RowVector3d v2_vec = coarse_dog.getV().row(v2);
    Eigen::RowVector3d w1_vec = coarse_dog.getV().row(w1);
    Eigen::RowVector3d w2_vec = coarse_dog.getV().row(w2);
    Eigen::RowVector3d va = t1*v1_vec + (1-t1)*v2_vec;
    Eigen::RowVector3d vb = t2*v1_vec + (1-t2)*v2_vec;
    Eigen::RowVector3d wa = t1*w1_vec + (1-t1)*w2_vec;
    Eigen::RowVector3d wb = t2*w1_vec + (1-t2)*w2_vec;
    //Eigen::RowVector3d dum;dum.setConstant(1);
    w_coords[submesh_v].block(sub_i[submesh_v],0, 1,3) = va *2;//coarse scale
    w_coords[submesh_v].block(sub_i[submesh_v],3, 1,3) = vb *2;//coarse scale
    w_coords[submesh_w].block(sub_i[submesh_w],0, 1,3) = wa *2;//coarse scale
    w_coords[submesh_w].block(sub_i[submesh_w],3, 1,3) = wb *2;//coarse scale
    sub_i[submesh_v] += 1; sub_i[submesh_w] +=1;
    //cout << "wcoords v  "<<w_coords[submesh_v].row(0)<<endl;
    //cout << "wcoords w  "<<w_coords[submesh_w].row(0)<<endl;
  }
  for(int i=0; i<sub_dogsolver.size(); ++i){
    cout << "wcoords "<<i<<"   "<<w_coords[i].row(0)<<endl;
    sub_dogsolver[i]->update_w_coords(w_coords[i]);
  }
}
