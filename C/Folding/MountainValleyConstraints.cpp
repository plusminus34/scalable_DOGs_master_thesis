#include "MountainValleyConstraints.h"

using namespace std;

MountainValleyConstraints::MountainValleyConstraints(const Eigen::MatrixXd& initV,const DogEdgeStitching& edgeStitching, 
																				  std::vector<bool> is_mountain) :
																					eS(edgeStitching), is_mountain(is_mountain) {
    // TODO: different handling of vertices (maybe no consatraints there?)
	const_n= 2*edgeStitching.edge_coordinates.size();
	IJV.resize(15*const_n);

	// initialize flip_binormal
	flip_binormal.resize(stitched_curves.size());
	for (int curve_i = 0; curve_i < stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];	
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < stitched_curves[curve_i].size()-1; edge_idx++) {
			// get the 3 points
			EdgePoint ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_f = foldingCurve[edge_idx+1];
			auto ep_p = ep.getPositionInMesh(initV); auto ep_b_p = ep_b.getPositionInMesh(initV); auto ep_f_p = ep_f.getPositionInMesh(initV);

			Eigen::RowVector3d b = (ep_f_p-ep_p).cross(ep_p-ep_b_p).normalized();
			Eigen::RowVector3d z_axis = Eigen::RowVector3d(0,0,1);
			// if closer to -1,0,0 than to 1,0,0 mark flip_binormals as "yes"
			// necessary to define M/V assignment correctly
			flip_binormal[curve_i].push_back( (b.dot(z_axis) < 0 ) );
		}
	}
}

Eigen::VectorXd MountainValleyConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
  	
  	for (int curve_i = 0; curve_i < stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < stitched_curves[curve_i].size()-1; edge_idx++) {
			EdgePoint ep = foldingCurve[edge_idx]; ep_b = foldingCurve[edge_idx-1]; ep_f = foldingCurve[edge_idx+1];
			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			bool flip_binormal_sign = (is_mountain[curve_i] ^ flip_binormal[curve_i][edge_idx-1]);

			int ep_v1_i(ep.edge.v1), ep_v2_i(ep.edge.v2); const double ep_t(ep.t);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double ep_v1_x(x(ep_v1_i)); const double ep_v1_y(x(ep_v1_i+1*vnum)); const double ep_v1_z(x(ep_v1_i+2*vnum));
			const double ep_v2_x(x(ep_v2_i)); const double ep_v2_y(x(ep_v2_i+1*vnum)); const double ep_v2_z(x(ep_v2_i+2*vnum));
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			for (int const_side_i = 0; const_side_i < 2 ; const_side_i++) {
				int fold_v_i = fold_v_indices[const_side_i];
				const double fold_v_x(x(fold_v_i)); const double fold_v_y(x(fold_v_i+1*vnum)); const double fold_v_z(x(fold_v_i+2*vnum));
			}
		}
	}


	//for (auto pair : v_equality_matchings) {
	for (int i = 0; i < eS.edge_const_1.size(); i++) {
		Edge e1 = eS.edge_const_1[i]; Edge e2 = eS.edge_const_2[i];
		double e1_c = eS.edge_coordinates[i]; double e2_c = 1-e1_c;

		constVals(const_cnt) = e1_c*x(e1.v1+0)+e2_c*x(e1.v2+0)-e1_c*x(e2.v1+0)-e2_c*x(e2.v2+0);
		constVals(const_cnt+1) = e1_c*x(e1.v1+1*vnum)+e2_c*x(e1.v2+1*vnum)-e1_c*x(e2.v1+1*vnum)-e2_c*x(e2.v2+1*vnum);
		constVals(const_cnt+2) = e1_c*x(e1.v1+2*vnum)+e2_c*x(e1.v2+2*vnum)-e1_c*x(e2.v1+2*vnum)-e2_c*x(e2.v2+2*vnum);

		const_cnt+=3;
  }
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  return constVals;
  
}


void MountainValleyConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {

	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0; int ijv_cnt = 0;
  	#pragma clang loop vectorize(enable)
	for (int i = 0; i < eS.edge_const_1.size(); i++) {
		Edge e1 = eS.edge_const_1[i]; Edge e2 = eS.edge_const_2[i];
		double e1_c = eS.edge_coordinates[i]; double e2_c = 1.-e1_c;

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,e1.v1,e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,e1.v2,e2_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,e2.v1,-e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,e2.v2,-e2_c);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+1,e1.v1+1*vnum,e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+1,e1.v2+1*vnum,e2_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+1,e2.v1+1*vnum,-e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+1,e2.v2+1*vnum,-e2_c);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+2,e1.v1+2*vnum,e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+2,e1.v2+2*vnum,e2_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+2,e2.v1+2*vnum,-e1_c);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt+2,e2.v2+2*vnum,-e2_c);

		const_cnt+=3;
  }
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}