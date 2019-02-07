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
		// Go through all the inner vertices in a curve
		for (int v_i = 1; v_i < stitched_curves[curve_i].size()-1; v_i++) {
			// get the 3 points
			// get the binormal
			// check the dot product
			// if closer to -1,0,0 than to 1,0,0 mark as "yes" by flip_binormal[curve_i].push_back(true/false)
		}
	}
}

Eigen::VectorXd MountainValleyConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
  	#pragma clang loop vectorize(enable)
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