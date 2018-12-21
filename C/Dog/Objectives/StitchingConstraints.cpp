#include "StitchingConstraints.h"

using namespace std;

StitchingConstraints::StitchingConstraints(const QuadTopology& quadTop,const DogEdgeStitching& edgeStitching) : quadTop(quadTop), 
																					eS(edgeStitching) {
	const_n= 3*edgeStitching.edge_coordinates.size();
	IJV.resize(2*const_n);
}
Eigen::VectorXd StitchingConstraints::Vals(const Eigen::VectorXd& x) const {
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

		constVals(const_cnt) = 0.5*(e1_c*x(e1.v1+0)+e2_c*x(e1.v2+0)-e1_c*x(e2.v1+0)-e2_c*x(e2.v2+0));
		constVals(const_cnt+1) = 0.5*(e1_c*x(e1.v1+1*vnum)+e2_c*x(e1.v2+1*vnum)-e1_c*x(e2.v1+1*vnum)-e2_c*x(e2.v2+1*vnum));
		constVals(const_cnt+2) = 0.5*(e1_c*x(e1.v1+2*vnum)+e2_c*x(e1.v2+2*vnum)-e1_c*x(e2.v1+2*vnum)-e2_c*x(e2.v2+2*vnum));

		const_cnt+=3;
  }
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  return constVals;
  
}


void StitchingConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {

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