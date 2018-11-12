#include "StitchingConstraints.h"

using namespace std;
/*
struct DogFoldingConstraints {
	std::vector<std::pair<int,int>> edge_const_1, edge_const_2;
	std::vector<double> edge_coordinates;
	// Use for cases when it's important to have a precise representation (usually it doesn't)
	std::vector<CGAL::Exact_predicates_exact_constructions_kernel::FT> edge_coordinates_precise;
};
*/
Eigen::VectorXd StitchingConstraints::Vals(const Eigen::VectorXd& x) const {
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


std::vector<Eigen::Triplet<double> > StitchingConstraints::JacobianIJV(const Eigen::VectorXd& x) const {
	std::vector<Eigen::Triplet<double> > IJV;
	IJV.reserve(4*const_n); // 2 vertices for equality constraints

	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
  	#pragma clang loop vectorize(enable)
	//for (auto pair : v_equality_matchings) {
	for (int i = 0; i < eS.edge_const_1.size(); i++) {
		Edge e1 = eS.edge_const_1[i]; Edge e2 = eS.edge_const_2[i];
		double e1_c = eS.edge_coordinates[i]; double e2_c = 1.-e1_c;

		IJV.push_back(Eigen::Triplet<double>(const_cnt,e1.v1,e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt,e1.v2,e2_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt,e2.v1,-e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt,e2.v2,-e2_c));

		IJV.push_back(Eigen::Triplet<double>(const_cnt+1,e1.v1+1*vnum,e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+1,e1.v2+1*vnum,e2_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+1,e2.v1+1*vnum,-e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+1,e2.v2+1*vnum,-e2_c));

		IJV.push_back(Eigen::Triplet<double>(const_cnt+2,e1.v1+2*vnum,e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+2,e1.v2+2*vnum,e2_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+2,e2.v1+2*vnum,-e1_c));
		IJV.push_back(Eigen::Triplet<double>(const_cnt+2,e2.v2+2*vnum,-e2_c));

		const_cnt+=3;
  }
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
  return IJV;
}