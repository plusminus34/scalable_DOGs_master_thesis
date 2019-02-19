#include "FoldingBinormalBiasConstraints.h"

using namespace std;

FoldingBinormalBiasConstraints::FoldingBinormalBiasConstraints(const Dog& dog, std::vector<bool> is_mountain) :
																					dog(dog), eS(dog.getEdgeStitching()) {
    // TODO: different handling of vertices (maybe no consatraints there?)
	//const_n= 2*(eS.edge_coordinates.size()-2);
	const_n = 2;//*(eS.edge_coordinates.size()-2);
	IJV.resize(21*const_n);
	lambda_hessian_IJV.resize(228*const_n);
}

Eigen::VectorXd FoldingBinormalBiasConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
  	
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		//for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			
			int v1,v2,w1,w2;
			dog->get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*wnum)); const double w1_z(x(w1_i+2*wnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*wnum)); const double w2_z(x(w2_i+2*wnum));

			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double t2 = ep_0_t-1.0;
			double t3 = ep_b_t-1.0;
			double t4 = t2*v2_y;
			double t5 = ep_f_t-1.0;
			double t6 = t2*v2_x;
			double t7 = ep_b_t*ep_b_v1_x;
			double t11 = ep_0_t*v1_x;
			double t20 = ep_b_v2_x*t3;
			double t8 = t6+t7-t11-t20;
			double t9 = t2*v2_z;
			double t10 = ep_f_t*ep_f_v1_x;
			double t12 = ep_b_t*ep_b_v1_y;
			double t19 = ep_0_t*v1_y;
			double t22 = ep_b_v2_y*t3;
			double t13 = t4+t12-t19-t22;
			double t14 = ep_f_t*ep_f_v1_z;
			double t17 = ep_0_t*v1_z;
			double t24 = ep_f_v2_z*t5;
			double t15 = t9+t14-t17-t24;
			double t16 = ep_b_t*ep_b_v1_z;
			double t18 = ep_f_t*ep_f_v1_y;
			double t30 = ep_f_v2_y*t5;
			double t21 = t4+t18-t19-t30;
			double t27 = ep_f_v2_x*t5;
			double t23 = t6+t10-t11-t27;
			double t25 = t8*t15;
			double t29 = ep_b_v2_z*t3;
			double t26 = t9+t16-t17-t29;
			double t28 = t13*t15;
			double t31 = t28-t21*t26;
			constVals(const_cnt++) = (t25-t23*t26)*(w1_y-w2_y)+(t25-t23*(t9+t16-ep_b_v2_z*t3-ep_0_t*v1_z))*(v1_y-v2_y)-t31*(v1_x-v2_x)-t31*(w1_x-w2_x)+(t13*(t6+t10-ep_f_v2_x*t5-ep_0_t*v1_x)-t8*(t4+t18-ep_f_v2_y*t5-ep_0_t*v1_y))*(v1_z-v2_z)-(w1_z-w2_z)*(t8*t21-t13*t23);
		}
	}
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  return constVals;
}


void FoldingBinormalBiasConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {

	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0; int ijv_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		//for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			
			int v1,v2,w1,w2;
			dog->get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*wnum)); const double w1_z(x(w1_i+2*wnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*wnum)); const double w2_z(x(w2_i+2*wnum));

			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double t2 = ep_f_t-1.0;
			double t3 = ep_0_t-1.0;
			double t4 = ep_f_t*ep_f_v1_y;
			double t5 = t3*v2_y;
			double t21 = ep_0_t*v1_y;
			double t22 = ep_f_v2_y*t2;
			double t6 = t4+t5-t21-t22;
			double t7 = ep_f_t*ep_f_v1_z;
			double t8 = t3*v2_z;
			double t11 = ep_0_t*v1_z;
			double t12 = ep_f_v2_z*t2;
			double t9 = t7+t8-t11-t12;
			double t10 = v1_z-v2_z;
			double t13 = w1_z-w2_z;
			double t14 = ep_f_t*ep_f_v1_x;
			double t15 = t3*v2_x;
			double t18 = ep_0_t*v1_x;
			double t19 = ep_f_v2_x*t2;
			double t16 = t14+t15-t18-t19;
			double t17 = v1_y-v2_y;
			double t20 = v1_x-v2_x;
			double t23 = w1_y-w2_y;
			double t24 = w1_x-w2_x;
			double t25 = ep_b_t-1.0;
			double t26 = ep_b_t*ep_b_v1_y;
			double t34 = ep_b_v2_y*t25;
			double t27 = t5-t21+t26-t34;
			double t28 = ep_b_t*ep_b_v1_z;
			double t30 = ep_b_v2_z*t25;
			double t29 = t8-t11+t28-t30;
			double t31 = ep_b_t*ep_b_v1_x;
			double t33 = ep_b_v2_x*t25;
			double t32 = t15-t18+t31-t33;
			double t35 = ep_0_t*t6;
			double t36 = ep_0_t*t9;
			double t38 = ep_0_t*t29;
			double t37 = t36-t38;
			double t39 = ep_0_t*t16;
			double t41 = ep_0_t*t32;
			double t40 = t39-t41;
			double t43 = ep_0_t*t27;
			double t42 = t35-t43;
			double t44 = t6*t29;
			double t45 = t3*t6;
			double t54 = t3*t27;
			double t46 = t45-t54;
			double t47 = t3*t9;
			double t48 = t9*t32;
			double t51 = t3*t29;
			double t49 = t47-t51;
			double t50 = t3*t16;
			double t52 = t16*t27;
			double t55 = t3*t32;
			double t53 = t50-t55;
			double t56 = t9*t27;
			double t57 = t16*t29;
			double t58 = t6*t32;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i,-ep_b_t*t6*t10-ep_b_t*t6*t13+ep_b_t*t9*t17+ep_b_t*t9*t23);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+vnum,ep_b_t*t10*t16-ep_b_t*t9*t20+ep_b_t*t13*t16-ep_b_t*t9*t24);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+2*vnum,ep_b_t*t6*t20+ep_b_t*t6*t24-ep_b_t*t16*t17-ep_b_t*t16*t23);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i,t6*t10*t25+t6*t13*t25-t9*t17*t25-t9*t23*t25);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+vnum,-t10*t16*t25+t9*t20*t25-t13*t16*t25+t9*t24*t25);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+2*vnum,-t6*t20*t25-t6*t24*t25+t16*t17*t25+t16*t23*t25);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i,ep_f_t*t10*t27+ep_f_t*t13*t27-ep_f_t*t17*t29-ep_f_t*t23*t29);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+vnum,-ep_f_t*t10*t32-ep_f_t*t13*t32+ep_f_t*t20*t29+ep_f_t*t24*t29);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+2*vnum,-ep_f_t*t20*t27+ep_f_t*t17*t32-ep_f_t*t24*t27+ep_f_t*t23*t32);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i,-t2*t10*t27-t2*t13*t27+t2*t17*t29+t2*t23*t29);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+vnum,t2*t10*t32+t2*t13*t32-t2*t20*t29-t2*t24*t29);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+2*vnum,t2*t20*t27-t2*t17*t32+t2*t24*t27-t2*t23*t32);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i,t44-t9*t27+t10*t42-t17*t37+t13*t42-t23*t37);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+vnum,t48-t16*t29-t10*t40-t13*t40+t20*(t36-t38)+t24*(t36-t38));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+2*vnum,t52-t6*t32-t20*t42-t24*t42+t17*(t39-t41)+t23*(t39-t41));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i,-t44+t56-t10*t46-t13*t46+t17*t49+t23*t49);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+vnum,-t48+t57+t10*t53+t13*t53-t20*t49-t24*t49);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+2*vnum,-t52+t58-t17*t53-t23*t53+t20*(t45-t54)+t24*(t45-t54));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i,t44-t56);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+vnum,t48-t57);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+2*vnum,t52-t58);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i,-t44+t56);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+vnum,-t48+t57);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+2*vnum,-t52+t58);
		}
	}
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}

void FoldingBinormalBiasConstraints::updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda_v) {
	int vnum = x.rows()/3; int ijv_idx = 0;
	int const_cnt = 0;

	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		//for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			bool flip_binormal_sign = (is_mountain[curve_i] ^ flip_binormal[curve_i][edge_idx-1]);
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if (flip_binormal_sign) {
				// Flip the edge direction
				std::swap(ep_b,ep_f);
			}
			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);

			int ep_0_v1_i(ep.edge.v1), ep_0_v2_i(ep.edge.v2); const double ep_0_t(ep.t);
			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
			const double ep_0_v2_x(x(ep_0_v2_i)); const double ep_0_v2_y(x(ep_0_v2_i+1*vnum)); const double ep_0_v2_z(x(ep_0_v2_i+2*vnum));
			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			for (int const_side_i = 0; const_side_i < 2 ; const_side_i++) {
				int folded_v_i = fold_v_indices[const_side_i];
				const double folded_v_x(x(folded_v_i)); const double folded_v_y(x(folded_v_i+1*vnum)); const double folded_v_z(x(folded_v_i+2*vnum));
				double lambda = lambda_v(const_cnt);

				double t2 = ep_0_t*ep_0_v1_z;
				double t3 = ep_0_t-1.0;
				double t4 = ep_b_t-1.0;
				double t5 = ep_b_v2_z*t4;
				double t7 = ep_0_v2_z*t3;
				double t21 = ep_b_t*ep_b_v1_z;
				double t6 = t2+t5-t7-t21;
				double t8 = ep_f_t-1.0;
				double t9 = ep_f_v2_z*t8;
				double t16 = ep_f_t*ep_f_v1_z;
				double t10 = t2-t7+t9-t16;
				double t11 = ep_0_t*ep_0_v1_y;
				double t12 = ep_b_v2_y*t4;
				double t14 = ep_0_v2_y*t3;
				double t22 = ep_b_t*ep_b_v1_y;
				double t13 = t11+t12-t14-t22;
				double t15 = ep_f_v2_y*t8;
				double t20 = ep_f_t*ep_f_v1_y;
				double t17 = t11-t14+t15-t20;
				double t18 = folded_v_z-t2+t7;
				double t19 = folded_v_y-t11+t14;
				double t23 = ep_0_t*t6;
				double t26 = ep_0_t*t10;
				double t24 = t23-t26;
				double t25 = ep_0_t*t13;
				double t27 = t3*t24;
				double t28 = t3*t6;
				double t94 = t3*t10;
				double t29 = t28-t94;
				double t95 = ep_0_t*t29;
				double t30 = t27-t95;
				double t31 = ep_0_t*ep_0_v1_x;
				double t32 = ep_b_v2_x*t4;
				double t34 = ep_0_v2_x*t3;
				double t50 = ep_b_t*ep_b_v1_x;
				double t33 = t31+t32-t34-t50;
				double t35 = ep_f_v2_x*t8;
				double t40 = ep_f_t*ep_f_v1_x;
				double t36 = t31-t34+t35-t40;
				double t37 = ep_0_t*ep_b_t*t18;
				double t38 = ep_0_t*ep_b_t*t10;
				double t39 = t37+t38;
				double t41 = ep_0_t*t4*t18;
				double t42 = ep_0_t*t4*t10;
				double t43 = t41+t42;
				double t44 = lambda*t43;
				double t45 = folded_v_x-t31+t34;
				double t46 = ep_0_t*ep_f_t*t18;
				double t47 = ep_0_t*ep_f_t*t6;
				double t48 = t46+t47;
				double t49 = lambda*t48;
				double t51 = ep_0_t*t8*t18;
				double t52 = ep_0_t*t6*t8;
				double t53 = t51+t52;
				double t54 = ep_0_t*t33;
				double t59 = ep_0_t*t36;
				double t55 = t54-t59;
				double t92 = ep_0_t*t17;
				double t56 = t25-t92;
				double t57 = t3*t13;
				double t97 = t3*t17;
				double t58 = t57-t97;
				double t60 = t3*t55;
				double t61 = t3*t33;
				double t100 = t3*t36;
				double t62 = t61-t100;
				double t101 = ep_0_t*t62;
				double t63 = t60-t101;
				double t64 = ep_0_t*ep_b_t*t19;
				double t65 = ep_0_t*ep_b_t*t17;
				double t66 = t64+t65;
				double t67 = lambda*t66;
				double t68 = ep_0_t*ep_b_t*t45;
				double t69 = ep_0_t*ep_b_t*t36;
				double t70 = t68+t69;
				double t71 = ep_0_t*t4*t19;
				double t72 = ep_0_t*t4*t17;
				double t73 = t71+t72;
				double t74 = ep_0_t*t4*t45;
				double t75 = ep_0_t*t4*t36;
				double t76 = t74+t75;
				double t77 = lambda*t76;
				double t78 = ep_0_t*ep_f_t*t19;
				double t79 = ep_0_t*ep_f_t*t13;
				double t80 = t78+t79;
				double t81 = ep_0_t*ep_f_t*t45;
				double t82 = ep_0_t*ep_f_t*t33;
				double t83 = t81+t82;
				double t84 = lambda*t83;
				double t85 = ep_0_t*t8*t19;
				double t86 = ep_0_t*t8*t13;
				double t87 = t85+t86;
				double t88 = lambda*t87;
				double t89 = ep_0_t*t8*t45;
				double t90 = ep_0_t*t8*t33;
				double t91 = t89+t90;
				double t93 = lambda*t56;
				double t96 = lambda*t30;
				double t98 = ep_0_t*t58;
				double t118 = t3*t56;
				double t119 = t98-t118;
				double t99 = lambda*t119;
				double t102 = lambda*t63;
				double t103 = ep_b_t*t3*t18;
				double t104 = ep_b_t*t3*t10;
				double t105 = t103+t104;
				double t106 = lambda*t105;
				double t107 = t3*t4*t18;
				double t108 = t3*t4*t10;
				double t109 = t107+t108;
				double t110 = ep_f_t*t3*t18;
				double t111 = ep_f_t*t3*t6;
				double t112 = t110+t111;
				double t113 = t3*t8*t18;
				double t114 = t3*t6*t8;
				double t115 = t113+t114;
				double t116 = lambda*t115;
				double t117 = lambda*t29;
				double t120 = ep_b_t*t3*t19;
				double t121 = ep_b_t*t3*t17;
				double t122 = t120+t121;
				double t123 = ep_b_t*t3*t45;
				double t124 = ep_b_t*t3*t36;
				double t125 = t123+t124;
				double t126 = lambda*t125;
				double t127 = t3*t4*t19;
				double t128 = t3*t4*t17;
				double t129 = t127+t128;
				double t130 = lambda*t129;
				double t131 = t3*t4*t45;
				double t132 = t3*t4*t36;
				double t133 = t131+t132;
				double t134 = ep_f_t*t3*t19;
				double t135 = ep_f_t*t3*t13;
				double t136 = t134+t135;
				double t137 = lambda*t136;
				double t138 = ep_f_t*t3*t45;
				double t139 = ep_f_t*t3*t33;
				double t140 = t138+t139;
				double t141 = t3*t8*t19;
				double t142 = t3*t8*t13;
				double t143 = t141+t142;
				double t144 = t3*t8*t45;
				double t145 = t3*t8*t33;
				double t146 = t144+t145;
				double t147 = lambda*t146;
				double t148 = lambda*t62;
				double t149 = lambda*t39;
				double t150 = lambda*t122;
				double t151 = lambda*t70;
				double t152 = ep_b_t*lambda*t8*t18;
				double t153 = ep_b_t*ep_f_t*lambda*t19;
				double t154 = ep_b_t*lambda*t8*t45;
				double t155 = ep_b_t*lambda*t17;
				double t156 = lambda*t73;
				double t157 = lambda*t109;
				double t158 = lambda*t133;
				double t159 = ep_f_t*lambda*t4*t18;
				double t160 = lambda*t4*t10;
				double t161 = ep_f_t*lambda*t4*t45;
				double t162 = lambda*t4*t8*t19;
				double t163 = lambda*t4*t36;
				double t164 = lambda*t80;
				double t165 = lambda*t112;
				double t166 = ep_b_t*ep_f_t*lambda*t18;
				double t167 = ep_f_t*lambda*t4*t19;
				double t168 = lambda*t140;
				double t169 = ep_b_t*ep_f_t*lambda*t45;
				double t170 = ep_f_t*lambda*t6;
				double t171 = ep_f_t*lambda*t33;
				double t172 = lambda*t53;
				double t173 = lambda*t143;
				double t174 = ep_b_t*lambda*t8*t19;
				double t175 = lambda*t4*t8*t18;
				double t176 = lambda*t91;
				double t177 = lambda*t4*t8*t45;
				double t178 = lambda*t8*t13;
				double t179 = lambda*t24;
				double t180 = lambda*t58;
				double t181 = ep_b_t*lambda*t10;
				double t182 = lambda*t4*t17;
				double t183 = ep_f_t*lambda*t13;
				double t184 = lambda*t6*t8;
				double t185 = lambda*t55;
				double t186 = ep_b_t*lambda*t36;
				double t187 = lambda*t8*t33;

				//[ ep_0_v1_x, ep_0_v1_y, ep_0_v1_z, ep_0_v2_x, ep_0_v2_y, ep_0_v2_z, ep_b_v1_x, ep_b_v1_y, ep_b_v1_z, ep_b_v2_x, ep_b_v2_y, ep_b_v2_z, ep_f_v1_x, ep_f_v1_y, ep_f_v1_z, ep_f_v2_x, ep_f_v2_y, ep_f_v2_z, folded_v_x, folded_v_y, folded_v_z]
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_0_v2_i+vnum, -lambda*t30);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_0_v2_i+2*vnum, lambda*(t3*(t25-ep_0_t*(t11+t15-ep_f_t*ep_f_v1_y-ep_0_v2_y*t3))-ep_0_t*t58));
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_b_v1_i+vnum, -lambda*t39);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_b_v1_i+2*vnum, t67);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_b_v2_i+vnum, t44);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_b_v2_i+2*vnum, -lambda*t73);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_f_v1_i+vnum, t49);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_f_v1_i+2*vnum, -lambda*t80);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_f_v2_i+vnum, -lambda*t53);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,ep_f_v2_i+2*vnum, t88);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,folded_v_i+vnum, -lambda*t24);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i,folded_v_i+2*vnum, t93);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_0_v2_i, t96);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_0_v2_i+2*vnum, -lambda*t63);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_b_v1_i, t149);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_b_v1_i+2*vnum, -lambda*t70);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_b_v2_i, -t44);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_b_v2_i+2*vnum, t77);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_f_v1_i, -t49);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_f_v1_i+2*vnum, t84);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_f_v2_i, t172);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,ep_f_v2_i+2*vnum, -lambda*t91);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,folded_v_i, t179);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+vnum,folded_v_i+2*vnum, -lambda*t55);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_0_v2_i, t99);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_0_v2_i+vnum, t102);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_b_v1_i, -t67);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_b_v1_i+vnum, t151);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_b_v2_i, t156);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_b_v2_i+vnum, -t77);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_f_v1_i, t164);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_f_v1_i+vnum, -t84);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_f_v2_i, -t88);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,ep_f_v2_i+vnum, t176);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,folded_v_i, -t93);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v1_i+2*vnum,folded_v_i+vnum, t185);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_0_v1_i+vnum, t96);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_0_v1_i+2*vnum, t99);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_b_v1_i+vnum, t106);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_b_v1_i+2*vnum, -lambda*t122);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_b_v2_i+vnum, -lambda*t109);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_b_v2_i+2*vnum, t130);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_f_v1_i+vnum, -lambda*t112);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_f_v1_i+2*vnum, t137);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_f_v2_i+vnum, t116);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,ep_f_v2_i+2*vnum, -lambda*t143);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,folded_v_i+vnum, t117);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i,folded_v_i+2*vnum, -lambda*t58);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_0_v1_i, -t96);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_0_v1_i+2*vnum, t102);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_b_v1_i, -t106);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_b_v1_i+2*vnum, t126);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_b_v2_i, t157);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_b_v2_i+2*vnum, -lambda*t133);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_f_v1_i, t165);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_f_v1_i+2*vnum, -lambda*t140);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_f_v2_i, -t116);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,ep_f_v2_i+2*vnum, t147);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,folded_v_i, -t117);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+vnum,folded_v_i+2*vnum, t148);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_0_v1_i, -t99);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_0_v1_i+vnum, -t102);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_b_v1_i, t150);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_b_v1_i+vnum, -t126);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_b_v2_i, -t130);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_b_v2_i+vnum, t158);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_f_v1_i, -t137);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_f_v1_i+vnum, t168);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_f_v2_i, t173);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,ep_f_v2_i+vnum, -t147);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,folded_v_i, t180);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_0_v2_i+2*vnum,folded_v_i+vnum, -t148);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_0_v1_i+vnum, t149);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_0_v1_i+2*vnum, -t67);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_0_v2_i+vnum, -t106);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_0_v2_i+2*vnum, t150);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v1_i+vnum, -ep_b_t*ep_f_t*lambda*t18);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v1_i+2*vnum, t153);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v2_i+vnum, t152);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v2_i+2*vnum, -ep_b_t*lambda*t8*t19);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,folded_v_i+vnum, -ep_b_t*lambda*t10);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,folded_v_i+2*vnum, t155);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_0_v1_i, -t149);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_0_v1_i+2*vnum, t151);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_0_v2_i, t106);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_0_v2_i+2*vnum, -t126);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v1_i, t166);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v1_i+2*vnum, -ep_b_t*ep_f_t*lambda*t45);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v2_i, -t152);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v2_i+2*vnum, t154);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,folded_v_i, t181);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,folded_v_i+2*vnum, -ep_b_t*lambda*t36);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_0_v1_i, t67);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_0_v1_i+vnum, -t151);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_0_v2_i, -t150);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_0_v2_i+vnum, t126);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v1_i, -t153);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v1_i+vnum, t169);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v2_i, t174);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v2_i+vnum, -t154);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,folded_v_i, -t155);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,folded_v_i+vnum, t186);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_0_v1_i+vnum, -t44);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_0_v1_i+2*vnum, t156);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_0_v2_i+vnum, t157);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_0_v2_i+2*vnum, -t130);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v1_i+vnum, t159);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v1_i+2*vnum, -ep_f_t*lambda*t4*t19);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v2_i+vnum, -lambda*t4*t8*t18);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v2_i+2*vnum, t162);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,folded_v_i+vnum, t160);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,folded_v_i+2*vnum, -lambda*t4*t17);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_0_v1_i, t44);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_0_v1_i+2*vnum, -t77);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_0_v2_i, -t157);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_0_v2_i+2*vnum, t158);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v1_i, -t159);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v1_i+2*vnum, t161);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v2_i, t175);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v2_i+2*vnum, -lambda*t4*t8*t45);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,folded_v_i, -t160);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,folded_v_i+2*vnum, t163);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_0_v1_i, -t156);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_0_v1_i+vnum, t77);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_0_v2_i, t130);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_0_v2_i+vnum, -t158);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v1_i, t167);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v1_i+vnum, -t161);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v2_i, -t162);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v2_i+vnum, t177);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,folded_v_i, t182);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,folded_v_i+vnum, -t163);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_0_v1_i+vnum, -t49);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_0_v1_i+2*vnum, t164);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_0_v2_i+vnum, t165);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_0_v2_i+2*vnum, -t137);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v1_i+vnum, t166);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v1_i+2*vnum, -t153);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v2_i+vnum, -t159);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v2_i+2*vnum, t167);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,folded_v_i+vnum, t170);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,folded_v_i+2*vnum, -ep_f_t*lambda*t13);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_0_v1_i, t49);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_0_v1_i+2*vnum, -t84);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_0_v2_i, -t165);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_0_v2_i+2*vnum, t168);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v1_i, -t166);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v1_i+2*vnum, t169);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v2_i, t159);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v2_i+2*vnum, -t161);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,folded_v_i, -t170);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,folded_v_i+2*vnum, t171);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_0_v1_i, -t164);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_0_v1_i+vnum, t84);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_0_v2_i, t137);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_0_v2_i+vnum, -t168);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v1_i, t153);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v1_i+vnum, -t169);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v2_i, -t167);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v2_i+vnum, t161);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,folded_v_i, t183);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,folded_v_i+vnum, -t171);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_0_v1_i+vnum, t172);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_0_v1_i+2*vnum, -t88);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_0_v2_i+vnum, -t116);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_0_v2_i+2*vnum, t173);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v1_i+vnum, -t152);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v1_i+2*vnum, t174);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v2_i+vnum, t175);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v2_i+2*vnum, -t162);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,folded_v_i+vnum, -lambda*t6*t8);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,folded_v_i+2*vnum, t178);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_0_v1_i, -t172);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_0_v1_i+2*vnum, t176);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_0_v2_i, t116);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_0_v2_i+2*vnum, -t147);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v1_i, t152);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v1_i+2*vnum, -t154);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v2_i, -t175);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v2_i+2*vnum, t177);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,folded_v_i, t184);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,folded_v_i+2*vnum, -lambda*t8*t33);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_0_v1_i, t88);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_0_v1_i+vnum, -t176);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_0_v2_i, -t173);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_0_v2_i+vnum, t147);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v1_i, -t174);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v1_i+vnum, t154);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v2_i, t162);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v2_i+vnum, -t177);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,folded_v_i, -t178);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,folded_v_i+vnum, t187);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_0_v1_i+vnum, t179);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_0_v1_i+2*vnum, -t93);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_0_v2_i+vnum, -t117);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_0_v2_i+2*vnum, t180);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_b_v1_i+vnum, t181);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_b_v1_i+2*vnum, -t155);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_b_v2_i+vnum, -t160);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_b_v2_i+2*vnum, t182);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_f_v1_i+vnum, -t170);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_f_v1_i+2*vnum, t183);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_f_v2_i+vnum, t184);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i,ep_f_v2_i+2*vnum, -t178);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_0_v1_i, -t179);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_0_v1_i+2*vnum, t185);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_0_v2_i, t117);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_0_v2_i+2*vnum, -t148);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_b_v1_i, -t181);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_b_v1_i+2*vnum, t186);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_b_v2_i, t160);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_b_v2_i+2*vnum, -t163);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_f_v1_i, t170);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_f_v1_i+2*vnum, -t171);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_f_v2_i, -t184);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+vnum,ep_f_v2_i+2*vnum, t187);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_0_v1_i, t93);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_0_v1_i+vnum, -t185);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_0_v2_i, -t180);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_0_v2_i+vnum, t148);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_b_v1_i, t155);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_b_v1_i+vnum, -t186);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_b_v2_i, -t182);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_b_v2_i+vnum, t163);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_f_v1_i, -t183);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_f_v1_i+vnum, t171);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_f_v2_i, t178);
				lambda_hessian_IJV[ijv_idx++] = Eigen::Triplet<double>(folded_v_i+2*vnum,ep_f_v2_i+vnum, -t187);

				const_cnt++;
			}
		}
	}
	if (ijv_idx!= lambda_hessian_IJV.size()) {
		cout << "error, ijv_idx " << ijv_idx << " but lambda_hessian_IJV.size() = " << lambda_hessian_IJV.size() << endl;
		exit(1);		
	}
};