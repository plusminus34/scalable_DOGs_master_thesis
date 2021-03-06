#include "FoldingBinormalBiasConstraints.h"

using namespace std;

FoldingBinormalBiasConstraints::FoldingBinormalBiasConstraints(const Dog& dog) : 
	dog(dog), eS(dog.getEdgeStitching()), vnum(dog.getQuadTopology().v_n) {
    // TODO: different handling of vertices (maybe no consatraints there?)
	//const_n= 2*(eS.edge_coordinates.size()-2);
	//const_n = 1;//*(eS.edge_coordinates.size()-2);
	std::cout << "number of curves = " << eS.stitched_curves.size() << std::endl;
	const_n = 0;
	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
			EdgePoint ep = foldingCurve[edge_idx];
			if ((eS.get_vertex_edge_point_deg(ep.edge) == 1) && !dog.is_crease_vertex_flat(curve_i,edge_idx) ) const_n++;
			else {
				if (dog.is_crease_vertex_flat(curve_i,edge_idx)) {
					std::cout << "curve_i = " << curve_i << " edge_idx = " << edge_idx << " is flat = " <<  dog.is_crease_vertex_flat(curve_i,edge_idx) << endl;
				}
				if (eS.get_vertex_edge_point_deg(ep.edge) != 1) {
					std::cout << "curve_i = " << curve_i << " edge_idx = " << edge_idx << " is a vertex" << " with rank = " << eS.get_vertex_edge_point_deg(ep.edge) + 1 << std::endl;	
				}
			}
		}
	}
	//const_n = eS.edge_coordinates.size()-2*eS.stitched_curves.size();
	std::cout << "const_n = " << const_n << " and eS.edge_coordinates.size()-2*eS.stitched_curves.size() = " << eS.edge_coordinates.size()-2*eS.stitched_curves.size()  << std::endl;
	//exit(1);
	IJV.resize(24*const_n);
	//lambda_hessian_IJV.resize(300*const_n);
}

Eigen::VectorXd FoldingBinormalBiasConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int const_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		//for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if ((eS.get_vertex_edge_point_deg(ep.edge) != 1) || dog.is_crease_vertex_flat(curve_i,edge_idx) ) continue;
			int v1,v2,w1,w2;
			dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

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
			double t23 = ep_b_v2_x*t3;
			double t8 = t6+t7-t11-t23;
			double t9 = t2*v2_z;
			double t10 = ep_f_t*ep_f_v1_x;
			double t12 = ep_b_t*ep_b_v1_y;
			double t19 = ep_0_t*v1_y;
			double t26 = ep_b_v2_y*t3;
			double t13 = t4+t12-t19-t26;
			double t14 = ep_f_t*ep_f_v1_z;
			double t17 = ep_0_t*v1_z;
			double t29 = ep_f_v2_z*t5;
			double t15 = t9+t14-t17-t29;
			double t16 = ep_b_t*ep_b_v1_z;
			double t18 = ep_f_t*ep_f_v1_y;
			double t27 = ep_f_v2_x*t5;
			double t20 = t6+t10-t11-t27;
			double t31 = ep_b_v2_z*t3;
			double t21 = t9+t16-t17-t31;
			double t24 = ep_f_v2_y*t5;
			double t22 = t4+t18-t19-t24;
			double t25 = t8*t22;
			double t28 = t25-t13*t20;
			double t30 = t8*t15;
			double t32 = t30-t20*t21;
			double t33 = t13*t15;
			double t34 = t33-t21*t22;
			constVals(const_cnt++) = -tanh(t28*(v1_z-v2_z)*1.0E3-t32*(v1_y-v2_y)*1.0E3+t34*(v1_x-v2_x)*1.0E3)-tanh(t28*(w1_z-w2_z)*1.0E3-t32*(w1_y-w2_y)*1.0E3+t34*(w1_x-w2_x)*1.0E3);
		}
	}
  if (const_cnt != const_n) {
		cout << "error in Vals, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  //std::cout << "constVals.norm() = " << constVals.norm() << std::endl;
  return constVals;
}

void FoldingBinormalBiasConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {
	// Add curve fold constraints
	int const_cnt = 0; int ijv_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		//for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if ((eS.get_vertex_edge_point_deg(ep.edge) != 1) || dog.is_crease_vertex_flat(curve_i,edge_idx) ) continue;
			
			int v1,v2,w1,w2;
			dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));

			double t2 = ep_f_t-1.0;
			double t3 = ep_0_t-1.0;
			double t4 = v1_z-v2_z;
			double t5 = ep_f_t*ep_f_v1_y;
			double t6 = t3*v2_y;
			double t8 = ep_0_t*v1_y;
			double t25 = ep_f_v2_y*t2;
			double t7 = t5+t6-t8-t25;
			double t9 = ep_b_t-1.0;
			double t10 = t3*v2_x;
			double t11 = v1_y-v2_y;
			double t12 = ep_b_t*ep_b_v1_x;
			double t19 = ep_0_t*v1_x;
			double t28 = ep_b_v2_x*t9;
			double t13 = t10+t12-t19-t28;
			double t14 = ep_f_t*ep_f_v1_z;
			double t15 = t3*v2_z;
			double t17 = ep_0_t*v1_z;
			double t22 = ep_f_v2_z*t2;
			double t16 = t14+t15-t17-t22;
			double t18 = ep_f_t*ep_f_v1_x;
			double t20 = ep_b_t*ep_b_v1_y;
			double t30 = ep_b_v2_y*t9;
			double t21 = t6-t8+t20-t30;
			double t23 = ep_b_t*ep_b_v1_z;
			double t35 = ep_b_v2_z*t9;
			double t24 = t15-t17+t23-t35;
			double t29 = t7*t13;
			double t31 = ep_f_v2_x*t2;
			double t32 = t10+t18-t19-t31;
			double t34 = t13*t16;
			double t36 = t24*t32;
			double t37 = t34-t36;
			double t38 = t16*t21;
			double t39 = t7*t24;
			double t40 = t38-t39;
			double t42 = v1_x-v2_x;
			double t45 = t11*t37*1.0E3;
			double t46 = t40*t42*1.0E3;
			double t26 = tanh(-t45+t46+t4*(t29-t21*(t10+t18-t31-ep_0_t*v1_x))*1.0E3);
			double t27 = w1_z-w2_z;
			double t33 = w1_y-w2_y;
			double t43 = t21*t32;
			double t44 = t29-t43;
			double t48 = w1_x-w2_x;
			double t49 = t27*t44*1.0E3;
			double t50 = t33*t37*1.0E3;
			double t51 = t40*t48*1.0E3;
			double t52 = t49-t50+t51;
			double t41 = tanh(t52);
			double t55 = t4*t44*1.0E3;
			double t56 = -t45+t46+t55;
			double t47 = tanh(t56);
			double t53 = t41*t41;
			double t54 = t53-1.0;
			double t57 = t47*t47;
			double t58 = t57-1.0;
			double t59 = ep_0_t*t7;
			double t67 = ep_0_t*t21;
			double t60 = t59-t67;
			double t61 = ep_0_t*t16;
			double t63 = ep_0_t*t24;
			double t62 = t61-t63;
			double t64 = ep_0_t*t13;
			double t66 = ep_0_t*t32;
			double t65 = t64-t66;
			double t68 = t7*t24*1.0E3;
			double t69 = t3*t7;
			double t77 = t3*t21;
			double t70 = t69-t77;
			double t71 = t3*t16;
			double t72 = t13*t16*1.0E3;
			double t73 = t3*t13;
			double t76 = t3*t32;
			double t74 = t73-t76;
			double t75 = t7*t13*1.0E3;
			double t78 = t72-t24*t32*1.0E3;
			double t79 = t75-t21*t32*1.0E3;
			double t80 = t54*t79;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i,(t26*t26-1.0)*(ep_b_t*t4*t7*1.0E3-ep_b_t*t11*t16*1.0E3)+t54*(ep_b_t*t7*t27*1.0E3-ep_b_t*t16*t33*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+vnum,-t58*(ep_b_t*t4*t32*1.0E3-ep_b_t*t16*t42*1.0E3)-t54*(ep_b_t*t27*t32*1.0E3-ep_b_t*t16*t48*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+2*vnum,t58*(ep_b_t*t11*t32*1.0E3-ep_b_t*t7*t42*1.0E3)-t54*(ep_b_t*t7*t48*1.0E3-ep_b_t*t32*t33*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i,-t58*(t4*t7*t9*1.0E3-t9*t11*t16*1.0E3)-t54*(t7*t9*t27*1.0E3-t9*t16*t33*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+vnum,t58*(t4*t9*t32*1.0E3-t9*t16*t42*1.0E3)+t54*(t9*t27*t32*1.0E3-t9*t16*t48*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+2*vnum,-t58*(t9*t11*t32*1.0E3-t7*t9*t42*1.0E3)+t54*(t7*t9*t48*1.0E3-t9*t32*t33*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i,-t58*(ep_f_t*t4*t21*1.0E3-ep_f_t*t11*t24*1.0E3)-t54*(ep_f_t*t21*t27*1.0E3-ep_f_t*t24*t33*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+vnum,t58*(ep_f_t*t4*t13*1.0E3-ep_f_t*t24*t42*1.0E3)+t54*(ep_f_t*t13*t27*1.0E3-ep_f_t*t24*t48*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+2*vnum,-t58*(ep_f_t*t11*t13*1.0E3-ep_f_t*t21*t42*1.0E3)-t54*(ep_f_t*t13*t33*1.0E3-ep_f_t*t21*t48*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i,t58*(t2*t4*t21*1.0E3-t2*t11*t24*1.0E3)+t54*(t2*t21*t27*1.0E3-t2*t24*t33*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+vnum,-t58*(t2*t4*t13*1.0E3-t2*t24*t42*1.0E3)-t54*(t2*t13*t27*1.0E3-t2*t24*t48*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+2*vnum,t58*(t2*t11*t13*1.0E3-t2*t21*t42*1.0E3)+t54*(t2*t13*t33*1.0E3-t2*t21*t48*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i,-t58*(t68-t16*t21*1.0E3+t4*t60*1.0E3-t11*t62*1.0E3)-t54*(t27*t60*1.0E3-t33*t62*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+vnum,-t58*(t72-t24*t32*1.0E3+t4*t65*1.0E3+t42*t62*1.0E3)-t54*(t27*t65*1.0E3+t48*t62*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+2*vnum,t58*(t75-t21*t32*1.0E3+t11*t65*1.0E3+t42*(t59-t67)*1.0E3)+t54*(t33*t65*1.0E3+t48*(t59-t67)*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i,t58*(t68-t16*t21*1.0E3+t4*t70*1.0E3-t11*(t71-t3*t24)*1.0E3)+t54*(t27*t70*1.0E3-t33*(t71-t3*t24)*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+vnum,t58*(t72-t24*t32*1.0E3+t4*t74*1.0E3+t42*(t71-t3*t24)*1.0E3)+t54*(t27*t74*1.0E3+t48*(t71-t3*t24)*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+2*vnum,-t58*(t75-t21*t32*1.0E3+t11*t74*1.0E3+t42*t70*1.0E3)-t54*(t33*t74*1.0E3+t48*t70*1.0E3));

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i,-t54*(t68-t16*t21*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+vnum,-t54*t78);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+2*vnum,t80);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i,t54*(t68-t16*t21*1.0E3));
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+vnum,t54*t78);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+2*vnum,-t80);

			const_cnt++;
		}
	}
  if (const_cnt != const_n) {
		cout << "error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}
/*
void FoldingBinormalBiasConstraints::updateJacobianIJV_old(const Eigen::VectorXd& x) {
	// Add curve fold constraints
	int const_cnt = 0; int ijv_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		//for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if ((eS.get_vertex_edge_point_deg(ep.edge) != 1) || dog.is_crease_vertex_flat(curve_i,edge_idx) ) continue;
			
			int v1,v2,w1,w2;
			dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

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

			const_cnt++;
		}
	}
  if (const_cnt != const_n) {
		cout << "error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}


Eigen::VectorXd FoldingBinormalBiasConstraints::Vals_old(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int const_cnt = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		//for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if ((eS.get_vertex_edge_point_deg(ep.edge) != 1) || dog.is_crease_vertex_flat(curve_i,edge_idx) ) continue;
			int v1,v2,w1,w2;
			dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

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
		cout << "error in Vals, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  return constVals;
}

*/


void FoldingBinormalBiasConstraints::updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda_v) {
		// Add curve fold constraints
	int const_cnt = 0; int ijv_idx = 0;
  	for (int curve_i = 0; curve_i < eS.stitched_curves.size(); curve_i++) {
		const vector<EdgePoint>& foldingCurve = eS.stitched_curves[curve_i];
		// Go through all the inner vertices in a curve
		for (int edge_idx = 1; edge_idx < eS.stitched_curves[curve_i].size()-1; edge_idx++) {
		//for (int edge_idx = 1; edge_idx < 2; edge_idx++) {
			// Should flip the binormal only if is_mountain xor flip_binormal is true
			EdgePoint ep = foldingCurve[edge_idx], ep_b = foldingCurve[edge_idx-1], ep_f = foldingCurve[edge_idx+1];
			if ((eS.get_vertex_edge_point_deg(ep.edge) != 1) || dog.is_crease_vertex_flat(curve_i,edge_idx) ) continue;
			
			int v1,v2,w1,w2;
			dog.get_2_submeshes_vertices_from_edge(ep.edge, v1,v2,w1,w2);

			int fold_v_indices[2]; dog.get_2_inner_vertices_from_edge(ep.edge,fold_v_indices[0],fold_v_indices[1]);
			int v1_i(v1), v2_i(v2), w1_i(w1), w2_i(w2); const double ep_0_t(ep.t);

			int ep_b_v1_i(ep_b.edge.v1), ep_b_v2_i(ep_b.edge.v2); const double ep_b_t(ep_b.t);
			int ep_f_v1_i(ep_f.edge.v1), ep_f_v2_i(ep_f.edge.v2); const double ep_f_t(ep_f.t);

			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

			const double ep_b_v1_x(x(ep_b_v1_i)); const double ep_b_v1_y(x(ep_b_v1_i+1*vnum)); const double ep_b_v1_z(x(ep_b_v1_i+2*vnum));
			const double ep_b_v2_x(x(ep_b_v2_i)); const double ep_b_v2_y(x(ep_b_v2_i+1*vnum)); const double ep_b_v2_z(x(ep_b_v2_i+2*vnum));
			const double ep_f_v1_x(x(ep_f_v1_i)); const double ep_f_v1_y(x(ep_f_v1_i+1*vnum)); const double ep_f_v1_z(x(ep_f_v1_i+2*vnum));
			const double ep_f_v2_x(x(ep_f_v2_i)); const double ep_f_v2_y(x(ep_f_v2_i+1*vnum)); const double ep_f_v2_z(x(ep_f_v2_i+2*vnum));
			double lambda = lambda_v(const_cnt);

			double t2 = v1_z-v2_z;
			double t3 = w1_z-w2_z;
			double t4 = ep_f_t-1.0;
			double t5 = v1_y-v2_y;
			double t6 = w1_y-w2_y;
			double t7 = ep_0_t-1.0;
			double t8 = ep_f_t*ep_f_v1_z;
			double t9 = t7*v2_z;
			double t16 = ep_0_t*v1_z;
			double t17 = ep_f_v2_z*t4;
			double t10 = t8+t9-t16-t17;
			double t11 = ep_b_t*t10;
			double t12 = ep_f_t*ep_f_v1_y;
			double t13 = t7*v2_y;
			double t18 = ep_0_t*v1_y;
			double t19 = ep_f_v2_y*t4;
			double t14 = t12+t13-t18-t19;
			double t15 = ep_b_t*t14;
			double t20 = ep_b_t*lambda*t10;
			double t21 = ep_b_t*ep_f_t*t2;
			double t22 = ep_b_t*ep_f_t*t3;
			double t23 = t21+t22;
			double t24 = ep_b_t*t2*t4;
			double t25 = ep_b_t*t3*t4;
			double t26 = t24+t25;
			double t27 = lambda*t26;
			double t28 = v1_x-v2_x;
			double t29 = w1_x-w2_x;
			double t30 = ep_0_t*ep_b_t*t2;
			double t31 = ep_0_t*ep_b_t*t3;
			double t32 = t11+t30+t31;
			double t33 = lambda*t32;
			double t34 = ep_b_t*t2*t7;
			double t35 = ep_b_t*t3*t7;
			double t36 = t11+t34+t35;
			double t37 = ep_f_t*ep_f_v1_x;
			double t38 = t7*v2_x;
			double t41 = ep_0_t*v1_x;
			double t42 = ep_f_v2_x*t4;
			double t39 = t37+t38-t41-t42;
			double t40 = ep_b_t*t39;
			double t43 = ep_b_t*lambda*t39;
			double t44 = ep_b_t*ep_f_t*t5;
			double t45 = ep_b_t*ep_f_t*t6;
			double t46 = t44+t45;
			double t47 = lambda*t46;
			double t48 = ep_b_t*ep_f_t*t28;
			double t49 = ep_b_t*ep_f_t*t29;
			double t50 = t48+t49;
			double t51 = ep_b_t*t4*t5;
			double t52 = ep_b_t*t4*t6;
			double t53 = t51+t52;
			double t54 = ep_b_t*t4*t28;
			double t55 = ep_b_t*t4*t29;
			double t56 = t54+t55;
			double t57 = lambda*t56;
			double t58 = ep_0_t*ep_b_t*t5;
			double t59 = ep_0_t*ep_b_t*t6;
			double t60 = t15+t58+t59;
			double t61 = ep_0_t*ep_b_t*t28;
			double t62 = ep_0_t*ep_b_t*t29;
			double t63 = t40+t61+t62;
			double t64 = lambda*t63;
			double t65 = ep_b_t*t5*t7;
			double t66 = ep_b_t*t6*t7;
			double t67 = t15+t65+t66;
			double t68 = lambda*t67;
			double t69 = ep_b_t*t7*t28;
			double t70 = ep_b_t*t7*t29;
			double t71 = t40+t69+t70;
			double t72 = ep_b_t*lambda*t14;
			double t73 = ep_b_t-1.0;
			double t74 = t10*t73;
			double t75 = t14*t73;
			double t76 = lambda*t14*t73;
			double t77 = ep_f_t*t2*t73;
			double t78 = ep_f_t*t3*t73;
			double t79 = t77+t78;
			double t80 = lambda*t79;
			double t81 = t2*t4*t73;
			double t82 = t3*t4*t73;
			double t83 = t81+t82;
			double t84 = ep_0_t*t2*t73;
			double t85 = ep_0_t*t3*t73;
			double t86 = t74+t84+t85;
			double t87 = t2*t7*t73;
			double t88 = t3*t7*t73;
			double t89 = t74+t87+t88;
			double t90 = lambda*t89;
			double t91 = t39*t73;
			double t92 = lambda*t10*t73;
			double t93 = ep_f_t*t5*t73;
			double t94 = ep_f_t*t6*t73;
			double t95 = t93+t94;
			double t96 = ep_f_t*t28*t73;
			double t97 = ep_f_t*t29*t73;
			double t98 = t96+t97;
			double t99 = lambda*t98;
			double t100 = t4*t5*t73;
			double t101 = t4*t6*t73;
			double t102 = t100+t101;
			double t103 = lambda*t102;
			double t104 = t4*t28*t73;
			double t105 = t4*t29*t73;
			double t106 = t104+t105;
			double t107 = ep_0_t*t5*t73;
			double t108 = ep_0_t*t6*t73;
			double t109 = t75+t107+t108;
			double t110 = lambda*t109;
			double t111 = ep_0_t*t28*t73;
			double t112 = ep_0_t*t29*t73;
			double t113 = t91+t111+t112;
			double t114 = t5*t7*t73;
			double t115 = t6*t7*t73;
			double t116 = t75+t114+t115;
			double t117 = t7*t28*t73;
			double t118 = t7*t29*t73;
			double t119 = t91+t117+t118;
			double t120 = lambda*t119;
			double t121 = lambda*t39*t73;
			double t122 = lambda*t23;
			double t123 = lambda*t95;
			double t124 = ep_b_t*ep_b_v1_z;
			double t130 = ep_b_v2_z*t73;
			double t125 = t9-t16+t124-t130;
			double t126 = ep_f_t*t125;
			double t127 = ep_b_t*ep_b_v1_y;
			double t131 = ep_b_v2_y*t73;
			double t128 = t13-t18+t127-t131;
			double t129 = ep_f_t*t128;
			double t132 = ep_f_t*lambda*t128;
			double t133 = lambda*t50;
			double t134 = ep_0_t*ep_f_t*t2;
			double t135 = ep_0_t*ep_f_t*t3;
			double t136 = t126+t134+t135;
			double t137 = ep_f_t*t2*t7;
			double t138 = ep_f_t*t3*t7;
			double t139 = t126+t137+t138;
			double t140 = lambda*t139;
			double t141 = ep_b_t*ep_b_v1_x;
			double t145 = ep_b_v2_x*t73;
			double t142 = t38-t41+t141-t145;
			double t143 = ep_f_t*t142;
			double t144 = ep_f_t*lambda*t125;
			double t146 = ep_0_t*ep_f_t*t5;
			double t147 = ep_0_t*ep_f_t*t6;
			double t148 = t129+t146+t147;
			double t149 = lambda*t148;
			double t150 = ep_0_t*ep_f_t*t28;
			double t151 = ep_0_t*ep_f_t*t29;
			double t152 = t143+t150+t151;
			double t153 = ep_f_t*t5*t7;
			double t154 = ep_f_t*t6*t7;
			double t155 = t129+t153+t154;
			double t156 = ep_f_t*t7*t28;
			double t157 = ep_f_t*t7*t29;
			double t158 = t143+t156+t157;
			double t159 = lambda*t158;
			double t160 = ep_f_t*lambda*t142;
			double t161 = lambda*t53;
			double t162 = lambda*t83;
			double t163 = t4*t125;
			double t164 = t4*t128;
			double t165 = lambda*t4*t125;
			double t166 = lambda*t106;
			double t167 = ep_0_t*t2*t4;
			double t168 = ep_0_t*t3*t4;
			double t169 = t163+t167+t168;
			double t170 = lambda*t169;
			double t171 = t2*t4*t7;
			double t172 = t3*t4*t7;
			double t173 = t163+t171+t172;
			double t174 = t4*t142;
			double t175 = lambda*t4*t142;
			double t176 = ep_0_t*t4*t5;
			double t177 = ep_0_t*t4*t6;
			double t178 = t164+t176+t177;
			double t179 = ep_0_t*t4*t28;
			double t180 = ep_0_t*t4*t29;
			double t181 = t174+t179+t180;
			double t182 = lambda*t181;
			double t183 = t4*t5*t7;
			double t184 = t4*t6*t7;
			double t185 = t164+t183+t184;
			double t186 = lambda*t185;
			double t187 = t4*t7*t28;
			double t188 = t4*t7*t29;
			double t189 = t174+t187+t188;
			double t190 = lambda*t4*t128;
			double t191 = lambda*t60;
			double t192 = lambda*t86;
			double t193 = lambda*t136;
			double t194 = lambda*t178;
			double t195 = ep_0_t*t10;
			double t196 = ep_0_t*t14;
			double t209 = ep_0_t*t128;
			double t197 = t196-t209;
			double t198 = lambda*t197;
			double t199 = lambda*t113;
			double t200 = lambda*t152;
			double t201 = t7*t125;
			double t204 = ep_0_t*t125;
			double t220 = t7*t10;
			double t202 = t195+t201-t204-t220;
			double t203 = lambda*t202;
			double t205 = ep_0_t*t39;
			double t206 = t195-t204;
			double t207 = lambda*t206;
			double t208 = t7*t128;
			double t210 = t7*t142;
			double t213 = ep_0_t*t142;
			double t229 = t7*t39;
			double t211 = t205+t210-t213-t229;
			double t212 = lambda*t211;
			double t214 = t205-t213;
			double t215 = lambda*t214;
			double t216 = lambda*t36;
			double t217 = lambda*t116;
			double t218 = lambda*t155;
			double t219 = lambda*t173;
			double t223 = t7*t14;
			double t221 = t196+t208-t209-t223;
			double t222 = lambda*t221;
			double t224 = t201-t220;
			double t225 = t208-t223;
			double t226 = lambda*t225;
			double t227 = lambda*t71;
			double t228 = lambda*t189;
			double t230 = lambda*t224;
			double t231 = t210-t229;
			double t232 = lambda*t231;
			double t233 = lambda*(t196-t209);
			double t234 = lambda*(t195-t204);
			double t235 = lambda*(t205-t213);

			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v1_i+vnum, -lambda*t23);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v1_i+2*vnum, t47);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v2_i+vnum, t27);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,ep_f_v2_i+2*vnum, -lambda*t53);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,v1_i+vnum, t33);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,v1_i+2*vnum, -lambda*t60);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,v2_i+vnum, -lambda*t36);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,v2_i+2*vnum, t68);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,w1_i+vnum, t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,w1_i+2*vnum, -ep_b_t*lambda*t14);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,w2_i+vnum, -t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i,w2_i+2*vnum, t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v1_i, t122);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v1_i+2*vnum, -lambda*t50);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v2_i, -t27);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,ep_f_v2_i+2*vnum, t57);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,v1_i, -t33);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,v1_i+2*vnum, t64);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,v2_i, t216);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,v2_i+2*vnum, -lambda*t71);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,w1_i, -t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,w1_i+2*vnum, t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,w2_i, t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+vnum,w2_i+2*vnum, -t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v1_i, -t47);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v1_i+vnum, t133);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v2_i, t161);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,ep_f_v2_i+vnum, -t57);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,v1_i, t191);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,v1_i+vnum, -t64);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,v2_i, -t68);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,v2_i+vnum, t227);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,w1_i, t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,w1_i+vnum, -t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,w2_i, -t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v1_i+2*vnum,w2_i+vnum, t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v1_i+vnum, t80);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v1_i+2*vnum, -lambda*t95);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v2_i+vnum, -lambda*t83);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,ep_f_v2_i+2*vnum, t103);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,v1_i+vnum, -lambda*t86);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,v1_i+2*vnum, t110);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,v2_i+vnum, t90);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,v2_i+2*vnum, -lambda*t116);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,w1_i+vnum, -lambda*t10*t73);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,w1_i+2*vnum, t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,w2_i+vnum, t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i,w2_i+2*vnum, -t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v1_i, -t80);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v1_i+2*vnum, t99);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v2_i, t162);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,ep_f_v2_i+2*vnum, -lambda*t106);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,v1_i, t192);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,v1_i+2*vnum, -lambda*t113);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,v2_i, -t90);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,v2_i+2*vnum, t120);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,w1_i, t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,w1_i+2*vnum, -lambda*t39*t73);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,w2_i, -t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+vnum,w2_i+2*vnum, t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v1_i, t123);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v1_i+vnum, -t99);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v2_i, -t103);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,ep_f_v2_i+vnum, t166);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,v1_i, -t110);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,v1_i+vnum, t199);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,v2_i, t217);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,v2_i+vnum, -t120);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,w1_i, -t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,w1_i+vnum, t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,w2_i, t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_b_v2_i+2*vnum,w2_i+vnum, -t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v1_i+vnum, t122);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v1_i+2*vnum, -t47);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v2_i+vnum, -t80);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,ep_b_v2_i+2*vnum, t123);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,v1_i+vnum, -lambda*t136);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,v1_i+2*vnum, t149);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,v2_i+vnum, t140);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,v2_i+2*vnum, -lambda*t155);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,w1_i+vnum, -ep_f_t*lambda*t125);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,w1_i+2*vnum, t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,w2_i+vnum, t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i,w2_i+2*vnum, -t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v1_i, -t122);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v1_i+2*vnum, t133);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v2_i, t80);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,ep_b_v2_i+2*vnum, -t99);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,v1_i, t193);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,v1_i+2*vnum, -lambda*t152);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,v2_i, -t140);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,v2_i+2*vnum, t159);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,w1_i, t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,w1_i+2*vnum, -ep_f_t*lambda*t142);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,w2_i, -t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+vnum,w2_i+2*vnum, t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v1_i, t47);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v1_i+vnum, -t133);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v2_i, -t123);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,ep_b_v2_i+vnum, t99);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,v1_i, -t149);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,v1_i+vnum, t200);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,v2_i, t218);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,v2_i+vnum, -t159);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,w1_i, -t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,w1_i+vnum, t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,w2_i, t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v1_i+2*vnum,w2_i+vnum, -t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v1_i+vnum, -t27);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v1_i+2*vnum, t161);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v2_i+vnum, t162);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,ep_b_v2_i+2*vnum, -t103);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,v1_i+vnum, t170);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,v1_i+2*vnum, -lambda*t178);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,v2_i+vnum, -lambda*t173);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,v2_i+2*vnum, t186);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,w1_i+vnum, t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,w1_i+2*vnum, -lambda*t4*t128);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,w2_i+vnum, -t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i,w2_i+2*vnum, t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v1_i, t27);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v1_i+2*vnum, -t57);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v2_i, -t162);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,ep_b_v2_i+2*vnum, t166);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,v1_i, -t170);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,v1_i+2*vnum, t182);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,v2_i, t219);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,v2_i+2*vnum, -lambda*t189);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,w1_i, -t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,w1_i+2*vnum, t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,w2_i, t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+vnum,w2_i+2*vnum, -t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v1_i, -t161);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v1_i+vnum, t57);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v2_i, t103);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,ep_b_v2_i+vnum, -t166);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,v1_i, t194);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,v1_i+vnum, -t182);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,v2_i, -t186);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,v2_i+vnum, t228);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,w1_i, t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,w1_i+vnum, -t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,w2_i, -t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(ep_f_v2_i+2*vnum,w2_i+vnum, t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_b_v1_i+vnum, -t33);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_b_v1_i+2*vnum, t191);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_b_v2_i+vnum, t192);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_b_v2_i+2*vnum, -t110);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_f_v1_i+vnum, t193);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_f_v1_i+2*vnum, -t149);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_f_v2_i+vnum, -t170);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,ep_f_v2_i+2*vnum, t194);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,v2_i+vnum, t203);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,v2_i+2*vnum, -lambda*(t196+t208-ep_0_t*t128-t7*t14));
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,w1_i+vnum, -lambda*(t195-ep_0_t*t125));
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,w1_i+2*vnum, t198);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,w2_i+vnum, lambda*(t195-ep_0_t*t125));
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i,w2_i+2*vnum, -t198);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_b_v1_i, t33);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_b_v1_i+2*vnum, -t64);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_b_v2_i, -t192);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_b_v2_i+2*vnum, t199);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_f_v1_i, -t193);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_f_v1_i+2*vnum, t200);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_f_v2_i, t170);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,ep_f_v2_i+2*vnum, -t182);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,v2_i, -t203);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,v2_i+2*vnum, t212);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,w1_i, t207);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,w1_i+2*vnum, -lambda*(t205-ep_0_t*t142));
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,w2_i, -t207);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+vnum,w2_i+2*vnum, lambda*(t205-ep_0_t*t142));
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_b_v1_i, -t191);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_b_v1_i+vnum, t64);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_b_v2_i, t110);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_b_v2_i+vnum, -t199);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_f_v1_i, t149);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_f_v1_i+vnum, -t200);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_f_v2_i, -t194);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,ep_f_v2_i+vnum, t182);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,v2_i, t222);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,v2_i+vnum, -t212);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,w1_i, -t198);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,w1_i+vnum, t215);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,w2_i, t233);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v1_i+2*vnum,w2_i+vnum, -t215);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_b_v1_i+vnum, t216);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_b_v1_i+2*vnum, -t68);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_b_v2_i+vnum, -t90);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_b_v2_i+2*vnum, t217);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_f_v1_i+vnum, -t140);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_f_v1_i+2*vnum, t218);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_f_v2_i+vnum, t219);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,ep_f_v2_i+2*vnum, -t186);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,v1_i+vnum, -t203);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,v1_i+2*vnum, t222);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,w1_i+vnum, -lambda*t224);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,w1_i+2*vnum, t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,w2_i+vnum, t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i,w2_i+2*vnum, -t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_b_v1_i, -t216);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_b_v1_i+2*vnum, t227);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_b_v2_i, t90);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_b_v2_i+2*vnum, -t120);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_f_v1_i, t140);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_f_v1_i+2*vnum, -t159);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_f_v2_i, -t219);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,ep_f_v2_i+2*vnum, t228);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,v1_i, t203);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,v1_i+2*vnum, -t212);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,w1_i, t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,w1_i+2*vnum, -lambda*t231);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,w2_i, -t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+vnum,w2_i+2*vnum, t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_b_v1_i, t68);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_b_v1_i+vnum, -t227);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_b_v2_i, -t217);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_b_v2_i+vnum, t120);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_f_v1_i, -t218);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_f_v1_i+vnum, t159);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_f_v2_i, t186);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,ep_f_v2_i+vnum, -t228);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,v1_i, -t222);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,v1_i+vnum, t212);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,w1_i, -t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,w1_i+vnum, t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,w2_i, t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(v2_i+2*vnum,w2_i+vnum, -t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_b_v1_i+vnum, -t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_b_v1_i+2*vnum, t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_b_v2_i+vnum, t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_b_v2_i+2*vnum, -t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_f_v1_i+vnum, t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_f_v1_i+2*vnum, -t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_f_v2_i+vnum, -t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,ep_f_v2_i+2*vnum, t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,v1_i+vnum, t234);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,v1_i+2*vnum, -t198);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,v2_i+vnum, t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i,v2_i+2*vnum, -t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_b_v1_i, t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_b_v1_i+2*vnum, -t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_b_v2_i, -t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_b_v2_i+2*vnum, t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_f_v1_i, -t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_f_v1_i+2*vnum, t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_f_v2_i, t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,ep_f_v2_i+2*vnum, -t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,v1_i, -t207);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,v1_i+2*vnum, t235);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,v2_i, -t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+vnum,v2_i+2*vnum, t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_b_v1_i, -t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_b_v1_i+vnum, t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_b_v2_i, t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_b_v2_i+vnum, -t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_f_v1_i, t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_f_v1_i+vnum, -t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_f_v2_i, -t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,ep_f_v2_i+vnum, t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,v1_i, t233);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,v1_i+vnum, -t215);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,v2_i, t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w1_i+2*vnum,v2_i+vnum, -t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_b_v1_i+vnum, t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_b_v1_i+2*vnum, -t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_b_v2_i+vnum, -t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_b_v2_i+2*vnum, t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_f_v1_i+vnum, -t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_f_v1_i+2*vnum, t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_f_v2_i+vnum, t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,ep_f_v2_i+2*vnum, -t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,v1_i+vnum, -t207);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,v1_i+2*vnum, t233);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,v2_i+vnum, -t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i,v2_i+2*vnum, t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_b_v1_i, -t20);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_b_v1_i+2*vnum, t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_b_v2_i, t92);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_b_v2_i+2*vnum, -t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_f_v1_i, t144);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_f_v1_i+2*vnum, -t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_f_v2_i, -t165);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,ep_f_v2_i+2*vnum, t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,v1_i, t234);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,v1_i+2*vnum, -t215);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,v2_i, t230);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+vnum,v2_i+2*vnum, -t232);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_b_v1_i, t72);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_b_v1_i+vnum, -t43);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_b_v2_i, -t76);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_b_v2_i+vnum, t121);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_f_v1_i, -t132);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_f_v1_i+vnum, t160);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_f_v2_i, t190);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,ep_f_v2_i+vnum, -t175);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,v1_i, -t198);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,v1_i+vnum, t235);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,v2_i, -t226);
			IJV[ijv_idx++] = Eigen::Triplet<double>(w2_i+2*vnum,v2_i+vnum, t232);

			const_cnt++;
		}
	}
	if (ijv_idx!= lambda_hessian_IJV.size()) {
		cout << "error, ijv_idx " << ijv_idx << " but lambda_hessian_IJV.size() = " << lambda_hessian_IJV.size() << endl;
		exit(1);		
	}
};