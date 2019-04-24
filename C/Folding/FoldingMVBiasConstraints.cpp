#include "FoldingMVBiasConstraints.h"

using namespace std;

FoldingMVBiasConstraints::FoldingMVBiasConstraints(const Dog& dog, int curve_i): dog(dog), eS(dog.getEdgeStitching()), curve_i(curve_i) {
	int inner_v_n = eS.stitched_curves[curve_i].size()-2;
	const_n = inner_v_n;
	IJV.resize(24*const_n);
}

Eigen::VectorXd FoldingMVBiasConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0;
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

		double l1 = dog.stitched_curves_l[curve_i][edge_idx-1]; double l2 = dog.stitched_curves_l[curve_i][edge_idx];

		double t2 = ep_0_t-1.0;
		double t3 = t2*v2_x;
		double t4 = 1.0/l1;
		double t5 = 1.0/l2;
		double t6 = ep_b_t-1.0;
		double t7 = ep_f_t-1.0;
		double t8 = t2*v2_y;
		double t9 = ep_b_t*ep_b_v1_x;
		double t12 = ep_0_t*v1_x;
		double t10 = t3+t9-t12-ep_b_v2_x*t6;
		double t11 = ep_f_t*ep_f_v1_x;
		double t13 = l2*t10;
		double t14 = w1_x-w2_x;
		double t15 = t2*v2_z;
		double t16 = w1_z-w2_z;
		double t17 = ep_b_t*ep_b_v1_y;
		double t20 = ep_0_t*v1_y;
		double t18 = t8+t17-t20-ep_b_v2_y*t6;
		double t19 = ep_f_t*ep_f_v1_y;
		double t21 = l2*t18;
		double t22 = w1_y-w2_y;
		double t23 = ep_b_t*ep_b_v1_z;
		double t26 = ep_0_t*v1_z;
		double t24 = t15+t23-t26-ep_b_v2_z*t6;
		double t25 = ep_f_t*ep_f_v1_z;
		double t27 = t8+t19-t20-ep_f_v2_y*t7;
		double t28 = t21-l1*t27;
		double t29 = t3+t11-t12-ep_f_v2_x*t7;
		double t30 = t13-l1*t29;
		double t31 = t15+t25-t26-ep_f_v2_z*t7;
		double t32 = l2*t24-l1*t31;
		constVals(const_cnt++) = tanh((v1_y-v2_y)*(t4*t5*t14*t32-t4*t5*t16*t30)*1.0E3+(v1_x-v2_x)*(t4*t5*t16*t28-t4*t5*t22*t32)*1.0E3-(v1_z-v2_z)*(t4*t5*t14*t28-t4*t5*t22*t30)*1.0E3)*(-1.0/2.0)+1.0/2.0;
		//std::cout << "constVals(const_cnt--) = " << constVals(const_cnt-1) << std::endl;
	}
  if (const_cnt != const_n) {
		cout << "error in Vals, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  //std::cout << "constVals.norm() = " << constVals.norm() << std::endl;
  return constVals;
}

void FoldingMVBiasConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {
	// Add curve fold constraints
	int vnum = x.rows()/3;
	int const_cnt = 0; int ijv_cnt = 0;
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

		double l1 = dog.stitched_curves_l[curve_i][edge_idx-1]; double l2 = dog.stitched_curves_l[curve_i][edge_idx];

		double t2 = 1.0/l1;
		double t3 = w1_y-w2_y;
		double t4 = ep_0_t-1.0;
		double t5 = t4*v2_x;
		double t6 = 1.0/l2;
		double t7 = ep_b_t-1.0;
		double t8 = ep_f_t-1.0;
		double t9 = t4*v2_y;
		double t10 = v1_z-v2_z;
		double t11 = w1_z-w2_z;
		double t12 = ep_b_t*ep_b_v1_x;
		double t15 = ep_0_t*v1_x;
		double t31 = ep_b_v2_x*t7;
		double t13 = t5+t12-t15-t31;
		double t14 = ep_f_t*ep_f_v1_x;
		double t16 = l2*t13;
		double t17 = w1_x-w2_x;
		double t18 = t4*v2_z;
		double t19 = v1_y-v2_y;
		double t20 = ep_b_t*ep_b_v1_y;
		double t23 = ep_0_t*v1_y;
		double t36 = ep_b_v2_y*t7;
		double t21 = t9+t20-t23-t36;
		double t22 = ep_f_t*ep_f_v1_y;
		double t24 = l2*t21;
		double t25 = ep_b_t*ep_b_v1_z;
		double t28 = ep_0_t*v1_z;
		double t42 = ep_b_v2_z*t7;
		double t26 = t18+t25-t28-t42;
		double t27 = ep_f_t*ep_f_v1_z;
		double t30 = v1_x-v2_x;
		double t32 = ep_f_v2_x*t8;
		double t33 = t5+t14-t15-t32;
		double t34 = l1*t33;
		double t35 = t16-t34;
		double t37 = ep_f_v2_y*t8;
		double t38 = t9+t22-t23-t37;
		double t39 = l1*t38;
		double t40 = t24-t39;
		double t41 = t2*t6*t11*t35;
		double t43 = l2*t26;
		double t44 = ep_f_v2_z*t8;
		double t45 = t18+t27-t28-t44;
		double t46 = l1*t45;
		double t47 = t2*t6*t11*t40;
		double t29 = tanh(t19*(t41-t2*t6*t17*(t43-l1*(t18+t27-t44-ep_0_t*v1_z)))*-1.0E3+t30*(t47+t2*t3*t6*(t46-l2*t26))*1.0E3+t10*(t2*t3*t6*(t16-l1*(t5+t14-t32-ep_0_t*v1_x))-t2*t6*t17*(t24-l1*(t9+t22-t37-ep_0_t*v1_y)))*1.0E3);
		double t48 = t43-t46;
		double t50 = t2*t3*t6*t35;
		double t51 = t2*t6*t17*t40;
		double t52 = t50-t51;
		double t53 = t10*t52*1.0E3;
		double t54 = t2*t6*t17*t48;
		double t55 = t41-t54;
		double t56 = t19*t55*1.0E3;
		double t57 = t2*t3*t6*t48;
		double t58 = t47-t57;
		double t59 = t30*t58*1.0E3;
		double t60 = t53-t56+t59;
		double t49 = tanh(t60);
		double t61 = t49*t49;
		double t62 = t61-1.0;
		double t63 = ep_0_t*l1;
		double t65 = ep_0_t*l2;
		double t64 = t63-t65;
		double t66 = t2*t3*t6*t48*1.0E3;
		double t67 = l1*t4;
		double t70 = l2*t4;
		double t68 = t67-t70;
		double t69 = t2*t6*t11*t35*1.0E3;
		double t71 = t2*t3*t6*t35*1.0E3;
		double t72 = t2*t6*t10*t40*1.0E3;
		double t73 = t72-t2*t6*t19*t48*1.0E3;
		double t74 = t2*t6*t10*t35*1.0E3;
		double t75 = t74-t2*t6*t30*t48*1.0E3;
		double t76 = t62*t75*(1.0/2.0);
		double t77 = t2*t6*t19*t35*1.0E3;
		double t78 = t77-t2*t6*t30*t40*1.0E3;

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i,(t29*t29-1.0)*(ep_b_t*t2*t3*t10*1.0E3-ep_b_t*t2*t11*t19*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+vnum,t62*(ep_b_t*t2*t10*t17*1.0E3-ep_b_t*t2*t11*t30*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v1_i+2*vnum,t62*(ep_b_t*t2*t3*t30*1.0E3-ep_b_t*t2*t17*t19*1.0E3)*(-1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i,t62*(t2*t3*t7*t10*1.0E3-t2*t7*t11*t19*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+vnum,t62*(t2*t7*t10*t17*1.0E3-t2*t7*t11*t30*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_b_v2_i+2*vnum,t62*(t2*t3*t7*t30*1.0E3-t2*t7*t17*t19*1.0E3)*(1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i,t62*(ep_f_t*t3*t6*t10*1.0E3-ep_f_t*t6*t11*t19*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+vnum,t62*(ep_f_t*t6*t10*t17*1.0E3-ep_f_t*t6*t11*t30*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v1_i+2*vnum,t62*(ep_f_t*t3*t6*t30*1.0E3-ep_f_t*t6*t17*t19*1.0E3)*(1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i,t62*(t3*t6*t8*t10*1.0E3-t6*t8*t11*t19*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+vnum,t62*(t6*t8*t10*t17*1.0E3-t6*t8*t11*t30*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,ep_f_v2_i+2*vnum,t62*(t3*t6*t8*t30*1.0E3-t6*t8*t17*t19*1.0E3)*(-1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i,t62*(t66-t2*t6*t11*t40*1.0E3-t2*t3*t6*t10*t64*1.0E3+t2*t6*t11*t19*t64*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+vnum,t62*(t69-t2*t6*t17*t48*1.0E3+t2*t6*t10*t17*t64*1.0E3-t2*t6*t11*t30*t64*1.0E3)*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v1_i+2*vnum,t62*(t71-t2*t6*t17*t40*1.0E3-t2*t3*t6*t30*t64*1.0E3+t2*t6*t17*t19*t64*1.0E3)*(1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i,t62*(t66-t2*t6*t11*t40*1.0E3-t2*t3*t6*t10*t68*1.0E3+t2*t6*t11*t19*t68*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+vnum,t62*(t69-t2*t6*t17*t48*1.0E3+t2*t6*t10*t17*t68*1.0E3-t2*t6*t11*t30*t68*1.0E3)*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,v2_i+2*vnum,t62*(t71-t2*t6*t17*t40*1.0E3-t2*t3*t6*t30*t68*1.0E3+t2*t6*t17*t19*t68*1.0E3)*(-1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i,t62*t73*(-1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+vnum,t76);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w1_i+2*vnum,t62*t78*(-1.0/2.0));

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i,t62*t73*(1.0/2.0));
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+vnum,-t76);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,w2_i+2*vnum,t62*t78*(1.0/2.0));

		const_cnt++;
	}
  if (const_cnt != const_n) {
		cout << "error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}

void FoldingMVBiasConstraints::updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda_v) {
	// empty on purpose
};