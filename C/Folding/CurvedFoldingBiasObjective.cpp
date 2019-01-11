#include "CurvedFoldingBiasObjective.h"

double CurvedFoldingBiasObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	int vnum = x.rows()/3;

	int h_cnt = 0;
	//for (int fold_i = 0; fold_i < curvedFoldBiases.size; fold_i++) {
	for (auto curvedFold : curvedFoldBiases) {
		int ep_0_v1_i(curvedFold.ep_0.edge.v1), ep_0_v2_i(curvedFold.ep_0.edge.v2); const double ep_0_t(curvedFold.ep_0.t);
		int ep_b_v1_i(curvedFold.ep_b.edge.v1), ep_b_v2_i(curvedFold.ep_b.edge.v2); const double ep_b_t(curvedFold.ep_b.t);
		int ep_f_v1_i(curvedFold.ep_f.edge.v1), ep_f_v2_i(curvedFold.ep_f.edge.v2); const double ep_f_t(curvedFold.ep_f.t);
		int v1_i(curvedFold.v1),v2_i(curvedFold.v2);


		const double ep_0_v1_x(x(ep_0_v1_i)); const double ep_0_v1_y(x(ep_0_v1_i+1*vnum)); const double ep_0_v1_z(x(ep_0_v1_i+2*vnum));
		/*
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		
  		*/
  		//e += t5*t5;
  }
  return e;
}

/*
ep_0_v1_x = sym('ep_0_v1_x', 'real');
ep_0_v1_y = sym('ep_0_v1_y', 'real');
ep_0_v1_z= sym('ep_0_v1_z', 'real');
ep_0_v1 = [ep_0_v1_x,ep_0_v1_y,ep_0_v1_z];

ep_b_v2_x = sym('ep_b_v2_x', 'real');
ep_b_v2_y = sym('ep_b_v2_y', 'real');
ep_b_v2_z= sym('ep_b_v2_z', 'real');
ep_b_v2 = [ep_b_v2_x,ep_b_v2_y,ep_b_v2_z];

ep_0_t = sym('ep_0_t', 'real'); assume(ep_0_t > 0);
ep_0 = ep_0_v1*ep_0_t+(1-ep_0_t)*ep_b_v2;

ep_b_v1_x = sym('ep_b_v1_x', 'real');
ep_b_v1_y = sym('ep_b_v1_y', 'real');
ep_b_v1_z= sym('ep_b_v1_z', 'real');
ep_b_v1 = [ep_b_v1_x,ep_b_v1_y,ep_b_v1_z];

ep_b_v2_x = sym('ep_b_v2_x', 'real');
ep_b_v2_y = sym('ep_b_v2_y', 'real');
ep_b_v2_z= sym('ep_b_v2_z', 'real');
ep_b_v2 = [ep_b_v2_x,ep_b_v2_y,ep_b_v2_z];

ep_b_t = sym('ep_b_t', 'real'); assume(ep_b_t > 0);

ep_b = ep_b_v1*ep_b_t+(1-ep_b_t)*ep_b_v2;

% Folded points
v1_x = sym('v1_x', 'real');
v1_y = sym('v1_y', 'real');
v1_z = sym('v1_z', 'real');
v1 = [v1_x,v1_y,v1_z];

v2_x = sym('v2_x', 'real');
v2_y = sym('v2_y', 'real');
v2_z = sym('v2_z', 'real');
v2 = [v2_x,v2_y,v2_z];

ep_f_v1_x = sym('ep_f_v1_x', 'real');
ep_f_v1_y = sym('ep_f_v1_y', 'real');
ep_f_v1_z= sym('ep_f_v1_z', 'real');
ep_f_v1 = [ep_f_v1_x,ep_f_v1_y,ep_f_v1_z];

ep_f_v2_x = sym('ep_f_v2_x', 'real');
ep_f_v2_y = sym('ep_f_v2_y', 'real');
ep_f_v2_z= sym('ep_f_v2_z', 'real');
*/