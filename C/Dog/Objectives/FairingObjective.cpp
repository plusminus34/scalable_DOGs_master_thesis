#include "FairingObjective.h"

FairingObjective::FairingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0) {
	// set up indices for std::vector<int> p_xb,p_0,p_xf,p_xff;
	// for every star vertex, add nearby vertices if they are stars as well
	int const_n = 0;
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
		// x-forward curve
		if (quadTop.vi_to_star[p_xf_i] != -1) {
			int nb_star = quadTop.vi_to_star[p_xf_i];
			p_xb.push_back(p_xb_i); p_0.push_back(p_0_i); p_xf.push_back(p_xf_i); p_xff.push_back(quadTop.stars(nb_star+1));
			const_n++;
		}
		// x-backward curve
		if (quadTop.vi_to_star[p_xb_i] != -1) {
			int nb_star = quadTop.vi_to_star[p_xb_i];
			p_xb.push_back(quadTop.stars(nb_star+3)); p_0.push_back(p_xb_i); p_xf.push_back(p_0_i); p_xff.push_back(p_xf_i);
			const_n++;
		}
		// y-upper curve
		if (quadTop.vi_to_star[p_yf_i] != -1) {
			int nb_star = quadTop.vi_to_star[p_yf_i];
			p_xb.push_back(p_yb_i); p_0.push_back(p_0_i); p_xf.push_back(p_yf_i); p_xff.push_back(quadTop.stars(nb_star+2));
			const_n++;
		}
		// y - lower curve
		if (quadTop.vi_to_star[p_yb_i] != -1) {
			int nb_star = quadTop.vi_to_star[p_yb_i];
			p_xb.push_back(quadTop.stars(nb_star+4)); p_0.push_back(p_yb_i); p_xf.push_back(p_0_i); p_xff.push_back(p_yf_i);
			const_n++;
		}
		std::cout << std::endl << std::endl;
	}
	//int wait; std::cin >> wait;
	set_ref(x0);
	// Then resize IJV
	IJV.resize(const_n*48);
}

void FairingObjective::set_ref(const Eigen::VectorXd& x) {
	int vnum = x.rows()/3;
	for (int i = 0; i < p_xb.size(); i++) {
		int p_xb_i(p_xb[i]), p_0_i(p_0[i]), p_xf_i(p_xf[i]), p_xff_i(p_xff[i]);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pxff_x(x(p_xff_i+0)); const double pxff_y(x(p_xff_i+1*vnum)); const double pxff_z(x(p_xff_i+2*vnum));

		// todo save lengths
		len_ex_b_v.push_back(sqrt(pow(p0_x-pxb_x,2)+pow(p0_y-pxb_y,2)+pow(p0_z-pxb_z,2)));
		len_ex_f_v.push_back(sqrt(pow(pxf_x-p0_x,2)+pow(pxf_y-p0_y,2)+pow(pxf_z-p0_z,2)));
		len_ex_ff_v.push_back(sqrt(pow(pxff_x-pxf_x,2)+pow(pxff_y-pxf_y,2)+pow(pxff_z-pxf_z,2)));
	}
}

double FairingObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	int vnum = x.rows()/3;
	for (int i = 0; i < p_xb.size(); i++) {
		int p_xb_i(p_xb[i]), p_0_i(p_0[i]), p_xf_i(p_xf[i]), p_xff_i(p_xff[i]);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pxff_x(x(p_xff_i+0)); const double pxff_y(x(p_xff_i+1*vnum)); const double pxff_z(x(p_xff_i+2*vnum));

		const double len_ex_b(len_ex_b_v[i]); const double len_ex_f(len_ex_f_v[i]); const double len_ex_ff(len_ex_ff_v[i]);
		std::cout << "len_ex_b = " << len_ex_b << std::endl;
		std::cout << "len_ex_f = " << len_ex_f << std::endl;
		std::cout << "len_ex_ff = " << len_ex_ff << std::endl;

		double t2 = len_ex_ff*p0_x-len_ex_ff*pxb_x-len_ex_b*pxf_x+len_ex_b*pxff_x;
		double t3 = 1.0/(len_ex_b*len_ex_b);
		double t4 = 1.0/(len_ex_ff*len_ex_ff);
		double t5 = len_ex_ff*p0_y-len_ex_ff*pxb_y-len_ex_b*pxf_y+len_ex_b*pxff_y;
		double t6 = len_ex_ff*p0_z-len_ex_ff*pxb_z-len_ex_b*pxf_z+len_ex_b*pxff_z;
		e += (t2*t2)*t3*t4*4.0+t3*t4*(t5*t5)*4.0+t3*t4*(t6*t6)*4.0;
	}
	std::cout << "FairingObjective::obj = " << e << std::endl;
	int wait; std::cin >> wait;
	return e;
}

Eigen::VectorXd FairingObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int vnum = x.rows()/3;
  int v_num = vnum;
 
  #pragma clang loop vectorize(enable)
	for (int i = 0; i < p_xb.size(); i++) {
		int p_xb_i(p_xb[i]), p_0_i(p_0[i]), p_xf_i(p_xf[i]), p_xff_i(p_xff[i]);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pxff_x(x(p_xff_i+0)); const double pxff_y(x(p_xff_i+1*vnum)); const double pxff_z(x(p_xff_i+2*vnum));

		const double len_ex_b(len_ex_b_v[i]); const double len_ex_f(len_ex_f_v[i]); const double len_ex_ff(len_ex_ff_v[i]);

		double t2 = 1.0/(len_ex_b*len_ex_b);
		double t3 = 1.0/len_ex_ff;
		double t4 = len_ex_ff*p0_x;
		double t5 = len_ex_b*pxff_x;
		double t16 = len_ex_ff*pxb_x;
		double t17 = len_ex_b*pxf_x;
		double t6 = t4+t5-t16-t17;
		double t7 = t2*t3*t6*8.0;
		double t8 = len_ex_ff*p0_y;
		double t9 = len_ex_b*pxff_y;
		double t20 = len_ex_ff*pxb_y;
		double t21 = len_ex_b*pxf_y;
		double t10 = t8+t9-t20-t21;
		double t11 = t2*t3*t10*8.0;
		double t12 = len_ex_ff*p0_z;
		double t13 = len_ex_b*pxff_z;
		double t22 = len_ex_ff*pxb_z;
		double t23 = len_ex_b*pxf_z;
		double t14 = t12+t13-t22-t23;
		double t15 = t2*t3*t14*8.0;
		double t18 = 1.0/len_ex_b;
		double t19 = 1.0/(len_ex_ff*len_ex_ff);

		grad(p_0_i) += t7;
		grad(p_0_i+v_num) += t11;
		grad(p_0_i+2*v_num) += t15;
		grad(p_xb_i) += -t7;
		grad(p_xb_i+v_num) += -t11;
		grad(p_xb_i+2*v_num) += -t15;
		grad(p_xf_i) += t6*t18*t19*-8.0;
		grad(p_xf_i+v_num) += t10*t18*t19*-8.0;
		grad(p_xf_i+2*v_num) += t14*t18*t19*-8.0;
		grad(p_xff_i) += t6*t18*t19*8.0;
		grad(p_xff_i+v_num) += t10*t18*t19*8.0;
		grad(p_xff_i+2*v_num) += t14*t18*t19*8.0;
  }
  return grad;
}

void FairingObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  int vnum = x.rows()/3;
  int v_num = vnum;
  int h_cnt = 0;

  int ijv_idx = 0;
  #pragma clang loop vectorize(enable)
	for (int i = 0; i < p_xb.size(); i++) {
		int p_xb_i(p_xb[i]), p_0_i(p_0[i]), p_xf_i(p_xf[i]), p_xff_i(p_xff[i]);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pxff_x(x(p_xff_i+0)); const double pxff_y(x(p_xff_i+1*vnum)); const double pxff_z(x(p_xff_i+2*vnum));

		const double len_ex_b(len_ex_b_v[i]); const double len_ex_f(len_ex_f_v[i]); const double len_ex_ff(len_ex_ff_v[i]);

		double t2 = 1.0/(len_ex_b*len_ex_b);
		double t3 = t2*8.0;
		double t4 = 1.0/len_ex_b;
		double t5 = 1.0/len_ex_ff;
		double t6 = t4*t5*8.0;
		double t7 = 1.0/(len_ex_ff*len_ex_ff);
		double t8 = t7*8.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t4*t5*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xff_i, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xff_i+vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xff_i+2*vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xff_i, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xff_i+vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, -t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t3);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xff_i+2*vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xff_i, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xff_i+vnum, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xff_i+2*vnum, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_0_i, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xb_i, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xf_i, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xff_i, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_0_i+vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xb_i+vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xf_i+vnum, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xff_i+vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_0_i+2*vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xb_i+2*vnum, -t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xf_i+2*vnum, -t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xff_i+2*vnum, t8);
  }
}