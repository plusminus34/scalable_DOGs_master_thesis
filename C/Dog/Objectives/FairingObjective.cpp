#include "FairingObjective.h"

FairingObjective::FairingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0) {
	// new strategy: Go over all stars
	// For all stars check each one of the 4 neighbours, and find the correct neighbour of that (check up to epsilon)
	// Do the same for bnd3, just with less neighbours
	// So place this in an inner method, and then just run it on stars and bnd3 vertices

	int vnum = x0.rows()/3;
	int const_n = 0;
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
		// x-forward curve
		add_indices_for_vertex(p_0_i, p_xb_i, p_xf_i, quadTop, x0, const_n);
		// x-backward curve
		add_indices_for_vertex(p_0_i, p_xf_i, p_xb_i, quadTop, x0, const_n);
		// y-forward curve
		add_indices_for_vertex(p_0_i, p_yb_i, p_yf_i, quadTop, x0, const_n);
		// y-backward curve
		add_indices_for_vertex(p_0_i, p_yf_i, p_yb_i, quadTop, x0, const_n);
		std::cout << std::endl << std::endl;
	}
	for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	   // x-forward curve
		add_indices_for_vertex(p_0_i, p_xb_i, p_xf_i, quadTop, x0, const_n);
		// x-backward curve
		add_indices_for_vertex(p_0_i, p_xf_i, p_xb_i, quadTop, x0, const_n);
	}
	set_ref(x0);
	// Then resize IJV
	IJV.resize(const_n*48);
	//IJV.resize(const_n*27);
}

void FairingObjective::add_indices_for_vertex(int p_0_i, int p_xb_i, int p_xf_i, 
				const QuadTopology& quadTop, const Eigen::VectorXd& x0, int& const_n) {
	int vnum = x0.rows()/3;
	Eigen::RowVector3d direction;
	direction << x0(p_0_i)-x0(p_xb_i),x0(p_0_i+vnum)-x0(p_xb_i+vnum),x0(p_0_i+2*+vnum)-x0(p_xb_i+2*+vnum); direction.normalize();
	std::cout << "direction = " << direction << std::endl;
	double eps = 1e-12;
	for (int j = 0; j < quadTop.A[p_xf_i].size(); j++) {
		int p_xff_i = quadTop.A[p_xf_i][j];
		//std::cout << "p_xb_i = " << p_xb_i << " p_xf_i = " << p_xf_i << " p_xff_i = " << p_xff_i << std::endl;
		Eigen::RowVector3d direction2; direction2 << x0(p_xff_i)-x0(p_xf_i),x0(p_xff_i+vnum)-x0(p_xf_i+vnum),x0(p_xff_i+2*+vnum)-x0(p_xf_i+2*+vnum); 
		direction2.normalize();
		//std::cout << "direction2 = " << direction2 << std::endl;
		if ((direction-direction2).norm() < eps) {
			p_xb.push_back(p_xb_i); p_0.push_back(p_0_i); p_xf.push_back(p_xf_i); p_xff.push_back(p_xff_i);
			const_n++;
			//std::cout << "found it" << std::endl; int wait; std::cin >> wait;
			return;
		}
	}
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

		
		double t3 = 1.0/len_ex_b;
		double t4 = 1.0/len_ex_f;
		double t5 = 1.0/len_ex_ff;
		double t2 = t3*(p0_x-pxb_x)*2.0+t4*(p0_x-pxf_x)*4.0-t5*(pxf_x-pxff_x)*2.0;
		double t6 = t3*(p0_y-pxb_y)*2.0+t4*(p0_y-pxf_y)*4.0-t5*(pxf_y-pxff_y)*2.0;
		double t7 = t3*(p0_z-pxb_z)*2.0+t4*(p0_z-pxf_z)*4.0-t5*(pxf_z-pxff_z)*2.0;
		e += t2*t2+t6*t6+t7*t7;
		
		/*
		double t2 = len_ex_b*p0_x+len_ex_f*p0_x-len_ex_f*pxb_x-len_ex_b*pxf_x;
		double t3 = 1.0/(len_ex_b*len_ex_b);
		double t4 = 1.0/(len_ex_f*len_ex_f);
		double t5 = len_ex_b*p0_y+len_ex_f*p0_y-len_ex_f*pxb_y-len_ex_b*pxf_y;
		double t6 = len_ex_b*p0_z+len_ex_f*p0_z-len_ex_f*pxb_z-len_ex_b*pxf_z;
		e += (t2*t2)*t3*t4*4.0+t3*t4*(t5*t5)*4.0+t3*t4*(t6*t6)*4.0;
		*/
	}
	std::cout << "FairingObjective::obj = " << e << std::endl;
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

		
		double t2 = 1.0/len_ex_b;
		double t3 = 1.0/len_ex_f;
		double t4 = t2*2.0;
		double t5 = t3*4.0;
		double t6 = t4+t5;
		double t7 = 1.0/len_ex_ff;
		double t8 = p0_x-pxb_x;
		double t9 = t2*t8*2.0;
		double t10 = p0_x-pxf_x;
		double t11 = t3*t10*4.0;
		double t12 = pxf_x-pxff_x;
		double t26 = t7*t12*2.0;
		double t13 = t9+t11-t26;
		double t14 = p0_y-pxb_y;
		double t15 = t2*t14*2.0;
		double t16 = p0_y-pxf_y;
		double t17 = t3*t16*4.0;
		double t18 = pxf_y-pxff_y;
		double t29 = t7*t18*2.0;
		double t19 = t15+t17-t29;
		double t20 = p0_z-pxb_z;
		double t21 = t2*t20*2.0;
		double t22 = p0_z-pxf_z;
		double t23 = t3*t22*4.0;
		double t24 = pxf_z-pxff_z;
		double t30 = t7*t24*2.0;
		double t25 = t21+t23-t30;
		double t27 = t7*2.0;
		double t28 = t5+t27;

		grad(p_0_i) += t6*t13*2.0;
		grad(p_0_i+v_num) += t6*t19*2.0;
		grad(p_0_i+2*v_num) += t6*t25*2.0;
		grad(p_xb_i) += t2*t13*-4.0;
		grad(p_xb_i+v_num) += t2*t19*-4.0;
		grad(p_xb_i+2*v_num) += t2*t25*-4.0;
		grad(p_xf_i) += t13*t28*-2.0;
		grad(p_xf_i+v_num) += t19*t28*-2.0;
		grad(p_xf_i+2*v_num) += t25*t28*-2.0;
		grad(p_xff_i) += t7*t13*4.0;
		grad(p_xff_i+v_num) += t7*t19*4.0;
		grad(p_xff_i+2*v_num) += t7*t25*4.0;
		
		/*
		double t2 = 1.0/(len_ex_b*len_ex_b);
		double t3 = 1.0/(len_ex_f*len_ex_f);
		double t4 = len_ex_b+len_ex_f;
		double t5 = len_ex_b*p0_x;
		double t6 = len_ex_f*p0_x;
		double t15 = len_ex_f*pxb_x;
		double t16 = len_ex_b*pxf_x;
		double t7 = t5+t6-t15-t16;
		double t8 = 1.0/len_ex_f;
		double t9 = len_ex_b*p0_y;
		double t10 = len_ex_f*p0_y;
		double t18 = len_ex_f*pxb_y;
		double t19 = len_ex_b*pxf_y;
		double t11 = t9+t10-t18-t19;
		double t12 = len_ex_b*p0_z;
		double t13 = len_ex_f*p0_z;
		double t20 = len_ex_f*pxb_z;
		double t21 = len_ex_b*pxf_z;
		double t14 = t12+t13-t20-t21;
		double t17 = 1.0/len_ex_b;

		grad(p_0_i) += t2*t3*t4*t7*8.0;
		grad(p_0_i+v_num) += t2*t3*t4*t11*8.0;
		grad(p_0_i+2*v_num) += t2*t3*t4*t14*8.0;
		grad(p_xb_i) += t2*t7*t8*-8.0;
		grad(p_xb_i+v_num) += t2*t8*t11*-8.0;
		grad(p_xb_i+2*v_num) += t2*t8*t14*-8.0;
		grad(p_xf_i) += t3*t7*t17*-8.0;
		grad(p_xf_i+v_num) += t3*t11*t17*-8.0;
		grad(p_xf_i+2*v_num) += t3*t14*t17*-8.0;
		*/
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
		
		double t3 = 1.0/len_ex_b;
		double t4 = t3*2.0;
		double t5 = 1.0/len_ex_f;
		double t6 = t5*4.0;
		double t2 = t4+t6;
		double t7 = 1.0/len_ex_ff;
		double t8 = t2*t2;
		double t9 = t8*2.0;
		double t10 = t7*2.0;
		double t11 = t6+t10;
		double t12 = t2*t7*4.0;
		double t13 = 1.0/(len_ex_b*len_ex_b);
		double t14 = t13*8.0;
		double t15 = t3*t11*4.0;
		double t16 = t11*t11;
		double t17 = t16*2.0;
		double t18 = 1.0/(len_ex_ff*len_ex_ff);
		double t19 = t18*8.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xff_i, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xff_i+vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xff_i+2*vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t14);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xff_i, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t14);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xff_i+vnum, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t14);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xff_i+2*vnum, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xff_i, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xff_i+vnum, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t2*t11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xff_i+2*vnum, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_0_i, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xb_i, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xf_i, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i,p_xff_i, t19);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_0_i+vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xb_i+vnum, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xf_i+vnum, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+vnum,p_xff_i+vnum, t19);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_0_i+2*vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xb_i+2*vnum, t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xf_i+2*vnum, t7*t11*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xff_i+2*vnum,p_xff_i+2*vnum, t19);
		/*
		double t2 = len_ex_b+len_ex_f;
		double t3 = 1.0/(len_ex_b*len_ex_b);
		double t4 = 1.0/(len_ex_f*len_ex_f);
		double t5 = t2*t2;
		double t6 = t3*t4*t5*8.0;
		double t7 = 1.0/len_ex_f;
		double t8 = 1.0/len_ex_b;
		double t9 = t3*8.0;
		double t10 = t7*t8*8.0;
		double t11 = t4*8.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t6);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, t2*t3*t7*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t2*t4*t8*-8.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t11);
		*/
  }
}