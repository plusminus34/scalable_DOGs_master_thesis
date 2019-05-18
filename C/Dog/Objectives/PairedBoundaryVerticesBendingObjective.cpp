#include "PairedBoundaryVerticesBendingObjective.h"

PairedBoundaryVerticesBendingObjective::PairedBoundaryVerticesBendingObjective(const QuadTopology& quadTop, 
					const std::vector<std::pair<int,int>>& pairs, const Eigen::VectorXd& x, const Eigen::Vector3d& axis_direction) : 
		quadTop(quadTop), vnum(quadTop.v_n) {

	// First locate the indices of the pairs. We have 3 vertices for each pair (the paired vertices are "the same" and then we need their neighbours along a curve)
	obj_vertices.resize(3*pairs.size()); int obj_v_cnt = 0;
	for (auto p: pairs) {
		int v1 = p.first, v2 = p.second;
		obj_vertices(obj_v_cnt++) = v1;

		// now locate the neihbours of v1,v2 that is on the given direction, i.e. if v1 and v2 have two y neighbours and 1 x neighbours, we want the x neighbours
		//	so that we will add a bending objective on that curve that will smooth everything
		// This is given as input now as there are also corner boundary vertices
		int nb1 = find_neighbour_in_axis_direction(v1,quadTop,x,axis_direction);
		int nb2 = find_neighbour_in_axis_direction(v2,quadTop,x,axis_direction);
		obj_vertices(obj_v_cnt++) = nb1;
		obj_vertices(obj_v_cnt++) = nb2;
	}
	//std::cout << "obj_vertices.size() == "<< obj_vertices.size() << std::endl;
	//std::cout << "arrived here" << std::endl; int wait; std::cin >> wait;

	// Number of hessian triplets
	IJV.resize(27*obj_vertices.rows()/3);
	init_edge_lengths.resize(2*obj_vertices.rows()/3); 
	int cnt = 0;
	for (int si = 0; si < obj_vertices.rows(); si+=3) {
	    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

	    double ex_f_l = sqrt( pow(p0_x-pxf_x,2) + pow(p0_y-pxf_y,2) + pow(p0_z-pxf_z,2) );
		double ex_b_l = sqrt( pow(p0_x-pxb_x,2) + pow(p0_y-pxb_y,2) + pow(p0_z-pxb_z,2) );

		init_edge_lengths[cnt++] = ex_f_l;
		init_edge_lengths[cnt++] = ex_b_l;
	}
}

int PairedBoundaryVerticesBendingObjective::find_neighbour_in_axis_direction(int v, const QuadTopology& quadTop, const Eigen::VectorXd& x, \
	const Eigen::Vector3d& axis_direction) {
	double eps = 1e-10;
	for (auto nb_i : quadTop.A[v]) {
		// S bit ugly but will do the job
		// x direction so y difference is 0
		if (axis_direction(1) == 0) {
			if (abs(x(v+vnum)-x(nb_i+vnum)) < eps) return nb_i;
		} else { // y direction (x difference is 0)
			if (abs(x(v)-x(nb_i)) < eps) return nb_i;
		}
	}
	// Getting here means error
	std::cout << "Error at PairedBoundaryVerticesBendingObjective: Could not find a neighbour in the direction of the axis" << std::endl; exit(1);
}

double PairedBoundaryVerticesBendingObjective::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  
  int cnt = 0;  
  
  for (int si = 0; si < obj_vertices.rows(); si+=3) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

    double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];
    
	double t3 = 1.0/len_ex_b;
	double t4 = 1.0/len_ex_f;
	double t2 = t3*(p0_x-pxb_x)*2.0+t4*(p0_x-pxf_x)*2.0;
	double t5 = t3*(p0_y-pxb_y)*2.0+t4*(p0_y-pxf_y)*2.0;
	double t6 = t3*(p0_z-pxb_z)*2.0+t4*(p0_z-pxf_z)*2.0;
	e += t2*t2+t5*t5+t6*t6;

  }
  //std::cout << "e = " << std::endl;
  return e;
}

Eigen::VectorXd PairedBoundaryVerticesBendingObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int v_num = vnum;

  int cnt = 0;
 
  for (int si = 0; si < obj_vertices.rows(); si+=3) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

	double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];

	double t2 = 1.0/len_ex_b;
	double t3 = 1.0/len_ex_f;
	double t4 = t2*2.0;
	double t5 = t3*2.0;
	double t6 = t4+t5;
	double t7 = p0_x-pxb_x;
	double t8 = t2*t7*2.0;
	double t9 = p0_x-pxf_x;
	double t10 = t3*t9*2.0;
	double t11 = t8+t10;
	double t12 = p0_y-pxb_y;
	double t13 = t2*t12*2.0;
	double t14 = p0_y-pxf_y;
	double t15 = t3*t14*2.0;
	double t16 = t13+t15;
	double t17 = p0_z-pxb_z;
	double t18 = t2*t17*2.0;
	double t19 = p0_z-pxf_z;
	double t20 = t3*t19*2.0;
	double t21 = t18+t20;

	grad(p_0_i) += t6*t11*2.0;
	grad(p_0_i+v_num) += t6*t16*2.0;
	grad(p_0_i+2*v_num) += t6*t21*2.0;
	grad(p_xb_i+0) += t2*t11*-4.0;
	grad(p_xb_i+v_num) += t2*t16*-4.0;
	grad(p_xb_i+2*v_num) += t2*t21*-4.0;
	grad(p_xf_i) += t3*t11*-4.0;
	grad(p_xf_i+v_num) += t3*t16*-4.0;
	grad(p_xf_i+2*v_num) += t3*t21*-4.0;

  }
  // TODO: maybe add corners (or 4 vertices boundaries for cuts)
  //cout << "new grad.norm() = " << grad.norm() << endl; //exit(1);
  return grad;
}


void PairedBoundaryVerticesBendingObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  // Number of ijv values
  int v_num = vnum;
  int ijv_idx = 0; int cnt = 0;
 
   for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
 		const int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];

		double t3 = 1.0/len_ex_b;
		double t4 = t3*2.0;
		double t5 = 1.0/len_ex_f;
		double t6 = t5*2.0;
		double t2 = t4+t6;
		double t7 = t2*t2;
		double t8 = t7*2.0;
		double t9 = 1.0/(len_ex_b*len_ex_b);
		double t10 = t9*8.0;
		double t11 = t3*t5*8.0;
		double t12 = 1.0/(len_ex_f*len_ex_f);
		double t13 = t12*8.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t13);
	}
}