#include "PointsRigidAlignmentObjective.h"
#include <igl/procrustes.h>

using namespace std;


PointsRigidAlignmentObjective::PointsRigidAlignmentObjective(std::vector<int>& src_points, std::vector<int>& target_points) :
		src_points(src_points), target_points(target_points) { 

	IJV.resize(src_points.size()*30);
}

void PointsRigidAlignmentObjective::update_rigid_motion(const Eigen::VectorXd& x, 
	const std::vector<int>& src_points, const std::vector<int>& target_points,
			Eigen::Matrix3d& R, Eigen::Vector3d& T) {
	if (src_points.size() < 3) return;
	int vnum = x.rows()/3;
	Eigen::MatrixXd src(src_points.size(),3); Eigen::MatrixXd target(target_points.size(),3);
	for (int i = 0; i < src_points.size(); i++) {
		src.row(i) << x(src_points[i]), x(src_points[i]+vnum), x(src_points[i]+2*vnum);
		target.row(i) << x(target_points[i]), x(target_points[i]+vnum), x(target_points[i]+2*vnum);
	}
	// call procrustes with includeScaling = false, includeReflections = false
	double scale_dummy; igl::procrustes(src, target, false, false, scale_dummy,R,T);
	//std::cout << "R = " << endl << R << endl << "T = " << endl << T << endl;
	//int wait; std::cin >> wait; exit(1);
}

double PointsRigidAlignmentObjective::obj(const Eigen::VectorXd& x) const {
	if (src_points.size() < 3) return 0;
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,src_points,target_points,R,T);
	double t1(T(0)),t2(T(1)),t3(T(2));
	double r11(R(0,0)),r12(R(0,1)),r13(R(0,2)),r21(R(1,0)),r22(R(1,1)),r23(R(1,2)),r31(R(2,0)),r32(R(2,1)),r33(R(2,2));

	double e = 0;
	int vnum = x.rows()/3;
	#pragma clang loop vectorize(enable)
	for (int i = 0; i < src_points.size(); i++) {
		int p_1_i = src_points[i], p_2_i = target_points[i];
		const double p1_x(x(p_1_i+0)); const double p1_y(x(p_1_i+1*vnum)); const double p1_z(x(p_1_i+2*vnum));
		const double p2_x(x(p_2_i+0)); const double p2_y(x(p_2_i+1*vnum)); const double p2_z(x(p_2_i+2*vnum));

		double t5 = -p2_x+t1+p1_x*r11+p1_y*r21+p1_z*r31;
 		double t6 = -p2_y+t2+p1_x*r12+p1_y*r22+p1_z*r32;
 		double t7 = -p2_z+t3+p1_x*r13+p1_y*r23+p1_z*r33;
		e += t5*t5+t6*t6+t7*t7;
	}
	std::cout << "e rigid = " << e << std::endl;
	return e;
}

Eigen::VectorXd PointsRigidAlignmentObjective::grad(const Eigen::VectorXd& x) const {
	Eigen::VectorXd grad;
  	grad.resize(x.rows(),1); grad.setZero();
	if (src_points.size() < 3) return grad;
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,src_points,target_points,R,T);
	double t1(T(0)),t2(T(1)),t3(T(2));
	double r11(R(0,0)),r12(R(0,1)),r13(R(0,2)),r21(R(1,0)),r22(R(1,1)),r23(R(1,2)),r31(R(2,0)),r32(R(2,1)),r33(R(2,2));

  	int vnum = x.rows()/3;

  	#pragma clang loop vectorize(enable)
  	for (int i = 0; i < src_points.size(); i++) {
		int p_1_i = src_points[i], p_2_i = target_points[i];
		const double p1_x(x(p_1_i+0)); const double p1_y(x(p_1_i+1*vnum)); const double p1_z(x(p_1_i+2*vnum));
		const double p2_x(x(p_2_i+0)); const double p2_y(x(p_2_i+1*vnum)); const double p2_z(x(p_2_i+2*vnum));

		double t5 = p1_x*r11;
		double t6 = p1_y*r21;
		double t7 = p1_z*r31;
		double t8 = -p2_x+t1+t5+t6+t7;
		double t9 = p1_x*r12;
		double t10 = p1_y*r22;
		double t11 = p1_z*r32;
		double t12 = -p2_y+t2+t9+t10+t11;
		double t13 = p1_x*r13;
		double t14 = p1_y*r23;
		double t15 = p1_z*r33;
		double t16 = -p2_z+t3+t13+t14+t15;

		grad(p_1_i) += r11*t8*2.0+r12*t12*2.0+r13*t16*2.0;
		grad(p_1_i + vnum) += r21*t8*2.0+r22*t12*2.0+r23*t16*2.0;
		grad(p_1_i + 2*vnum) += r31*t8*2.0+r32*t12*2.0+r33*t16*2.0;
		//grad(p_2_i) += p2_x*2.0-t1*2.0-p1_x*r11*2.0-p1_y*r21*2.0-p1_z*r31*2.0;
		//grad(p_2_i + vnum) += p2_y*2.0-t2*2.0-p1_x*r12*2.0-p1_y*r22*2.0-p1_z*r32*2.0;
		//grad(p_2_i + 2*vnum) += p2_z*2.0-t3*2.0-p1_x*r13*2.0-p1_y*r23*2.0-p1_z*r33*2.0;
	}
  	return grad;
}


void PointsRigidAlignmentObjective::updateHessianIJV(const Eigen::VectorXd& x) {
	if (src_points.size() < 3) return;
	Eigen::Matrix3d R; Eigen::Vector3d T;
	update_rigid_motion(x,src_points,target_points,R,T);
	double t1(T(0)),t2(T(1)),t3(T(2));
	double r11(R(0,0)),r12(R(0,1)),r13(R(0,2)),r21(R(1,0)),r22(R(1,1)),r23(R(1,2)),r31(R(2,0)),r32(R(2,1)),r33(R(2,2));

	int vnum = x.rows()/3;
	int ijv_idx = 0;
  	#pragma clang loop vectorize(enable)
  	for (int i = 0; i < src_points.size(); i++) {
		int p_1_i = src_points[i], p_2_i = target_points[i];
		const double p1_x(x(p_1_i+0)); const double p1_y(x(p_1_i+1*vnum)); const double p1_z(x(p_1_i+2*vnum));
		const double p2_x(x(p_2_i+0)); const double p2_y(x(p_2_i+1*vnum)); const double p2_z(x(p_2_i+2*vnum));

		double t2 = r11*r21*2.0;
		double t3 = r12*r22*2.0;
		double t4 = r13*r23*2.0;
		double t5 = t2+t3+t4;
		double t6 = r11*r31*2.0;
		double t7 = r12*r32*2.0;
		double t8 = r13*r33*2.0;
		double t9 = t6+t7+t8;
		double t10 = r21*r31*2.0;
		double t11 = r22*r32*2.0;
		double t12 = r23*r33*2.0;
		double t13 = t10+t11+t12;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_1_i, (r11*r11)*2.0+(r12*r12)*2.0+(r13*r13)*2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_1_i+vnum, t5);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_1_i+2*vnum, t9);
		/*
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_2_i, r11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_2_i+vnum, r12*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i,p_2_i+2*vnum, r13*-2.0);
		*/
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_1_i, t5);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_1_i+vnum, (r21*r21)*2.0+(r22*r22)*2.0+(r23*r23)*2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_1_i+2*vnum, t13);
		/*
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_2_i, r21*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_2_i+vnum, r22*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+vnum,p_2_i+2*vnum, r23*-2.0);
		*/
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_1_i, t9);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_1_i+vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_1_i+2*vnum, (r31*r31)*2.0+(r32*r32)*2.0+(r33*r33)*2.0);

		/*
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_2_i, r31*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_2_i+vnum, r32*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_1_i+2*vnum,p_2_i+2*vnum, r33*-2.0);

		
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i,p_1_i, r11*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i,p_1_i+vnum, r21*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i,p_1_i+2*vnum, r31*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i,p_2_i, 2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+vnum,p_1_i, r12*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+vnum,p_1_i+vnum, r22*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+vnum,p_1_i+2*vnum, r32*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+vnum,p_2_i+vnum, 2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+2*vnum,p_1_i, r13*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+2*vnum,p_1_i+vnum, r23*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+2*vnum,p_1_i+2*vnum, r33*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_2_i+2*vnum,p_2_i+2*vnum, 2.0);
		*/
	}
}
