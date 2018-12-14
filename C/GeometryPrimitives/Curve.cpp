#include "Curve.h"

#include <vector>
#include "LinearTransformations.h"

//#include "Utils.h"

Curve::Curve(const std::vector<double>& _len, const std::vector<double>& _k, const std::vector<double>& _t) : len(_len), k(_k), t(_t) {
	// empty on purpose
}

Curve::Curve(const Eigen::MatrixXd& coords) {
	int vn = coords.rows();
	len.resize(vn-1); k.resize(vn-2); t.resize(vn-3);

	for (int i = 0; i < vn-1; i++) {
		len[i] = (coords.row(i+1)-coords.row(i)).norm();
	}
	auto cur_frame = get_frame(coords,1);
	int k_sign = 1;
	for (int i = 0; i < vn-2; i++) {
		Eigen::RowVector3d e1 = (coords.row(i+1)-coords.row(i)).normalized();
		Eigen::RowVector3d e2 = (coords.row(i+2)-coords.row(i+1)).normalized();
		double angle = get_angle_and_orientation(e1,e2); //cout << "angle = " << angle << endl;
		k[i] = k_sign*2*sin(angle)/(coords.row(i+2)-coords.row(i)).norm();

		if ( (i > 1) && (k[i]!=0) ) {
			auto new_frame = get_frame(coords,i+1);
			auto old_b = cur_frame.col(2); auto new_b = new_frame.col(2);
			//std::cout << "old_b = " << old_b << " new_b = " << new_b << std::endl;
			if (old_b.dot(new_b) < 0) {
				k[i] = -1*k[i]; 
				k_sign = -1*k_sign;
			}
			cur_frame = new_frame;
		}
		/*
		// maybe here I have the tangents, normals, binormal
		// Can check the equality of the next edge by t,n with angle
		Eigen::RowVector3d ef = (coords.row(i+2)-coords.row(i+1)).normalized(), eb = (coords.row(i+1)-coords.row(i)).normalized();
		Eigen::RowVector3d t1 = (ef+eb).normalized(), n1 = (ef-eb).normalized();
		
		cout << "t,n way diff norm = " << 
				(coords.row(i+2)-coords.row(i+1) - len[i+1]*(cos(0.5*angle)*t1+sin(0.5*angle)*n1)).norm() << endl;
		
		// Can also check the equality by rotation around the binormal
		Eigen::RowVector3d b = e1.cross(e2).normalized(); cout << "b = " << b << endl;
		Eigen::RowVector3d zero3d; zero3d.setZero();
		Eigen::RowVector3d e_f = rotate_vec(e1, zero3d, b, angle).normalized();
		cout << "binormal way diff norm = " << 
				(coords.row(i+2)-coords.row(i+1) - len[i+1]*e_f).norm() << endl;
		*/
	}
	auto old_frame = get_frame(coords,1);
	Eigen::RowVector3d old_b = get_frame(coords,1).col(2);
	for (int i = 0; i < vn-3; i++) {
		Eigen::RowVector3d e1 = (coords.row(i+1)-coords.row(i)).normalized();
		Eigen::RowVector3d e2 = (coords.row(i+2)-coords.row(i+1)).normalized();
		Eigen::RowVector3d e3 = (coords.row(i+3)-coords.row(i+2)).normalized();

		if (k[i+1] == 0){
			t[i] = 0;
		} else {
			auto new_frame = get_frame(coords,i+2);
			auto new_b = get_frame(coords,i+2).col(2);
			double angle = get_angle_and_orientation(old_b,new_b);
			// This quantity, as opposed to curvature is defined on the edge and so normalized by it's length
			//int sign = 1;
			//if (old_frame.col(2).dot(new_frame.col(0)) < 0) {sign = -1;}
			t[i] = sin(angle)/(coords.row(i+2)-coords.row(i+1)).norm();

			Eigen::RowVector3d zero3d; zero3d.setZero();
			auto check1 = rotate_vec(old_b, zero3d, e2.normalized(), angle);
			auto check2 = rotate_vec(old_b, zero3d, e2.normalized(), -angle);
			if ((check2-new_b.transpose()).norm() < (check1-new_b.transpose()).norm()) t[i] = -t[i];
			//std::cout << "check1 = " << check1 << " check2 = " << check2 << " new_b = " << new_b <<std::endl;
			//t[i] = 0;
			old_b = new_b;
			new_frame = old_frame;
			//exit(1);
		}
	}
}

Curve::Curve(const SurfaceCurve& surfaceCurve, const Eigen::MatrixXd& V) :Curve(surfaceCurve.get_curve_coords(V)) {
	// empty on purpose
}

Curve::Curve(const std::vector<double>& len1, const std::vector<double>& k1, const std::vector<double>& t1,
		const std::vector<double>& len2, const std::vector<double>& k2, const std::vector<double>& t2, double time) {
	for (int i = 0; i < len1.size(); i++) {
		len.push_back(time*len2[i]+(1-time)*len1[i]);
	}
	for (int i = 0; i < k1.size(); i++) {
		k.push_back(time*k2[i]+(1-time)*k1[i]);
	}
	for (int i = 0; i < t1.size(); i++) {
		t.push_back(time*t2[i]+(1-time)*t1[i]);
	}

}

Curve::Curve(const Curve& curve1, const Curve& curve2, double time) : Curve(curve1.len,curve1.k,curve1.t,curve2.len,curve2.k,curve2.t,time) {
	// empty on purpose
}

Eigen::MatrixXd Curve::getCoords(const Eigen::RowVector3d& T, const Eigen::Matrix3d& F) {
	int vn = len.size()+1; Eigen::MatrixXd coords(vn,3);
	// As written in the header, the translation and rotation are for the second point on the curve (need a point with an osculation plane)
	coords.row(1) = T;
	double alpha = get_angle_from_lengths_and_k(len[0],len[1],k[0]);
	Eigen::RowVector3d T1 = F.col(0), N1 = F.col(1); // tangent and principle normal of the frame
	// Build backwards and forwards vertices (no torsion information for first 3 vertices)
	coords.row(2) = coords.row(1) + len[1]*(cos(0.5*alpha)*T1+sin(0.5*alpha)*N1);
	coords.row(0) = coords.row(1) + len[0]*(-cos(0.5*alpha)*T1+sin(0.5*alpha)*N1);

	Eigen::RowVector3d old_b = F.col(2);

	Eigen::RowVector3d zero3d; zero3d.setZero();
	for (int i = 3; i < vn; i++) {
		Eigen::RowVector3d e_bb = coords.row(i-2)-coords.row(i-3), e_b = coords.row(i-1)-coords.row(i-2);

		// Todo: here, needs a way to deal with straight lines.. shouldn't be hard. We can first decide on some binormal and always use this one while flat
		/*double angle_eb_ebb = get_angle(e_bb,e_b); 
		if ((angle_eb_ebb > 1e-10) && (angle_eb_ebb < M_PI- 1e-10)) (get_angle(e_b)) {

		}*/
		
		double euc_l = e_b.norm();
		double edge_torsion = t[i-3];
		double t_alpha = asin(clip(edge_torsion*euc_l,-1,1));
		// Rotate the binormal by the angle given by the torsion and lengths
		Eigen::RowVector3d b = old_b;
		//if (std::abs(t_alpha) > 1e-10) {
		if (t_alpha) {
			b = rotate_vec(old_b, zero3d, e_b.normalized(), t_alpha);
		}
		//cout << "b = " <<  b << endl;
		// Rotate the edge in axis of the normal to the new osculating plane, by the angle given by the curvature
		double k_alpha = get_angle_from_lengths_and_k(len[i-2],len[i-1],k[i-2]); //cout << "k_alpha = " << k_alpha << endl;
		Eigen::RowVector3d e_f = rotate_vec(e_b, zero3d, b, k_alpha);
		//Eigen::RowVector3d e_f = e_b*cos(k_alpha)+ (b.cross(e_b))*sin(k_alpha)+b*(b.dot(e_b))*(1-cos(k_alpha));//, zero3d, b, k_alpha);
		coords.row(i) = coords.row(i-1)+len[i-1]*(e_f.normalized());
		old_b = b;
	}
	return coords;
}

Eigen::Matrix3d Curve::get_frame(const Eigen::MatrixXd& coords, int i) {
	Eigen::RowVector3d ef = (coords.row(i+1)-coords.row(i)).normalized();
	Eigen::RowVector3d eb = (coords.row(i-1)-coords.row(i)).normalized();

	Eigen::RowVector3d T = (ef-eb).normalized();
	Eigen::RowVector3d N = (ef+eb).normalized();
	Eigen::RowVector3d B = T.cross(N);
	Eigen::Matrix3d frame; frame.col(0) = T; frame.col(1) = N; frame.col(2) = B;
	return frame;
}

double Curve::get_angle_and_orientation(Eigen::RowVector3d e1,Eigen::RowVector3d e2) {
	auto vec_norm_prod = e1.norm()*e2.norm();
	auto dot_prod = e1.dot(e2)/vec_norm_prod;
	auto cos_angle = clip(dot_prod,-1,1);
	if (cos_angle == 1) return acos(cos_angle);

	/*
	auto sin_angle = clip(e1.cross(e2).norm()/vec_norm_prod,-1,1);
	std::cout << "cos_angle = " << cos_angle << std::endl; 
	std::cout << "sin_angle = " << sin_angle << std::endl;
	auto tan_angle = sin_angle/cos_angle;
	return atan(tan_angle);
	*/

	auto angle = acos(cos_angle);
	auto b = e1.cross(e2).normalized();
	Eigen::RowVector3d zero3d; zero3d.setZero();
	auto diff1 = (e2-rotate_vec(e1, zero3d, b, angle)).norm();
	auto diff2 = (e2-rotate_vec(e1, zero3d, b, -angle)).norm();
	if (diff2 < diff1 ) angle = -angle;
	//std::cout << "angle = " << angle << std::endl;
	//std::cout << "diff1 = " << diff1 << " diff2 = " << diff2 << std::endl;
	//std::cout << "b = " << b << std::endl;
	return angle;

	//auto angle = acos(cos_angle);
	//std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
	//std::cout << "dot_prod = " << dot_prod << std::endl;
	


	// = rotate_vec(e_b, zero3d, b, k_alpha);
	
	//return angle;
}

void Curve::getTranslationAndFrameFromCoords(const Eigen::MatrixXd& coords, Eigen::RowVector3d& T, Eigen::Matrix3d& F) {
	T = coords.row(1);
	Eigen::RowVector3d ef = (coords.row(2)-coords.row(1)).normalized(), eb = (coords.row(1)-coords.row(0)).normalized();
	Eigen::RowVector3d t1 = (ef+eb).normalized(), n1 = (ef-eb).normalized();
	F.col(0) = t1; F.col(1) = n1; F.col(2) = (t1.cross(n1)).normalized();
	//cout << "coords.row(1) = " << coords.row(1) << endl;
	//cout << "coords.row(2) = " << coords.row(2) << endl;
	double len = (coords.row(2)-coords.row(1)).norm();
	double alpha = acos(clip((coords.row(2)-coords.row(1)).normalized().dot((coords.row(1)-coords.row(0)).normalized()),-1,1));
	//cout << "alpha = " << alpha << endl;
	//cout << "coords.row(1) + len[1]*(cos(0.5*alpha)*T1+sin(0.5*alpha)*N1) = " << coords.row(1) + len*(cos(0.5*alpha)*t1+sin(0.5*alpha)*(n1)) << endl;
	//exit(1);
}

double Curve::get_angle_from_lengths_and_k(double l1, double l2, double k) {
	if (std::abs(k) < 1e-12) return 0; // Todo: make sure that's correct and not pi
	double s1 = euc_to_arc(l1,k);
	double s2 = euc_to_arc(l2,k);
	double s= s1+s2; // Get the full arc of the circle with curvature k
	// now get the euc length from arc length
	double euc = arc_to_euc(s,k);
	return asin(clip(k*euc/2,-1,1));
}

double Curve::euc_to_arc(double arc, double k) {
	// see http://mathworld.wolfram.com/CircularSegment.html
	double R = 1./k;
	double theta = 2*asin(clip(arc/(2*R),-1,1));
	return R*theta;
}

double Curve::arc_to_euc(double s, double k) {
	// see http://mathworld.wolfram.com/CircularSegment.html
	double R = 1./k;
	double theta = s/R;
	return 2*R*sin(0.5*theta);
}

void Curve::print_geometric_represenation() {
	std::cout << "curve with " << 1+len.size() << " points" << std::endl;
	std::cout << "lengths: "; for (auto l: len) std::cout << " " << l << ","; std::cout<<std::endl;
	std::cout << "curvatures: "; for (auto ki: k) std::cout << " " << ki << ","; std::cout<<std::endl;
	std::cout << "torsions: "; for (auto ti: t) std::cout << " " << ti << ","; std::cout<<std::endl;
	//std::vector<double> len, k, t;
}

double Curve::clip(double n, double lower, double upper) {
  return std::max(lower, std::min(n, upper));
}