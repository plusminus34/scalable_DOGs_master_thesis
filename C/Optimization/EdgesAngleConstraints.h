#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

class EdgesAngleConstraints : public Constraints {
public:
	EdgesAngleConstraints(const QuadTopology& quadTop, const Eigen::VectorXd& x, const std::vector<std::pair<Edge,Edge>>& edge_pairs, const std::vector<double> cos_angles) :
			vnum(quadTop.v_n), edge_pairs(edge_pairs), cos_angles(cos_angles) {
		// Every edge pairs consist 1 angle constraint
		const_n = edge_pairs.size();
		IJV.resize(12*const_n);

		for (int i = 0; i < edge_pairs.size(); i++) {
			int v1_i(edge_pairs[i].first.v1),v2_i(edge_pairs[i].first.v2),w1_i(edge_pairs[i].second.v1),w2_i(edge_pairs[i].second.v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

			double l1 = sqrt(pow(v1_x-v2_x,2)+pow(v1_y-v2_y,2)+pow(v1_z-v2_z,2));
			double l2 = sqrt(pow(w1_x-w2_x,2)+pow(w1_y-w2_y,2)+pow(w1_z-w2_z,2));
			e_lens.push_back(l1*l2);
		};
	};

	void set_angles(const std::vector<double> cos_angles_i) { cos_angles = cos_angles_i;}
	virtual EdgesAngleConstraints* clone() const {return new EdgesAngleConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(edge_pairs.size()); 
		for (int i = 0; i < edge_pairs.size(); i++) {
			int v1_i(edge_pairs[i].first.v1),v2_i(edge_pairs[i].first.v2),w1_i(edge_pairs[i].second.v1),w2_i(edge_pairs[i].second.v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

			double e_len = e_lens[i];

			double cos_angle = cos_angles[i];
			vals(i) = -cos_angle+((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z))/e_len;
			//std::cout << "dot prod = " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z)) << std::endl;
			//std::cout << "cos_angle = " << cos_angle << std::endl;
			//std::cout << "normalized it's: " << ((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z))/e_len << std::endl;
			//exit(1);
		};
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;

		Eigen::VectorXd vals(edge_pairs.size()); 
		for (int i = 0; i < edge_pairs.size(); i++) {
			int v1_i(edge_pairs[i].first.v1),v2_i(edge_pairs[i].first.v2),w1_i(edge_pairs[i].second.v1),w2_i(edge_pairs[i].second.v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(x(w1_i)); const double w1_y(x(w1_i+1*vnum)); const double w1_z(x(w1_i+2*vnum));
			const double w2_x(x(w2_i)); const double w2_y(x(w2_i+1*vnum)); const double w2_z(x(w2_i+2*vnum));

			double cos_angle = cos_angles[i];
			double e_len = e_lens[i];
			//std::cout << "e_len = " << e_len << std::endl;

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i, (w1_x-w2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+vnum, (w1_y-w2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+2*vnum, (w1_z-w2_z)/e_len);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i, (-w1_x+w2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+vnum, (-w1_y+w2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+2*vnum, (-w1_z+w2_z)/e_len);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i, (v1_x-v2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i+vnum, (v1_y-v2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w1_i+2*vnum, (v1_z-v2_z)/e_len);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i, (-v1_x+v2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i+vnum, (-v1_y+v2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, w2_i+2*vnum, (-v1_z+v2_z)/e_len);

  			const_cnt++;
  		}
  		if (const_cnt != const_n) {
			std::cout << "EdgesAngleConstraints: Error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << std::endl;
			exit(1);
		}
		if (ijv_cnt != IJV.size()) {
			std::cout << "EdgesAngleConstraints: Error! ijv_cnt = " << ijv_cnt << " but IJV.size() " << IJV.size() << std::endl;
			exit(1);
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// The constraints are not linear but for now we use here gauss-newton..
	};

	int vnum;
	std::vector<std::pair<Edge,Edge>> edge_pairs;
	std::vector<double> cos_angles;
	std::vector<double> e_lens;
};