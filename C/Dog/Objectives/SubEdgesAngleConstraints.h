#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "../../Optimization/Constraints.h"

class SubEdgesAngleConstraints : public Constraints {
public:
	SubEdgesAngleConstraints() : vnum(0) {
		const_n = 0;
		edges.clear();
		outside_points.resize(0,6);
		cos_angles.clear();
		IJV.resize(0);
	}

	SubEdgesAngleConstraints(int v_num, const Eigen::VectorXd& x, const std::vector<std::pair<Edge,Edge>>& edge_pairs, const std::vector<double> cos_angles) :
			vnum(v_num), cos_angles(cos_angles) {
		const_n = edge_pairs.size();
		IJV.resize(6*const_n);
		outside_points.resize(const_n, 6);
		edges.resize(edge_pairs.size());
		e_lens.resize(edge_pairs.size());
		for(int i=0; i<edge_pairs.size(); ++i){
			edges[i] = edge_pairs[i].first;
			e_lens[i] = 1.0;// placeholder
		}
	}

	void set_angles(const std::vector<double> cos_angles_i) { cos_angles = cos_angles_i;}
	virtual SubEdgesAngleConstraints* clone() const {return new SubEdgesAngleConstraints(*this);}

	void init_outside_points(const Eigen::MatrixXd& V, const Eigen::VectorXd& x){
		outside_points = V;
		e_lens.clear();
		for (int i = 0; i < edges.size(); i++) {
			int v1_i(edges[i].v1),v2_i(edges[i].v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(outside_points(i,0)); const double w1_y(outside_points(i,1)); const double w1_z(outside_points(i,2));
			const double w2_x(outside_points(i,3)); const double w2_y(outside_points(i,4)); const double w2_z(outside_points(i,5));

			double l1 = sqrt(pow(v1_x-v2_x,2)+pow(v1_y-v2_y,2)+pow(v1_z-v2_z,2));
			double l2 = sqrt(pow(w1_x-w2_x,2)+pow(w1_y-w2_y,2)+pow(w1_z-w2_z,2));
			e_lens.push_back(l1*l2);
		}
	}
	virtual void update_outside_points(const Eigen::MatrixXd& V){ outside_points = V; }

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(edges.size());
		for (int i = 0; i < edges.size(); i++) {
			int v1_i(edges[i].v1),v2_i(edges[i].v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(outside_points(i,0)); const double w1_y(outside_points(i,1)); const double w1_z(outside_points(i,2));
			const double w2_x(outside_points(i,3)); const double w2_y(outside_points(i,4)); const double w2_z(outside_points(i,5));

			double e_len = e_lens[i];

			double cos_angle = cos_angles[i];
			vals(i) = -cos_angle+((v1_x-v2_x)*(w1_x-w2_x)+(v1_y-v2_y)*(w1_y-w2_y)+(v1_z-v2_z)*(w1_z-w2_z))/e_len;
			//std::cout << "angles("<<i<<"): is " << vals(i)+cos_angle << "\tshould be "<<cos_angle<<"\n";
		};
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_cnt = 0; int ijv_cnt = 0;

		Eigen::VectorXd vals(edges.size());
		for (int i = 0; i < edges.size(); i++) {
			int v1_i(edges[i].v1),v2_i(edges[i].v2);
			const double v1_x(x(v1_i)); const double v1_y(x(v1_i+1*vnum)); const double v1_z(x(v1_i+2*vnum));
			const double v2_x(x(v2_i)); const double v2_y(x(v2_i+1*vnum)); const double v2_z(x(v2_i+2*vnum));
			const double w1_x(outside_points(i,0)); const double w1_y(outside_points(i,1)); const double w1_z(outside_points(i,2));
			const double w2_x(outside_points(i,3)); const double w2_y(outside_points(i,4)); const double w2_z(outside_points(i,5));

			double cos_angle = cos_angles[i];
			double e_len = e_lens[i];

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i, (w1_x-w2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+vnum, (w1_y-w2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v1_i+2*vnum, (w1_z-w2_z)/e_len);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i, (-w1_x+w2_x)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+vnum, (-w1_y+w2_y)/e_len);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt, v2_i+2*vnum, (-w1_z+w2_z)/e_len);

  		const_cnt++;
  	}
  	if (const_cnt != const_n) {
			std::cout << "SubEdgesAngleConstraints: Error in jacobian, const_cnt = " << const_cnt << " but const_n = " << const_n << std::endl;
			exit(1);
		}
		if (ijv_cnt != IJV.size()) {
			std::cout << "SubEdgesAngleConstraints: Error! ijv_cnt = " << ijv_cnt << " but IJV.size() " << IJV.size() << std::endl;
			exit(1);
		}
	}

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// nothing to see here
	};

	int vnum;
	std::vector<Edge> edges;
	Eigen::MatrixXd outside_points;// size: edges.size x 6
	std::vector<double> cos_angles;
	std::vector<double> e_lens;
};
