#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

class PointPairConstraints : public Constraints {
public:
	PointPairConstraints(const std::vector<std::pair<int,int>>& pairs) : pairs(pairs) {
		// Here the pairs are assumed to be defined on the vector pairs (so pairing up points means having 3 index pairs in the c'tor here)
		const_n = pairs.size();
		IJV.resize(2*const_n);
	};

	virtual PointPairConstraints* clone() const {return new PointPairConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd vals(pairs.size()); 
		for (int i = 0; i < pairs.size(); i++) { vals(i) = x(pairs[i].first)-x(pairs[i].second); };
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_n = 0; int ijv_cnt = 0;
		for (int b_i = 0; b_i < pairs.size(); b_i++ ) {
			int var_const_idx1 = pairs[b_i].first; int var_const_idx2 = pairs[b_i].second;
			// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_n, var_const_idx1, 1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_n, var_const_idx2, -1);
			const_n++;
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

	std::vector<std::pair<int,int>> pairs;
};