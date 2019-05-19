#pragma once

#include "../Dog/Dog.h"
#include "../Optimization/Constraints.h"

class FoldingMVBiasConstraints: public Constraints {
public:
	FoldingMVBiasConstraints(const Dog& dog, bool flip_sign, int curve_i = 0);
	FoldingMVBiasConstraints(const Dog& dog, const std::vector<bool>& flip_sign, const std::vector<int>& curve_i);
	virtual FoldingMVBiasConstraints* clone() const {return new FoldingMVBiasConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda);

	bool is_mv_assignment_correct(const Eigen::VectorXd& x);

private:
	const Dog& dog;
	const DogEdgeStitching& eS;

	std::vector<int> curve_indices; std::vector<double> flip_signs_mult;
	int vnum;
	double delta = 0.00001;
};
