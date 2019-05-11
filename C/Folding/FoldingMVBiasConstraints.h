#pragma once

#include "../Dog/Dog.h"
#include "../Optimization/Constraints.h"

class FoldingMVBiasConstraints: public Constraints {
public:
	FoldingMVBiasConstraints(const Dog& dog, bool flip_sign, int curve_i = 0);
	virtual FoldingMVBiasConstraints* clone() const {return new FoldingMVBiasConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda);

private:
	const Dog& dog;
	const DogEdgeStitching& eS;
	int vnum;
	int curve_i;
	double flip_sign_mult = 1;
	double delta = 0.00001;
};
