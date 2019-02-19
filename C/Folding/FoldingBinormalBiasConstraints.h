#pragma once

#include "../Dog/Dog.h"
#include "../Optimization/Constraints.h"

class FoldingBinormalBiasConstraints: public Constraints {
public:
	FoldingBinormalBiasConstraints(const Dog& dog);
	virtual FoldingBinormalBiasConstraints* clone() const {return new FoldingBinormalBiasConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda);

private:
	const Dog& dog;
	const DogEdgeStitching& eS;
};
