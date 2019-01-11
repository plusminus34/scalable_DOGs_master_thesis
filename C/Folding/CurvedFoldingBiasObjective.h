#pragma once

#include "../Dog/Dog.h"
#include "../Optimization/Objective.h"

struct CurvedFoldBias {
	EdgePoint ep_0, ep_b, ep_f;
	int v1,v2;
};

// comparing to squared length
class CurvedFoldingBiasObjective: public Objective {
public:
	CurvedFoldingBiasObjective(const Dog& dog, std::vector<CurvedFoldBias>& curvedFoldBiases) : dog(dog), curvedFoldBiases(curvedFoldBiases) {}
	virtual CurvedFoldingBiasObjective* clone() const {return new CurvedFoldingBiasObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x) {}; // For now no hessian
	
	const Dog& dog;
	std::vector<CurvedFoldBias> curvedFoldBiases;
};