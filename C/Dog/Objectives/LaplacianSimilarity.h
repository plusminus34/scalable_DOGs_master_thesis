#pragma once

#include "../../Optimization/Objective.h"

#include "../Dog.h"

// comparing to squared length
class LaplacianSimilarity: public Objective {
  
public:
	LaplacianSimilarity(const Dog& dog, const Eigen::VectorXd& x0);
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;
	//virtual void set_ref(const Eigen::VectorXd& x0);
	virtual Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const;

	//Eigen::VectorXd refA,refB;
private:
	Eigen::SparseMatrix<double> L;
	Eigen::VectorXd lRef;
};