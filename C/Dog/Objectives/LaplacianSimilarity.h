#pragma once

#include "../../Optimization/Objective.h"

#include "../Dog.h"

// comparing to squared length
class LaplacianSimilarity: public Objective {
  
public:
	LaplacianSimilarity(const Dog& dog, const Eigen::VectorXd& x0);
	virtual LaplacianSimilarity* clone() const {return new LaplacianSimilarity(*this);}

	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

	virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x);

	//virtual void set_ref(const Eigen::VectorXd& x0);
	
	//virtual std::vector<Eigen::Triplet<double> > hessianIJV(const Eigen::VectorXd& x) const;

	//Eigen::VectorXd refA,refB;
private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);

	Eigen::SparseMatrix<double> L;
	Eigen::VectorXd lRef;
};