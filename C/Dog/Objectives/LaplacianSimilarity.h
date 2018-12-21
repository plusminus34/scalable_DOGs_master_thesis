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
	//virtual void set_ref(const Eigen::VectorXd& x0);
	virtual Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const;
	//virtual std::vector<Eigen::Triplet<double> > hessianIJV(const Eigen::VectorXd& x) const;

	//Eigen::VectorXd refA,refB;
private:
	//std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M);

	Eigen::SparseMatrix<double> L;
	//std::vector<Eigen::Triplet<double>> L_hessian_IJV;
	Eigen::VectorXd lRef;
};