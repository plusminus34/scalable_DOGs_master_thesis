#pragma once

#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class VSADMMConstraints: public Constraints {
public:
	VSADMMConstraints();
	VSADMMConstraints(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i);
	virtual VSADMMConstraints* clone() const {return new VSADMMConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

	virtual void setstuff(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i);

private:
	/*TODO what do I store?*/
	//or: where and how are A,z,etc stored?
	//Eigen::MatrixXd A;
	bool initd;
	void initnow();
	const Eigen::SparseMatrix<double> *A;
	const Eigen::VectorXd *z;
	const Eigen::VectorXd *lambda;
	double rho;
	/*
	inherited from constraints
	protected:
		int const_n;
		std::vector<Eigen::Triplet<double> > IJV;
		std::vector<Eigen::Triplet<double> > lambda_hessian_IJV;
	*/
};
