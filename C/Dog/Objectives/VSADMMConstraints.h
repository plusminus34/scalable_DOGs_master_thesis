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

	//Give addresses of A,z,lambda (and also the double rho)
	virtual void set_pointers(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i, const Eigen::VectorXd &lambda_i, double rho_i);

private:
	//Delayed initialization because A and others aren't ready during construction
	void initialize();

	bool ready;
	const Eigen::SparseMatrix<double> *A;
	const Eigen::VectorXd *z;
	const Eigen::VectorXd *lambda;
	double rho;
};
