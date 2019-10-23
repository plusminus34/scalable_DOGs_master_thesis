#pragma once

#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class LinearConstraints: public Constraints {
public:
	LinearConstraints();
	LinearConstraints(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i);
	virtual LinearConstraints* clone() const {return new LinearConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

	//Give addresses of A and z
	virtual void set_pointers(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &z_i);

private:
	//Delayed initialization because A and others aren't ready during construction
	void initialize();

	bool ready;
	const Eigen::SparseMatrix<double> *A;
	const Eigen::VectorXd *z;
};
