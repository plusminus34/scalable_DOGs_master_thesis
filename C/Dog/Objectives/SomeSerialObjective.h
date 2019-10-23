#pragma once

#include "../../Optimization/Objective.h"

class SomeSerialObjective: public Objective {

public:
  SomeSerialObjective();
	SomeSerialObjective(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &lambda_i);
	virtual SomeSerialObjective* clone() const {return new SomeSerialObjective(*this);}

  void set_pointers(const Eigen::SparseMatrix<double> &A_i, const Eigen::VectorXd &lambda_i);

	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);

  void initialize();

	bool ready;
	const Eigen::SparseMatrix<double> *A;
	const Eigen::VectorXd *lambda;
};
