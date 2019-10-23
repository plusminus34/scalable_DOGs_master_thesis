#pragma once

#include "../../Optimization/Objective.h"

class ProximalObjective: public Objective {

public:
  ProximalObjective();
	ProximalObjective(const Eigen::SparseMatrix<double> &P_i, const Eigen::VectorXd &x_old_i);
	virtual ProximalObjective* clone() const {return new ProximalObjective(*this);}

  void set_pointers(const Eigen::SparseMatrix<double> &P_i, const Eigen::VectorXd &x_old_i);

	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);

  void initialize();

	bool ready;
	const Eigen::SparseMatrix<double> *P;
	const Eigen::VectorXd *x_old;
};
