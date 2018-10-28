#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class SumOfQuadraticConstraints : public Objective {
  
public:
  SumOfQuadraticConstraints(Constraints& constraints) : constraints(constraints), Objective() {};

  virtual double obj(const Eigen::VectorXd& x) {return constraints.Vals(x).squaredNorm();};
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) {
  	return 2*constraints.Jacobian(x).transpose()*constraints.Vals(x);
  }
 private:
 	Constraints& constraints;
};