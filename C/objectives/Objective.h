#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Objective {
  
public:
  Objective(){};

  virtual double obj(const Eigen::VectorXd& x) = 0;
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) = 0;

  // return 0 sparse matrix if not implemented
  //virtual Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) {return Eigen::SparseMatrix<double>(x.rows(),x.rows());}

  void check_grad(const Eigen::VectorXd& x);

 private:

  void finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy = 1);
};