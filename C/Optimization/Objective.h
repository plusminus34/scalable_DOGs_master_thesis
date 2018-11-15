#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Objective {
  
public:
  virtual Objective* clone() const = 0;

  virtual double obj(const Eigen::VectorXd& x) const = 0;
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const = 0;

  // return 0 sparse matrix if not implemented
  virtual Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const {return Eigen::SparseMatrix<double>(x.rows(),x.rows());}

  // For objectives that depeends on a state/reference shape (which might also be updated at various times)
  virtual void set_ref(const Eigen::VectorXd& x0) {};

  void check_grad(const Eigen::VectorXd& x) const;

 private:

  void finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy = 1) const;
};