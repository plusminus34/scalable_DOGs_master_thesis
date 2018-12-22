#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>

class Objective {
  
public:
  virtual Objective* clone() const = 0;

  virtual double obj(const Eigen::VectorXd& x) const = 0;
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const = 0;

  // builds the hessian from an IJV
  // Can be overloaded by a function that build it by another way, not through IJV
  virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x) {
  	//Eigen::SparseMatrix<double> H(x.rows(),x.rows());
    updateHessianIJV(x);
    if ( cachedH.rows() == 0) {
      cachedH =  Eigen::SparseMatrix<double>(x.rows(),x.rows());
      igl::sparse_cached_precompute(IJV, cached_ijv_data, cachedH);
    } else {
      igl::sparse_cached(IJV, cached_ijv_data, cachedH);
    }
  	return cachedH;
  }
  
  const std::vector<Eigen::Triplet<double> >& update_and_get_hessian_ijv(const Eigen::VectorXd& x) {
    updateHessianIJV(x); return IJV; 
  }
  // For objectives that depeends on a state/reference shape (which might also be updated at various times)
  virtual void set_ref(const Eigen::VectorXd& x0) {};

  void check_grad(const Eigen::VectorXd& x) const;

 protected:
   std::vector<Eigen::Triplet<double> > IJV;
   Eigen::SparseMatrix<double> cachedH;
 private:
  // return 0 sparse matrix if not implemented
  virtual void updateHessianIJV(const Eigen::VectorXd& x) { /*empty on purpose */ } 

  Eigen::VectorXi cached_ijv_data;

  void finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy = 1) const;
};