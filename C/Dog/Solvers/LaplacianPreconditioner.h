#pragma once

#include "../../Optimization/Solver.h"
#include "../DogLaplacian.h"

class LaplacianPreconditioner : public Solver {
  
public:
	LaplacianPreconditioner(Dog& dog) {
        Eigen::SparseMatrix<double> preconditioner = DOG_laplacian(dog.getV(),dog.getF());
        Eigen::SparseMatrix<double> id(preconditioner.rows(),preconditioner.rows()); id.setIdentity();
        preconditioner = preconditioner - (1e-8)*id;    
    
        solver.compute(preconditioner);
    }
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) {
        // todo do a linesearch
        x = x0;
        for (int i = 0; i < 1000; i++) {
            double alpha = 1.;
            x = x0+alpha*solver.solve(obj.grad(x0));
        }
    }
private:
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
};