#pragma once

#include "../../Optimization/Solver.h"
#include "../DogLaplacian.h"

class LaplacianPreconditioner : public Solver {
  
public:
	LaplacianPreconditioner(Dog& dog) {
        Eigen::SparseMatrix<double> preconditioner = DOG_laplacian(dog.getV(),dog.getF());
        Eigen::SparseMatrix<double> id(preconditioner.rows(),preconditioner.rows()); id.setIdentity();
        preconditioner = 1e-8*id- DOG_laplacian(dog.getV(),dog.getF());///*preconditioner*/ - 1*id;
        dogLaplacianId = preconditioner;
    
        solver.compute(preconditioner);
        if(solver.info()!=Eigen::Success) {
            std::cout << "Eigen Failure!" << std::endl;
            exit(1);
        }
    }
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) {

        Eigen::SparseMatrix<double> preconditioner = dogLaplacianId + obj.hessian(x);

        // todo do a linesearch
        x = x0; double new_e; const double init_alpha = 1.;
        for (int i = 0; i < 1000; i++) {
            auto d = -1*solver.solve(obj.grad(x));
            //std::cout << "d.norm() = " << d.norm() << std::endl;
            new_e = line_search(x,d,init_alpha,obj);
        }
        return new_e;
    }
private:
    Eigen::SparseMatrix<double> dogLaplacianId;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
};