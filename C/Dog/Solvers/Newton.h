#pragma once

#include "../../Optimization/Solver.h"
#include "../DogLaplacian.h"

class Newton : public Solver {
  
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) {
        x = x0;
        Eigen::SparseMatrix<double> id(x.rows(),x.rows()); id.setIdentity();
        Eigen::SparseMatrix<double> H = obj.hessian(x) + 1e-2*id;
        //Eigen::SparseMatrix<double> H = 1e-2*id;
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(H);
        if(solver.info()!=Eigen::Success) {
            std::cout << "Eigen Failure!" << std::endl;
            exit(1);
        }
        auto d = -1*solver.solve(obj.grad(x));
        double init_step_size = 1;
        double new_e = line_search(x,d,init_step_size,obj);
        return new_e;
    }
    
};