#pragma once

#include "../../Optimization/Solver.h"
#include "../DogLaplacian.h"

// This just used the linearized constraints but the lagrangian but solves [H,J;J^t,0] 
//  meaning H is just the obj Hessian without the second order parts of the constraints
class NewtonKKT : public ConstrainedSolver {
  
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
    x = x0;
    int vnum = x.rows()/3;
    double new_e;

    // Get metric
    Eigen::MatrixXd V_x; vec_to_mat2(x,V_x); 
    Eigen::SparseMatrix<double> H = obj.hessian(x) + 1e-2*id;
    
    //energy->check_grad(x);
    double old_e = f.obj(x);

    Eigen::VectorXd g(f.grad(x));
    Eigen::VectorXd d(g.rows());
    Eigen::SparseMatrix<double> J = constraints.Jacobian(x);
    
    Eigen::SparseMatrix<double> Jt = J.transpose();
    Eigen::SparseMatrix<double> L_jt; igl::cat(2,metric,Jt, L_jt);          

    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    Eigen::SparseMatrix<double> A; igl::cat(1, L_jt, J_0, A);
    
    A.makeCompressed();
    Eigen::SparseMatrix<double> id_all(A.rows(),A.rows()); id_all.setIdentity();
    A = A + 0*id_all; // todo: stupid but I want to add zeros explicitly

    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    cout << "factorizing" << endl;
    //solver.factorize(A);
    solver.compute(H);
    if(solver.info()!=Eigen::Success) {
        cout << "Eigen Failure!" << endl;
        exit(1);
    }
    //m_solver.factorize();

    Eigen::VectorXd zeroV(J.rows()); zeroV.setZero();
    Eigen::VectorXd constraints_deviation = -1*constraints.Vals(x);
    Eigen::VectorXd g_const; igl::cat(1, g, constraints_deviation, g_const);
    
    Eigen::VectorXd res;
    //cout << "solving!" << endl;
    res = solver.solve(g_const);
    
    for (int d_i = 0; d_i < g.rows(); d_i++) {
        d[d_i] = -1*res[d_i];
    }
    new_e = line_search(x,d,flow_t,f);
    
    old_e = f.obj(x);
    
    return new_e;
    }


    
};