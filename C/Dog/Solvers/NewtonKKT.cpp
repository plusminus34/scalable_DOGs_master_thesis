#include "NewtonKKT.h"

#include "../../QuadMesh/Quad.h"

#include "igl/cat.h"

using namespace std;

double NewtonKKT::solve_constrained(const Eigen::VectorXd& x0, Objective& f, const Constraints& constraints, Eigen::VectorXd& x) {
	x = x0;
    int vnum = x.rows()/3;
    double new_e;

    // Get hessian
    Eigen::MatrixXd V_x; vec_to_mat2(x,V_x); 
    Eigen::SparseMatrix<double> id(x.rows(),x.rows()); id.setIdentity();
    Eigen::SparseMatrix<double> H = f.hessian(x) + 1e-2*id;
    //Eigen::SparseMatrix<double> H = 1e-2*id;//;f.hessian(x) + 1e-2*id;
    
    //energy->check_grad(x);
    double old_e = f.obj(x);

    Eigen::VectorXd g(f.grad(x));
    Eigen::VectorXd d(g.rows());
    Eigen::SparseMatrix<double> J = constraints.Jacobian(x);
    
    Eigen::SparseMatrix<double> Jt = J.transpose();
    Eigen::SparseMatrix<double> H_jt; igl::cat(2,H,Jt, H_jt);          

    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    Eigen::SparseMatrix<double> A; igl::cat(1, H_jt, J_0, A);
    
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
    double init_t = 1;
    new_e = line_search(x,d,init_t,f);
    
    old_e = f.obj(x);
    
    return new_e;
}