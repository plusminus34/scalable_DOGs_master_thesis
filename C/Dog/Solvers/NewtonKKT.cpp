#include "NewtonKKT.h"

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/line_search.h"

#include "igl/Timer.h"

#include "igl/cat.h"
#include "igl/matlab_format.h"

using namespace std;
//m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
                //m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
                //m_solver.iparm[20] = 1;
double NewtonKKT::solve_constrained(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x) {
	x = x0;
    int vnum = x.rows()/3;
    double new_e;

    if (id.rows() == 0) {
        id = Eigen::SparseMatrix<double>(x.rows(),x.rows()); id.setIdentity();
        eps_id = 1e-7*id;
    }

    igl::Timer timer; auto init_time = timer.getElapsedTime(); auto t = init_time;
    // Get hessian
    Eigen::SparseMatrix<double> H = -f.hessian(x) - eps_id;

    auto hessian_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    
    //energy->check_grad(x);
    double old_e = f.obj(x);

    Eigen::VectorXd g(f.grad(x));
    Eigen::VectorXd d(g.rows());
    auto J_time = timer.getElapsedTime();
    Eigen::SparseMatrix<double> J = constraints.Jacobian(x);
    J_time = timer.getElapsedTime()-J_time;
    auto transpose_time = timer.getElapsedTime();
    Eigen::SparseMatrix<double> Jt = J.transpose();
    transpose_time = timer.getElapsedTime()-transpose_time;
    auto kkt_init = timer.getElapsedTime();
    Eigen::SparseMatrix<double> H_jt; igl::cat(2,H,Jt, H_jt);

    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    Eigen::SparseMatrix<double> A; igl::cat(1, H_jt, J_0, A);
    auto kkt_without_J = timer.getElapsedTime()-kkt_init;

    //A.makeCompressed();
    if (!id_KKT.rows()) {id_KKT = Eigen::SparseMatrix<double>(A.rows(),A.rows()); id_KKT.setIdentity(); id_KKT = 0*id_KKT;}
    auto id_add_time = timer.getElapsedTime();
    A = A + id_KKT; // todo: stupid but Paradiso wants to add zeros explicitly
    id_add_time = timer.getElapsedTime()-id_add_time;
    auto kkt_time = timer.getElapsedTime()-t;

    t = timer.getElapsedTime();
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    //cout << "analayzing pattern" << endl;
    //solver.analyzePattern(A);
    if (first_solve) {
        m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
        m_solver.set_pattern();
        m_solver.analyze_pattern();
        first_solve = false;
    } else {
        m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    auto analyze_pattern_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    m_solver.factorize();
    /*solver.factorize(A);
    if(solver.info()!=Eigen::Success) {
        cout << "Eigen Failure!" << endl;
        exit(1);
    }*/
    auto factorize_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    
    //m_solver.factorize();

    Eigen::VectorXd zeroV(J.rows()); zeroV.setZero();
    Eigen::VectorXd constraints_deviation = -1*constraints.Vals(x);
    Eigen::VectorXd g_const; igl::cat(1, g, constraints_deviation, g_const);
    //g_const = -1*g_const;
    
    Eigen::VectorXd res;
    //cout << "solving!" << endl;
    //res = solver.solve(g_const);
    m_solver.solve(g_const,res);
    
    for (int d_i = 0; d_i < g.rows(); d_i++) {
        d[d_i] = res[d_i];
    }
    auto solve_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    double init_timestep = 1;
    //new_e = line_search(x,d,init_t,f);
    new_e = exact_l2_merit_linesearch(x,d,init_timestep,f,constraints,merit_p);
    auto linesearch_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();


    double total_time = timer.getElapsedTime()-init_time;
    
    cout << endl << endl << "total kkt system time  = " << total_time << endl;
    cout << "hessian compute time  = " << hessian_time << endl;
    cout << "total kkt_system_build_time  = " << kkt_time << endl;
    cout << "\t Out of it J time was  = " << J_time << endl;
    cout << "\t Transposing J time was  = " << transpose_time << endl;
    cout << "\t kkt_system_build_time_without jacobian  = " << kkt_without_J << endl;
    cout << "\t id_add_time time was  = " << id_add_time << endl;
    cout << "hessian analyze_pattern_time = " << analyze_pattern_time << endl;
    cout << "factorize_time = " << factorize_time << endl;
    cout << "solve_time  = " << solve_time << endl;
    cout << "linesearch_time  = " << linesearch_time << endl;
    
    //std::ofstream out_file(std::string("KKT_mat"));
    //out_file << igl::matlab_format(A,"A");
    //out_file << igl::matlab_format(g_const,"g_const");
    //std::ofstream out_file(std::string("null_space.m"));
    //out_file << igl::matlab_format(Jt,"Jt");
    
    old_e = f.obj(x);
    
    return new_e;
}