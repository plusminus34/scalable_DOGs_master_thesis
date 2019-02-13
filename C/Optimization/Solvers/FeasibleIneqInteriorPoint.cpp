#include "FeasibleIneqInteriorPoint.h"

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/line_search.h"

#include "igl/Timer.h"

#include "igl/cat.h"
#include "igl/matlab_format.h"

using namespace std;

double FeasibleIneqInteriorPoint::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints,
            Constraints& ineq_constraints, Eigen::VectorXd& x) {
    x = x0;
    if (!is_feasible) {
        // check feasibility
        get_feasible_point(x,obj,eq_constraints,ineq_constraints,x);
    }

    double mu = 1; double mu_sigma = 0.1;
    while (kkt_mu_error(x,obj, eq_constraints, ineq_constraints, 0 ) < tol) {
        while (kkt_mu_error(x,obj, eq_constraints, ineq_constraints, mu ) < mu) { // tol_mu = mu as in Nocedal
            single_homotopy_iter(x, obj, eq_constraints, ineq_constraints, mu, x, current_merit);
        }
        mu = mu*mu_sigma;
    }
    /*
    auto prev_x = x0; double current_merit_p = merit_p;
    double ret = compute_step_direction(x0,f,constraints,x,current_merit_p); int iter = 1;
    double const_dev = constraints.Vals(x).norm();
    if (const_dev > infeasability_filter) { x = prev_x; current_merit_p *=2;} prev_x = x;
    while ((const_dev > infeasability_epsilon) && iter < max_newton_iters) {
        std::cout << "const_dev = " << const_dev << " after " << iter  << " iters" << std::endl;
        ret = compute_step_direction(x,f,constraints,x,current_merit_p);
        const_dev = constraints.Vals(x).norm();
        if (const_dev > infeasability_filter) { x = prev_x; current_merit_p *=2;} prev_x = x;
        iter++;
    }
    */
    return obj.obj(x);
}
double FeasibleIneqInteriorPoint::get_feasible_point(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints,
            Constraints& ineq_constraints, Eigen::VectorXd& x) {
    const double MAX_FEASIBLE_ITER = 10;
    double mu = 0.5; int iter = 0;
    auto ineq_vals = ineq_constraints.Vals(x);
    is_feasible = ineq_vals.minCoeff() > 0;
    while ((!is_feasible) && iter < MAX_FEASIBLE_ITER) {
        mu *=2; iter++;
        cout << "Trying to find a feasible point: iter = " << iter << ", mu = " << mu << endl;
        // Retain feasiblity, currently assume inequalities are 0 or very close to it
        single_homotopy_iter(x, obj, eq_constraints, ineq_constraints, mu, x, current_merit);
        ineq_vals = ineq_constraints.Vals(x);
        is_feasible = ineq_vals.minCoeff() > 0;
    }
    if (is_feasible) cout << "Feasible!" << endl;
    else cout << "Failure: could not reach feasible point!" << endl;
}

double FeasibleIneqInteriorPoint::single_homotopy_iter(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
                                double mu, Eigen::VectorXd& x, double current_merit) {
    Eigen::VectorXd step_d;
    compute_step_direction(x,obj,eq_constraints,ineq_constraints, mu, step_d, current_merit);
    double alpha = get_max_alpha(x,d);
    ineq_linesearch(x,d,alpha,mu,obj,eq_constraints,ineq_constraints,current_merit);
    update_variables(x,d,alpha);
    return alpha;
}

double FeasibleIneqInteriorPoint::get_max_alpha(const Eigen::VectorXd& x, const Eigen::VectorXd& d) {
    const thetha = 0.995;
    double alpha = 1;
    int s_base = x.rows()+lambda.rows();
    int z_base = s_base + s.rows();
    // fraction to the boundary rule (19.9 in Nocedal), make sure the s and z remain positive
    for (int i = 0; i < s.rows(); i++) if (d(s_base+i) < 0) alpha = std::min(alpha, ((1-thetha)*s(i)-s(i))/d(s_base+i) );
    for (int i = 0; i < z.rows(); i++) if (d(z_base+i) < 0) alpha = std::min(alpha, ((1-thetha)*z(i)-z(i))/d(z_base+i) );
    return alpha;
}

void FeasibleIneqInteriorPoint::update_variables(Eigen::VectorXd& x,const Eigen::VectorXd& d, Constraints& ineq_constraints,
            double alpha) {

    int lambda_base = x.rows();
    int s_base = lambda_base+lambda.rows();
    int z_base = s_base + s.rows();
    for (int i = 0; i < x.rows(); i++) x(i) += alpha*d(i);
    for (int i = 0; i < lambda.rows(); i++) lambda(i) += alpha*d(lambda_base+i);
    //for (int i = 0; i < s.rows(); i++) s(i) += alpha*d(s_base+i);
    s = ineq_constraints.Vals(); // instead of line search, so that the merit function will get infinite
    for (int i = 0; i < z.rows(); i++) z(i) += alpha*d(z_base+i);
}

double FeasibleIneqInteriorPoint::merit_func(Eigen::VectorXd& x, double mu, Objective& f, 
        Constraints& eq_constraints, Constraints& ineq_constraints, double merit) {
    // Get the log of the inequalities
    auto ineq_log_vals = ineq_constraints.Vals(x);
    for (int i = 0; i <ineq_vals.rows(); i++) {
        ineq_log_vals(i) = log(max(0,ineq_log_vals(i)));
    }

    return f.obj(x) -mu*ineq_log_vals.sum() + merit*eq_constraints.Vals(x).norm()+merit*(ineq_constraints.Vals(x)-s).norm();
}

double FeasibleIneqInteriorPoint::ineq_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, double mu, Objective& f, 
        Constraints& eq_constraints, Constraints& ineq_constraints, double merit) {
  double old_energy = merit_func(x, mu, f, eq_constraints, ineq_constraints, merit);
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = 22;

  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::VectorXd new_x = x + step_size * d;

    double cur_e = merit_func(new_x, mu, f, eq_constraints, ineq_constraints, merit);
    
    //cout << "cur_e = " << cur_e << endl;
    if (cur_e >= old_energy) {
      step_size /= 2;
      //cout << "step_size = " << step_size << endl;
    } else {
      x = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  if (cur_iter < MAX_STEP_SIZE_ITER) {
    //cout << "ls success!" << endl;
  } else {
    cout << "ls failure, is it a local minimum?" << endl;
    return old_energy;
  }
  //cout << "step = " << step_size << endl;
  return new_energy;
}

double FeasibleIneqInteriorPoint::kkt_mu_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
        double mu) {

    auto jacobian = eq_constraints.Jacobian(x);
    auto jacobian_ineq = ineq_constraints.Jacobian(x);
    auto g = obj.grad(x0)
    double grad_error = (g-jacobian.transpose()*lambda-jacobian_ineq.transpose()*z).norm();

    double sz_minus_m_error = 0;
    for (int i = 0; i < s.rows(); i++) sz_minus_m_error += sqrt(pow(s(i)*z(i)-mu,2));
    double eq_const_error = eq_constraints.Vals().norm();
    double ineq_const_error = (ineq_constraints.Vals()-s).norm();
    return std::max({grad_error, sz_minus_m_error, eq_const_error, ineq_const_error});
}
void FeasibleIneqInteriorPoint::build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, 
                                     const std::vector<Eigen::Triplet<double> >& const_lambda_hessian, int var_n,
                                     const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n,
                                     const std::vector<Eigen::Triplet<double> >& ineq_lambda_hessian,
                                     const std::vector<Eigen::Triplet<double> >& ineq_jacobian_ijv, int ineq_const_n) {
    // build an IJV which contains the hessian, diag matrix along the hessian, jacobian, jacobian transpose,
    //  and diagonal matrix along the J to make sure the whole diagonal is nnz (needed for Pardiso solver)
    int ijv_idx = 0;
    if (kkt_IJV.size() == 0) {
        int kkt_size = hessian_IJV.size() + const_lambda_hessian.size() + 2*jacobian_IJV.size() + var_n + const_n;
        kkt_size += ineq_lambda_hessian.size() + 2*ineq_jacobian_ijv.size() + s.size() + 2*z.size();
        kkt_IJV.resize(kkt_size); int ijv_idx = 0;
    }

    // add hessian+const_lambda_hessian+id*eps_id
    double eps = 1e-8;
    for (int i = 0; i < hessian_IJV.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(hessian_IJV[i].row(),hessian_IJV[i].col(),-hessian_IJV[i].value());
    }
    for (int i = 0; i < const_lambda_hessian.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(const_lambda_hessian[i].row(),const_lambda_hessian[i].col(),const_lambda_hessian[i].value());
    }
    for (int i = 0; i < ineq_lambda_hessian.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(ineq_lambda_hessian[i].row(),ineq_lambda_hessian[i].col(),ineq_lambda_hessian[i].value());
    }
    for (int i = 0; i < var_n; i++) { kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,-eps);}

    // Add both J and J transpose
    for (int i = 0; i < jacobian_IJV.size(); i++) {
        int row = jacobian_IJV[i].row(), col = jacobian_IJV[i].col(); double val = jacobian_IJV[i].value();
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(col,row+var_n, val); // Jt columed offseted at var_n
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(row+var_n, col, val); // J offseted in var_n
    }
    // add zeros along var_n+i,var_n_i
    for (int i = var_n; i < var_n+const_n; i++) kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,0);

    if ( A.rows() == 0) {
      A =  Eigen::SparseMatrix<double>(var_n+const_n,var_n+const_n);
      igl::sparse_cached_precompute(kkt_IJV, cached_ijv_data, A);
    } else {
      igl::sparse_cached(kkt_IJV, cached_ijv_data, A);
    }
}

void FeasibleIneqInteriorPoint::build_kkt_system(const Eigen::SparseMatrix<double>& hessian,
                        const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& KKT) {
    Eigen::SparseMatrix<double> id (hessian.rows(),hessian.cols()); id.setIdentity();
    auto eps_id = 1e-8*id;

    Eigen::SparseMatrix<double> H = -hessian - eps_id;
    
    Eigen::SparseMatrix<double> Jt = J.transpose();

    Eigen::SparseMatrix<double> H_jt; igl::cat(2,H,Jt, H_jt);
    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    igl::cat(1, H_jt, J_0, KKT);

    //A.makeCompressed();
    Eigen::SparseMatrix<double> id_KKT(KKT.rows(),KKT.rows()); id_KKT.setIdentity(); id_KKT = 0*id_KKT;
    KKT = KKT + id_KKT; // todo: stupid but Paradiso wants to add zeros explicitly
}

double FeasibleIneqInteriorPoint::compute_step_direction(const Eigen::VectorXd& x0, Objective& f, Constraints& eq_constraints, Constraints& ineq_constraints, 
        double mu, Eigen::VectorXd& x, double current_merit) {
    x = x0;
    int vnum = x.rows()/3;
    double new_e;

    if (lambda.rows()!= eq_constraints.getConstNum()) {
        lambda.resize(eq_constraints.getConstNum());
        lambda.setZero();
        // init ineq slacks 's'
        s.resize(ineq_constraints.getConstNum()); for (int i = 0; i < s.rows(); i++) s(i) = 1e-3;
        // init slacks multipliers 'z'
        z.resize(s.rows()); for (int i = 0; i < z.rows(); i++) z(i) = mu/s(i);
    }
    std::cout << "lagrange multipliers norm = " << lambda.norm() << std::endl;

    igl::Timer timer; auto init_time = timer.getElapsedTime(); auto t = init_time;
    // Get Hessian
    //auto hessian = f.hessian(x); 
    auto hessian_ijv = f.update_and_get_hessian_ijv(x);
    auto lambda_hessian_ijv = eq_constraints.update_and_get_lambda_hessian(x,lambda);
    auto lambda_ineq_hessian_ijv = ineq_constraints.update_and_get_lambda_hessian(x,lambda);
    auto hessian_time = timer.getElapsedTime()-t;

    // Get Jacobian
    t = timer.getElapsedTime();
    auto jacobian = eq_constraints.Jacobian(x);
    auto jacobian_ineq = ineq_constraints.Jacobian(x);
    auto jacobian_ijv = eq_constraints.update_and_get_jacobian_ijv(x);
    auto ineq_jacobian_ijv = ineq_constraints.update_and_get_jacobian_ijv(x);
    auto jacobian_time = timer.getElapsedTime()-t;

    t = timer.getElapsedTime();
    //Eigen::SparseMatrix<double> A; build_kkt_system(hessian,jacobian,A);
    int var_n = x.rows(); int const_n = eq_constraints.getConstNum(); int ineq_const_n = ineq_constraints.getConstNum();
    build_kkt_system_from_ijv(hessian_ijv, lambda_hessian_ijv, var_n,
                        jacobian_ijv, const_n, lambda_ineq_hessian_ijv, ineq_jacobian_ijv, ineq_const_n);

    auto kkt_time = timer.getElapsedTime()-t;  

    t = timer.getElapsedTime();

    //energy->check_grad(x);
    double old_e = f.obj(x);
    Eigen::VectorXd g(f.grad(x));
    Eigen::VectorXd d(g.rows()); Eigen::VectorXd lambda_d(lambda.rows());
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    //cout << "analayzing pattern" << endl;
    //solver.analyzePattern(A);
    if (first_solve) {
        m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
        m_solver.set_pattern();
        m_solver.analyze_pattern();
        //m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
        //m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
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

    Eigen::VectorXd constraints_deviation = -1*eq_constraints.Vals(x);
    Eigen::VectorXd rhs_upper = g-jacobian.transpose()*lambda-jacobian_ineq.transpose()*z;
    Eigen::VectorXd g_const; igl::cat(1, rhs_upper, constraints_deviation, g_const);
    //g_const = -1*g_const;
    
    Eigen::VectorXd res;
    //cout << "solving!" << endl;
    //res = solver.solve(g_const);
    m_solver.solve(g_const,res);
    
    for (int d_i = 0; d_i < g.rows(); d_i++) {d[d_i] = res[d_i];}
    for (int d_i = g.rows(); d_i < res.rows(); d_i++) {lambda_d[d_i-g.rows()] = res[d_i];}
    auto solve_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    double step_size = 1;
    //new_e = line_search(x,d,init_t,f);
    new_e = exact_l2_merit_linesearch(x,d,step_size,f,eq_constraints,current_merit);
    // update lagrange multipliers
    lambda = lambda + step_size*lambda_d;
    auto linesearch_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();


    double total_time = timer.getElapsedTime()-init_time;
    /*
    cout << endl << endl << "total_kkt_system_time  = " << total_time << endl;
    cout << "Hessian_compute_time  = " << hessian_time << endl;
    cout << "Jacobian_compute_time  = " << jacobian_time << endl;
    cout << "kkt_system_build_time  = " << kkt_time << endl;
    
    cout << "hessian analyze_pattern_time = " << analyze_pattern_time << endl;
    cout << "factorize_time = " << factorize_time << endl;
    cout << "solve_time  = " << solve_time << endl;
    cout << "linesearch_time  = " << linesearch_time << endl;
    */
    //std::ofstream out_file(std::string("KKT_mat"));
    //out_file << igl::matlab_format(A,"A");
    //out_file << igl::matlab_format(g_const,"g_const");
    //std::ofstream out_file(std::string("null_space.m"));
    //out_file << igl::matlab_format(Jt,"Jt");
    
    //std::cout << "NewtonKKT: old_e = " << old_e << " new_e = " << new_e << std::endl;
    
    return new_e;
}