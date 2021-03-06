#include "FeasibleIneqInteriorPoint.h"

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/line_search.h"

#include "igl/Timer.h"

#include "igl/cat.h"
#include "igl/matlab_format.h"

#include <limits>

using namespace std;

double FeasibleIneqInteriorPoint::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints,
            Constraints& ineq_constraints, Eigen::VectorXd& x) {
    double current_merit_p = merit_p;
    x = x0;

    double mu = 1; double mu_sigma = 0.1;
    if (lambda.rows()!= eq_constraints.getConstNum()) {
        lambda.resize(eq_constraints.getConstNum());
        lambda.setZero();
        // init ineq slacks 's'
        s.resize(ineq_constraints.getConstNum()); for (int i = 0; i < s.rows(); i++) s(i) = 1;
        // init slacks multipliers 'z'
        z.resize(s.rows()); for (int i = 0; i < z.rows(); i++) z(i) = mu/s(i);
    }

    /*
    if (!is_feasible) {
        // check feasibility
        get_feasible_point(x,current_merit_p,obj,eq_constraints,ineq_constraints,x);
    }
    */
    cout << "Optimizing!" << endl;
    
    double opt_error = kkt_mu_error(x,obj, eq_constraints, ineq_constraints, 0 );
    cout << "opt_error = " << opt_error << endl;
    while ( opt_error > tol) {
        std::cout << "mu = " << mu << std::endl;
        double current_kkt_mu_error = kkt_mu_error(x,obj, eq_constraints, ineq_constraints,mu );
        while ( current_kkt_mu_error > mu) { // tol_mu = mu as in Nocedal
            std::cout << "current_kkt_mu_error before = " << current_kkt_mu_error << " performing another iter" << std::endl;
            //int wait; cin >> wait;
            single_homotopy_iter(x, obj, eq_constraints, ineq_constraints, mu, x, current_merit_p);
            current_kkt_mu_error = kkt_mu_error(x,obj, eq_constraints, ineq_constraints,mu );
            std::cout << "current_kkt_mu_error after = " << current_kkt_mu_error << endl;
        }
        exit(1);
        mu = mu*mu_sigma;
        opt_error = kkt_mu_error(x,obj, eq_constraints, ineq_constraints, 0 );
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
void FeasibleIneqInteriorPoint::get_feasible_point(const Eigen::VectorXd& x0, double current_merit_p, Objective& obj, Constraints& eq_constraints,
            Constraints& ineq_constraints, Eigen::VectorXd& x) {
    const double MAX_FEASIBLE_ITER = 10;
    double mu = 0.5; int iter = 0;
    auto ineq_vals = ineq_constraints.Vals(x);
    is_feasible = ineq_vals.minCoeff() > 0;
    while ((!is_feasible) && iter < MAX_FEASIBLE_ITER) {
        // updating mu also requires updating the multipliers
        mu *=2; for (int i = 0; i < z.rows(); i++) z(i) = mu/s(i); iter++;
        cout << "Trying to find a feasible point: iter = " << iter << ", mu = " << mu << endl;
        int wait; cin >> wait;
        // Retain feasiblity, currently assume inequalities are 0 or very close to it
        single_homotopy_iter(x, obj, eq_constraints, ineq_constraints, mu, x, current_merit_p);
        ineq_vals = ineq_constraints.Vals(x);
        is_feasible = ineq_vals.minCoeff() > 0;
    }
    std::cout << "Final ineq_vals = " << ineq_vals << std::endl;
    if (is_feasible) cout << "Feasible!" << endl;
    else cout << "Failure: could not reach feasible point!" << endl;
}

double FeasibleIneqInteriorPoint::single_homotopy_iter(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
                                double mu, Eigen::VectorXd& x, double current_merit) {
    Eigen::VectorXd step_d;
    //cout << "computing_step_direction" << endl;
    compute_step_direction(x,obj,eq_constraints,ineq_constraints, mu, step_d, current_merit);
    //cout << "getting alpha" << endl;
    double alpha = get_max_alpha(x,step_d);
    std::cout << "max alpha = " << alpha << std::endl;
    ineq_linesearch(x,step_d,alpha,mu,obj,eq_constraints,ineq_constraints,current_merit);
    if (alpha) update_variables(x,step_d,ineq_constraints, alpha);
    std::cout << "alpha = " << alpha << std::endl;
    return alpha;
}

double FeasibleIneqInteriorPoint::get_max_alpha(const Eigen::VectorXd& x, const Eigen::VectorXd& d) {
    const double thetha = 0.995;
    double alpha = 1;
    int s_base = x.rows();
    int z_base = s_base + lambda.rows();
    // fraction to the boundary rule (19.9 in Nocedal), make sure the s and z remain positive
    for (int i = 0; i < s.rows(); i++) if (d(s_base+i) < 0) alpha = std::min(alpha, ((1-thetha)*s(i)-s(i))/d(s_base+i) );
    for (int i = 0; i < z.rows(); i++) if (d(z_base+i) < 0) alpha = std::min(alpha, ((1-thetha)*z(i)-z(i))/d(z_base+i) );
    return alpha;
}

void FeasibleIneqInteriorPoint::update_variables(Eigen::VectorXd& x,const Eigen::VectorXd& d, Constraints& ineq_constraints,
            double alpha) {

    int s_base = x.rows();
    int lambda_base = s_base+s.rows();
    int z_base = s_base + lambda.rows();
    for (int i = 0; i < x.rows(); i++) x(i) += alpha*d(i);
    for (int i = 0; i < s.rows(); i++) s(i) += alpha*d(s_base+i);
    //s = ineq_constraints.Vals(x); // instead of line search, so that the merit function will get infinite
    for (int i = 0; i < lambda.rows(); i++) lambda(i) += alpha*d(lambda_base+i);
    for (int i = 0; i < z.rows(); i++) z(i) += alpha*d(z_base+i);
    cout << "alpha was " << alpha << endl;
    //cout << "s = " << s << endl;
    int wait; cin >> wait;
}

double FeasibleIneqInteriorPoint::merit_func(Eigen::VectorXd& new_x, Eigen::VectorXd& new_s, double mu, Objective& f, 
        Constraints& eq_constraints, Constraints& ineq_constraints, double merit) {
    // Get the log of the inequalities
    //auto ineq_log_vals = ineq_constraints.Vals(x);
    auto ineq_log_vals = new_s;
    //cout << "ineq_log_vals min coeff = " << ineq_log_vals.minCoeff() << endl;
    for (int i = 0; i < ineq_log_vals.rows(); i++) {
        ineq_log_vals(i) = log(max(0.,ineq_log_vals(i)));
    }
    
    std::cout << "f.obj(x) = " << f.obj(new_x) << " -mu*ineq_log_vals.sum() = " << -mu*ineq_log_vals.sum() << 
        " merit*eq_constraints.Vals(new_x).norm() = " << merit*eq_constraints.Vals(new_x).norm() << 
        " merit*(ineq_constraints.Vals(new_x)-new_s).norm() = " << merit*(ineq_constraints.Vals(new_x)-new_s).norm() << endl;
    return f.obj(new_x) -mu*ineq_log_vals.sum() + merit*eq_constraints.Vals(new_x).norm()+merit*(ineq_constraints.Vals(new_x)-new_s).norm();
}

double FeasibleIneqInteriorPoint::ineq_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, double mu, Objective& f, 
        Constraints& eq_constraints, Constraints& ineq_constraints, double merit) {
  double old_energy = merit_func(x,s,mu, f, eq_constraints, ineq_constraints, merit);
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = 22;

  cout << "d.norm() = " << d.norm() << endl;

  int x_rows = x.rows();
  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::VectorXd new_x = x; Eigen::VectorXd new_s = s;
    for (int i = 0; i < new_x.rows();i++) new_x(i) += step_size*d(i);
    for (int i = 0; i < new_s.rows();i++) new_s(i) += step_size*d(x_rows+i);

    double cur_e = merit_func(new_x, new_s, mu, f, eq_constraints, ineq_constraints, merit);
    cout << "step size = " << step_size << " cur_e = " << cur_e << " old_energy = " << old_energy << endl;
    
    //cout << "cur_e = " << cur_e << endl;
    //cout << " cur_e = " << cur_e << " old_energy = " << old_energy << endl;
    if (cur_e >= old_energy) {
      step_size /= 2;
      //cout << "step_size = " << step_size << endl;
    } else {
      //x = new_x; s = new_s;
      cout << "cur_e = " << cur_e << " old_energy = " << old_energy << endl;
      new_energy = cur_e;
      //exit(1);
    }
    cur_iter++;
  }
  if (cur_iter < MAX_STEP_SIZE_ITER) {
    //cout << "ls success!" << endl;
  } else {
    cout << "ls failure, is it a local minimum?" << endl; exit(1);
    //std::cout << "ineq_constraints.Vals(x): " << endl << ineq_constraints.Vals(x) << endl; exit(1);
    step_size = 0;
    return old_energy;
  }
  //cout << "step = " << step_size << endl;
  return new_energy;
}

double FeasibleIneqInteriorPoint::kkt_mu_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints, Constraints& ineq_constraints,
        double mu) {

    auto jacobian = eq_constraints.Jacobian(x);
    auto jacobian_ineq = ineq_constraints.Jacobian(x);
    auto g = obj.grad(x);
    double grad_error = (g-jacobian.transpose()*lambda-jacobian_ineq.transpose()*z).norm();

    double sz_minus_m_error = 0;
    for (int i = 0; i < s.rows(); i++) sz_minus_m_error += sqrt(pow(s(i)*z(i)-mu,2));
    double eq_const_error = eq_constraints.Vals(x).norm();
    double ineq_const_error = (ineq_constraints.Vals(x)-s).norm();
    if (mu == 0) {
        if (ineq_constraints.Vals(x).minCoeff() < 0) {ineq_const_error = numeric_limits<double>::infinity();}
    }
    std::cout << "\tgrad_error = " << grad_error << std::endl;
    std::cout << "\tsz_minus_m_error = " << sz_minus_m_error << std::endl;
    std::cout << "\teq_const_error = " << eq_const_error << std::endl;
    std::cout << "\tineq_const_error = " << ineq_const_error << std::endl;
    cout << "ineq_constraints.Vals(x).minCoeff() = " << ineq_constraints.Vals(x).minCoeff() << endl;
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
        int kkt_size = hessian_IJV.size() + const_lambda_hessian.size() + 2*jacobian_IJV.size() + var_n + const_n + ineq_const_n;
        kkt_size += ineq_lambda_hessian.size() + 2*ineq_jacobian_ijv.size() + s.size() + 2*z.size();
        kkt_IJV.resize(kkt_size);
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
        //kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(ineq_lambda_hessian[i].row(),ineq_lambda_hessian[i].col(),0);
    }
    for (int i = 0; i < var_n; i++) { kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,-eps);}

    // Add both J and J transpose
    int s_rows = s.rows();
    for (int i = 0; i < jacobian_IJV.size(); i++) {
        int row = jacobian_IJV[i].row(), col = jacobian_IJV[i].col(); double val = jacobian_IJV[i].value();
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(col,row+var_n+s_rows, val); // Jt columed offseted at var_n + s_rows
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(row+var_n+s_rows, col, val); // J offseted in var_n
    }
    // Add inequalities jacobian and transpose
    int lambda_rows = lambda.rows();
    for (int i = 0; i < ineq_jacobian_ijv.size(); i++) {
        int row = ineq_jacobian_ijv[i].row(), col = ineq_jacobian_ijv[i].col(); double val = ineq_jacobian_ijv[i].value();
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(col,row+var_n+s_rows+lambda_rows, val); // Jt columed offseted at var_n + s_rows + lambda_rows
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(row+var_n+s_rows+lambda_rows, col, val); // J offseted in var_n
    }

    // add diagonal matrix with -1*S^-1Z
    for (int i = 0; i < s.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(var_n+i,var_n+i, -z(i)/s(i));
    }
    // add minus Identity at size z
    for (int i = 0; i < s.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(var_n+i,var_n+s_rows+lambda_rows+i, -1);
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(var_n+s_rows+lambda_rows+i,var_n+i, -1);
    }

    // Add zeros along the rest of the diagonal
    // This is important since Pardiso requires the diagonal to be explicitly included in the matrix, even if it's zero
    for (int i = var_n+s_rows; i < var_n+s_rows+const_n+ineq_const_n; i++) kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,0);
    if ( A.rows() == 0) {
      A =  Eigen::SparseMatrix<double>(var_n+const_n+2*s_rows,var_n+const_n+2*s_rows);
      cout << "precomuting sparse matrix" << endl;
      igl::sparse_cached_precompute(kkt_IJV, cached_ijv_data, A);
    } else {
      cout << "building sparse matrix from cache" << endl;
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

void FeasibleIneqInteriorPoint::compute_step_direction(const Eigen::VectorXd& x, Objective& f, Constraints& eq_constraints, Constraints& ineq_constraints, 
        double mu, Eigen::VectorXd& step_d, double current_merit) {
    cout << "computing_step_direction" << std::endl;
    int vnum = x.rows()/3;
    double new_e;

    //std::cout << "lagrange multipliers norm = " << lambda.norm() << std::endl;

    igl::Timer timer; auto init_time = timer.getElapsedTime(); auto t = init_time;
    // Get Hessian
    //auto hessian = f.hessian(x); 
    cout << "getting hessian info" << endl;
    auto hessian_ijv = f.update_and_get_hessian_ijv(x);
    cout << "eq hessian info" << endl;
    auto lambda_hessian_ijv = eq_constraints.update_and_get_lambda_hessian(x,lambda);
    cout << "ineq hessian info" << endl;
    auto lambda_ineq_hessian_ijv = ineq_constraints.update_and_get_lambda_hessian(x,z);
    auto hessian_time = timer.getElapsedTime()-t;

    // Get Jacobian
    t = timer.getElapsedTime();
    cout << "getting jacobian info" << endl;
    auto jacobian = eq_constraints.Jacobian(x);
    auto jacobian_ineq = ineq_constraints.Jacobian(x);
    auto jacobian_ijv = eq_constraints.update_and_get_jacobian_ijv(x);
    auto ineq_jacobian_ijv = ineq_constraints.update_and_get_jacobian_ijv(x);
    auto jacobian_time = timer.getElapsedTime()-t;

    cout << "building kkt system" << endl;
    t = timer.getElapsedTime();
    //Eigen::SparseMatrix<double> A; build_kkt_system(hessian,jacobian,A);
    int var_n = x.rows(); int const_n = eq_constraints.getConstNum(); int ineq_const_n = ineq_constraints.getConstNum();
    build_kkt_system_from_ijv(hessian_ijv, lambda_hessian_ijv, var_n,
                        jacobian_ijv, const_n, lambda_ineq_hessian_ijv, ineq_jacobian_ijv, ineq_const_n);

    auto kkt_time = timer.getElapsedTime()-t;  

    t = timer.getElapsedTime();

    //energy->check_grad(x);
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    //cout << "analayzing pattern" << endl;
    //solver.analyzePattern(A);
    if (first_solve) {
        cout << "analyzing pattern" << endl;
        m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
        m_solver.set_pattern();
        m_solver.analyze_pattern();
        m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
        m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
        //m_solver.iparm[7]  = 5;       /* Max numbers of iterative refinement steps. */
        first_solve = false;
    } else {
        cout << "updating system matrix" << endl;
        m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    auto analyze_pattern_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    cout << "factorizing" << endl;
    m_solver.factorize();
    /*solver.factorize(A);
    if(solver.info()!=Eigen::Success) {
        cout << "Eigen Failure!" << endl;
        exit(1);
    }*/
    auto factorize_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();
    
    //m_solver.factorize();

    cout << "building rhs" << endl;
    // build rhs as the negative of the rhs Nocedal 19.12
    Eigen::VectorXd g = f.grad(x);
    Eigen::VectorXd rhs_upper = g-jacobian.transpose()*lambda-jacobian_ineq.transpose()*z;
    Eigen::VectorXd rhs_upper_below(s.rows()); for (int i = 0; i < s.rows(); i++) rhs_upper_below(i) = z(i)-mu/s(i);
    Eigen::VectorXd constraints_deviation = -1*eq_constraints.Vals(x);
    Eigen::VectorXd ineq_constraints_deviation = s-1*ineq_constraints.Vals(x);
    Eigen::VectorXd rhs(rhs_upper.rows()+rhs_upper_below.rows()+constraints_deviation.rows()+ineq_constraints_deviation.rows());
    rhs << rhs_upper,rhs_upper_below,constraints_deviation,ineq_constraints_deviation;
    
    cout << "solving" << endl;
    m_solver.solve(rhs,step_d);
    cout << "solved" << endl;
    cout << "result precision: (A*step_d-rhs).norm() = " << (A*step_d-rhs).norm() << endl;

    cout << "checking results" << endl << endl;
    Eigen::VectorXd px(x.rows()), ps(s.rows()),py(lambda.rows()),pz(z.rows());
    for (int i = 0; i < x.rows(); i++) px(i) = step_d(i);
    for (int i = 0; i < s.rows(); i++) ps(i) = step_d(x.rows()+i);
    for (int i = 0; i < lambda.rows(); i++) py(i) = step_d(x.rows()+s.rows()+i);
    for (int i = 0; i < z.rows(); i++) pz(i) = step_d(x.rows()+s.rows()+lambda.rows()+i);
    cout << "ps = " << ps << endl;
    int wait; cin >> wait;
    //cout << "ineq_constraints_deviation = " << endl << ineq_constraints_deviation << endl;
    //cout << "rhs = " << rhs << endl;
    //cout << "last equation norm: ((jacobian_ineq*px-ps) - (ineq_constraints_deviation)) = " << ((jacobian_ineq*px-ps) - (ineq_constraints_deviation)) << endl;
    //cout << "last equation norm: ((jacobian_ineq*px-ps) - (ineq_constraints_deviation)).norm() = " << ((jacobian_ineq*px-ps) - (ineq_constraints_deviation)).norm() << endl;
    //Eigen::VectorXd BigEps(s.rows()); for (int i = 0; i < s.rows(); i++) BigEps(i) = z(i)/s(i);
    //cout << "equation norm: ((BigEps.asDiagonal()*ps+pz) - (-1*rhs_upper_below)).norm() = " << ((BigEps.asDiagonal()*ps+pz) - (-1*rhs_upper_below)).norm() << endl;

    //exit(1);

    //g_const = -1*g_const;
    
    //cout << "solving!" << endl;
    //res = solver.solve(g_const);
    
    /*
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
    
    //return new_e;

}