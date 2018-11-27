#include "DOGFlowAndProject.h"

#include "../DogLaplacian.h"
#include "../Objectives/LaplacianSimilarity.h"

#include "igl/cat.h"

using namespace std;

DOGFlowAndProject::DOGFlowAndProject(const Dog& dog, double flow_t, const int& max_flow_project_iter, const int& max_lbfgs_proj_iter, 
								const int& penalty_repetitions): dog_init(dog), flow_t(flow_t), max_flow_project_iter(max_flow_project_iter),
									 max_lbfgs_proj_iter(max_lbfgs_proj_iter), penalty_repetitions(penalty_repetitions),
									 /*m_solver(ai,aj,K),*/ lbfgsWithPenalty(max_lbfgs_proj_iter,penalty_repetitions) {
	first_solve = true;
	//m_solver.set_type(-2);
}

double DOGFlowAndProject::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
	x = x0;
	// TODO: add stopping criteria
	double f;
	for (int iter = 0; iter < max_flow_project_iter; iter++) {
		f = solve_single_iter(x, obj, constraints, x);
	}
	return f;
}

double DOGFlowAndProject::solve_single_iter(const Eigen::VectorXd& x0, Objective& f, const Constraints& constraints, Eigen::VectorXd& x) {
	cout << "obj before flow = " << f.obj(x0) << endl;
	cout << "const deviation before flow = " << constraints.deviation(x0) << endl;
	auto e_after_flow = flow(x0, f, constraints, x);
	cout << "obj after flow = " << f.obj(x) << endl;
	//return e_after_flow;
	cout << "const deviation after flow = " << constraints.deviation(x) << endl;
	project(x, f, constraints, x);
	auto e_after_proj = f.obj(x);
	cout << "obj after project = " << f.obj(x) << endl;
	cout << "const deviation after project = " << constraints.deviation(x) << endl;
	return e_after_proj;
	
}
double DOGFlowAndProject::flow(const Eigen::VectorXd& x0, Objective& f, const Constraints& constraints, Eigen::VectorXd& x) {
	x = x0;
	int vnum = x.rows()/3;
	double new_e;

	// Get metric
	Eigen::MatrixXd V_x; vec_to_mat2(x,V_x); 
	auto metric = DOG_laplacian(V_x,dog_init.getF());
	Eigen::SparseMatrix<double> id(x.rows(),x.rows()); id.setIdentity();
	metric = metric - (1e-8)*id;
	metric = metric - f.hessian(x0); // adding hessian
	
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
	
    if (first_solve) {
    	/*
		m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
	    m_solver.set_pattern();
    	m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
    	m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
    	m_solver.iparm[20] = 1;
	    m_solver.analyze_pattern();
	    */
	    first_solve = false;
    } else {
		//m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    //m_solver.factorize();

	Eigen::VectorXd zeroV(J.rows()); zeroV.setZero();
	Eigen::VectorXd constraints_deviation = -1*constraints.Vals(x);
	Eigen::VectorXd g_const; igl::cat(1, g, constraints_deviation, g_const);
	
	Eigen::VectorXd res;
	//cout << "solving!" << endl;
	//m_solver.solve(g_const,res);

	for (int d_i = 0; d_i < g.rows(); d_i++) {
		d[d_i] = res[d_i];
	}
	new_e = line_search(x,d,flow_t,f);
	
	old_e = f.obj(x);
	
	return new_e;
}
void DOGFlowAndProject::project(const Eigen::VectorXd& x0, Objective& /*f*/, const Constraints& constraints, Eigen::VectorXd& x) {
	// We don't optimize for the same objective, but just project the given mesh to the closest DOG
	//	and we define "closest" in terms of the laplacian normals
	// Use a smoothness objective on DOG, together with the same constraints
	LaplacianSimilarity laplacianNormalsObj(dog_init, x0);
	lbfgsWithPenalty.solve_constrained(x0, laplacianNormalsObj, constraints, x);
	//lbfgsWithPenalty.solve_single_iter_with_fixed_p(x0, laplacianNormalsObj, constraints, x);
}

double DOGFlowAndProject::line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, double cur_energy) {
	double old_energy;
  if (cur_energy >= 0)
  {
    old_energy = cur_energy;
  }
  else
  {
    old_energy = f.obj(x); // no energy was given -> need to compute the current energy
  }
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = 22;

  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::VectorXd new_x = x + step_size * d;

    double cur_e = f.obj(new_x);
    
    //cout << "cur_e = " << cur_e << endl;
    if (cur_e >= old_energy)
    {
      step_size /= 2;
      //cout << "step_size = " << step_size << endl;
    }
    else
    {
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