#include "DOGFlowAndProject.h"

#include "../DogLaplacian.h"

DOGFlowAndProject::DOGFlowAndProject(const Dog& dog, int max_iter): dog_init(dog), max_iter(max_iter) {
		first_solve = true;
}

double DOGFlowAndProject::solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
	x = x0;
	// TODO: add stopping criteria
	double f;
	for (int iter = 0; iter < max_iter; iter++) {
		f = solve_single_iter(x, obj, constraints, x);
	}
	return f;
}

double DOGFlowAndProject::solve_single_iter(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
	flow(x0, obj, constraints, x);
	project(x, obj, constraints, x);
	return obj.obj(x);
}
void DOGFlowAndProject::flow(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {
	/*
	x = x0;
	int vnum = x.rows()/3;
	double new_e;

	igl::Timer timer; double time;

	// Get metric
	Eigen::MatrixXd V_x; vec_to_mat2(x,V_x); 
	auto metric = DOG_laplacian(V_x,dog_init.F);
	Eigen::SparseMatrix<double> id(3*V.rows(),3*V.rows()); id.setIdentity();
	metric = metric - (1e-8)*id;
	
	//energy->check_grad(x);
	double old_e = f.obj(x);

	double g_time = timer.getElapsedTime();
	//cout << "old e = " << old_e << endl;

	Eigen::VectorXd g(f.grad(x));
	
	//cout << "g.rows() = " << g.rows() << endl;
	Eigen::VectorXd d(g.rows());
	Eigen::SparseMatrix<double> J = constraints.Jacobian(x);
	
	Eigen::SparseMatrix<double> Jt = J.transpose();
	Eigen::SparseMatrix<double> L_jt; igl::cat(2,Metric,Jt, L_jt);			

	Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
	Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
	Eigen::SparseMatrix<double> A; igl::cat(1, L_jt, J_0, A);
	
	A.makeCompressed();
	Eigen::SparseMatrix<double> id_all(A.rows(),A.rows()); id_all.setIdentity();
	A = A + 0*id_all; // todo: stupid but I want to add zeros explicitly
	
    if (first_solve) {
    	time = timer.getElapsedTime();
		m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
	    m_solver.set_pattern();
    	m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices
    	m_solver.iparm[12] = 2; // imporved accuracy for highly indefinite symmetric matrices
    	m_solver.iparm[20] = 1;
	    m_solver.analyze_pattern();
	    first_solve = false;
    } else {
		m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    m_solver.factorize();

	Eigen::VectorXd zeroV(J.rows()); zeroV.setZero();
	Eigen::VectorXd constraints_deviation = -1*constraints.Vals();
	Eigen::VectorXd g_const; igl::cat(1, g, constraints_deviation, g_const);
	
	Eigen::VectorXd res;
	//cout << "solving!" << endl;
	m_solver.solve(g_const,res);

	for (int d_i = 0; d_i < g.rows(); d_i++) {
		d[d_i] = res[d_i];
	}
	new_e = line_search(x,d,flow_t,f);
	
	old_e = f.obj(x);
	
	return new_e;
	*/
}
void DOGFlowAndProject::project(const Eigen::VectorXd& x0, Objective& f, const Constraints& constraints, Eigen::VectorXd& x) {

}