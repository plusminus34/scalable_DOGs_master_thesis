#include "LBFGS.h"

double LBFGS::solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) {
	//param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;

	igl::Timer timer;
	LBFGS_obj_grad_interface fun(obj);

	//fun.check_grad(x0);
	double fx;
	//double time = timer.getElapsedTime();
	x = x0;
	int niter = solver->minimize(fun, x, fx);
	//cout << "Total solve took " << timer.getElapsedTime() - time << endl;
	//cout << "obj count = " << fun.obj_cnt << " grad count = " << fun.g_cnt << "both_cnt = " << fun.both_cnt << endl;
	//cout << "geo energy = " << compute_orth_dev_e(x0) << endl;

	return fx;
}