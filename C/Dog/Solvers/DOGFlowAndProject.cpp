#include "DOGFlowAndProject.h"

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

}
void DOGFlowAndProject::project(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, Eigen::VectorXd& x) {

}