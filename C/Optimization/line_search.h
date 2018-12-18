#pragma once

#include "Objective.h"
#include "Constraints.h"

double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, double cur_energy = -1);

double exact_l2_merit_lineserach(Eigen::VectorXd& x, const Eigen::VectorXd& d, double step_size, Objective& f, Constraints& constraints,
				const double& merit_penalty,
				double cur_energy = -1);


class ExactL2MeritObjective: public Objective {
  
public:
	ExactL2MeritObjective(const Objective& obj, const Constraints& constraints, const double& merit_p)  : inner_obj(obj), inner_constraints(constraints), merit_p(merit_p) {};
	virtual ExactL2MeritObjective* clone() const {return new ExactL2MeritObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const {
		return inner_obj.obj(x)+merit_p*inner_constraints.Vals(x).norm();
	};
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// Should not get here..
		Eigen::VectorXd g(x.rows());
		std::cout << "Error, L2Merit function is not differentiable" << std::endl; exit(1);
		return g;
	}

private:
	const Objective& inner_obj; 
	const Constraints& inner_constraints;
	const double& merit_p;
};