#pragma once

#include "Objective.h"
#include <igl/Timer.h>

class CompositeObjective: public Objective {
  
public:
	CompositeObjective(){/*empty on purpose*/}
	CompositeObjective(const std::vector<Objective*>& objectives_i, std::vector<double> weights) : weights(weights) {
		objectives.resize(objectives_i.size());
		for (int i = 0; i < objectives.size(); i++) objectives[i] = objectives_i[i]->clone();
		use_hessian.resize(objectives.size()); for (auto b : use_hessian) b = false;
	}
	virtual CompositeObjective* clone() const {return new CompositeObjective(*this);}
	void add_objective(Objective* e, double w = 1., bool c_use_hessian = false) {objectives.push_back(e->clone()); 
								weights.push_back(w); use_hessian.push_back(c_use_hessian);}

	void add_objective_permanent(Objective& e, double w = 1., bool c_use_hessian = false) {objectives.push_back(&e); 
								weights.push_back(w); use_hessian.push_back(c_use_hessian);}

	virtual double obj(const Eigen::VectorXd& x) const {
			double obj = 0;
			for (int i = 0; i < objectives.size(); i++) {obj+=weights[i]*objectives[i]->obj(x);}
			return obj;
	}
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		Eigen::VectorXd grad(x.rows()); grad.setZero();
		for (int i = 0; i < objectives.size(); i++) {grad+=weights[i]*objectives[i]->grad(x);}
		return grad;
	};

	virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x) {
		igl::Timer timer; double init_t = timer.getElapsedTime();
		cachedH = Eigen::SparseMatrix<double>(x.rows(),x.rows());
		for (int i = 0; i < objectives.size(); i++) {
			if (use_hessian[i]) {
				cachedH+=weights[i]*objectives[i]->hessian(x);
			};
		}
		return cachedH;
	};

	virtual void set_ref(const Eigen::VectorXd& x0) {for (auto obj : objectives) {obj->set_ref(x0);}};

private:
	std::vector<Objective*> objectives;
	std::vector<double> weights;
	std::vector<bool> use_hessian;
};

//This hessian time = 0.000199551
//This hessian time = 0.00027908