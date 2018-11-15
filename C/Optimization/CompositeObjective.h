#pragma once

#include "Objective.h"

class CompositeObjective: public Objective {
  
public:
	CompositeObjective(const std::vector<Objective*>& objectives_i, std::vector<double> weights) : weights(weights) {
		objectives.resize(objectives_i.size());
		for (int i = 0; i < objectives.size(); i++) objectives[i] = objectives_i[i]->clone();
	}
	virtual CompositeObjective* clone() const {return new CompositeObjective(*this);}
	void add_objective(Objective* e, double w = 1.) {objectives.push_back(e->clone()); weights.push_back(w);}

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

	Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const {
		Eigen::SparseMatrix<double> hessian(x.rows(),x.rows());
		for (int i = 0; i < objectives.size(); i++) {hessian+=weights[i]*objectives[i]->hessian(x);}
		return hessian;
	};

	virtual void set_ref(const Eigen::VectorXd& x0) {for (auto obj : objectives) {obj->set_ref(x0);}};

private:
	std::vector<Objective*> objectives;
	std::vector<double> weights;
};