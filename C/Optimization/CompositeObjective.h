#pragma once

#include "Objective.h"
#include <igl/Timer.h>

class CompositeObjective: public Objective {
  
public:
	CompositeObjective(){/*empty on purpose*/}
	CompositeObjective(const std::vector<Objective*>& objectives_i, const std::vector<double>& weights) : weights(weights) {
		objectives.resize(objectives_i.size());
		for (int i = 0; i < objectives.size(); i++) {
			objectives[i] = objectives_i[i]->clone();
			ijv_size += objectives[i]->get_hessian_IJV_size();
		}
		IJV.resize(ijv_size);
	}

	void update_weights(const std::vector<double>& weights_i) { weights = weights_i;}

	virtual CompositeObjective* clone() const {return new CompositeObjective(*this);}
	void add_objective(Objective* e, double w = 1.) {
		objectives.push_back(e->clone()); 
		weights.push_back(w);
		ijv_size += e->get_hessian_IJV_size();
		IJV.resize(ijv_size);
	}

	void add_objective_permanent(Objective& e, double w = 1.) {
		objectives.push_back(&e); 
		weights.push_back(w);
		ijv_size += e.get_hessian_IJV_size();
		IJV.resize(ijv_size);
	}

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

	virtual void set_ref(const Eigen::VectorXd& x0) {for (auto obj : objectives) {obj->set_ref(x0);}};

		virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x) {
		igl::Timer timer; double init_t = timer.getElapsedTime();
		cachedH = Eigen::SparseMatrix<double>(x.rows(),x.rows());
		for (int i = 0; i < objectives.size(); i++) {
			cachedH+=weights[i]*objectives[i]->hessian(x);
		}
		return cachedH;
	};


private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x) {
		int ijv_idx = 0;
		for (int i = 0; i < objectives.size(); i++) {
			const std::vector<Eigen::Triplet<double> >& obj_IJV = objectives[i]->update_and_get_hessian_ijv(x);
			for (auto val : obj_IJV) {IJV[ijv_idx++] = Eigen::Triplet<double>(val.row(),val.col(),weights[i]*val.value());}
		}
	}

	std::vector<Objective*> objectives;
	std::vector<double> weights;
	int ijv_size = 0;
};

//This hessian time = 0.000199551
//This hessian time = 0.00027908