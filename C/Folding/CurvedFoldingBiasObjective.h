#pragma once

#include "../Dog/Dog.h"
#include "../Optimization/Objective.h"

struct CurvedFoldBias {
	EdgePoint ep_b, ep_f;
	// Assumes v1,w1 refers to the "vertices duplicated from the same vertex" such that the edge point ep_0 = edge_t*v1+(1-edge_t)*v2 = edge_t*w1+(1-edge_t)*w2
	int v1,v2; int w1,w2;
	double edge_t;
};

// comparing to squared length
class CurvedFoldingBiasObjective: public Objective {
public:
	CurvedFoldingBiasObjective(bool sign = false, double alpha = 100, bool dbg = false) : sign(sign), alpha(alpha), dbg(dbg) {}
	void add_fold_bias(const CurvedFoldBias& foldBias);
	void reset_folds();
	void set_use_sign(bool i_sign) {sign = i_sign;}
	
	virtual CurvedFoldingBiasObjective* clone() const {return new CurvedFoldingBiasObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);
	std::vector<CurvedFoldBias> curvedFoldBiases;

	virtual double obj_l2(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad_l2(const Eigen::VectorXd& x) const;
	virtual void updateHessianIJV_l2(const Eigen::VectorXd& x);

	virtual double obj_sign(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad_sign(const Eigen::VectorXd& x) const;
	virtual void updateHessianIJV_sign(const Eigen::VectorXd& x);

	bool sign;
	double alpha;
	bool dbg;
};