#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class FairingObjective: public Objective {
public:
	FairingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0);
	virtual FairingObjective* clone() const {return new FairingObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;
	virtual void set_ref(const Eigen::VectorXd& x0);

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);

	void add_indices_for_vertex(int p_0_i, int p_xb_i, int p_xf_i, const QuadTopology& quadTop, const Eigen::VectorXd& x0, int& const_n);
	
	// Indices of curve points to compare curvature normals
	std::vector<int> p_xb,p_0,p_xf,p_xff;
	std::vector<double> len_ex_b_v,len_ex_f_v,len_ex_ff_v;
};