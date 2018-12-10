#pragma once

#include "../Solver.h"
/*
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

class IpOptSolver : public ConstrainedSolver {
  
public:
	//IpOptSolver(): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
    virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, 
            Eigen::VectorXd& x);
    
};

class IpOptVariables : public VariableSet {
public:
  IpOptVariables(Eigen::VectorXd& x_in) : ExVariables("IpOptVariables"), x(x_in) {};
  void SetVariables(const VectorXd& x_in) override {x = x_in;};

  VectorXd GetValues() const override {return x;};
private:
  Eigen::VectorXd x;
};

class IpOptObjective: public CostTerm {
public:
  IpOptObjective(Objective& obj_i) : CostTerm("IpOptCost"), obj(obj_i) {}

  double GetCost() const override {
    VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
    return obj.objective(x);
  };

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "IpOptVariables") {
      VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
      VectorXd g = obj.grad(x);
      jac.block(0,0,1,g.rows()) = g.transpose(); // should be one row apparently (so constant row idx)
    }
  }
private:
     Objective& obj;
};

class IpOptConstraints : public ConstraintSet {
public:
  IpOptConstraints(Constraints& constraints_i) : ConstraintSet("IpOptConstraints"), constraints(constraints_i) {}
  // The constraint value minus the constant value "1", moved to bounds.
  VectorXd GetValues() const override{
    VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
    return constraints->Vals(x);
  };
  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {
    if (var_set == "IpOptVariables") {
      VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
      Jacobian = constraints->Jacobian(x);
    }
  }
private:
    Constraints& constraints;
};
*/