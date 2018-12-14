/*
#pragma once

#include "../Solver.h"

#include "/Users/michaelrabinovich/ifopt/ifopt_core/include/ifopt/variable_set.h"
#include "/Users/michaelrabinovich/ifopt/ifopt_core/include/ifopt/constraint_set.h"
#include "/Users/michaelrabinovich/ifopt/ifopt_core/include/ifopt/cost_term.h"
//#include "../../ifopt/variable_set.h"
//#include <ifopt/constraint_set.h>
//#include <ifopt/cost_term.h>

class IpOptSolver : public ConstrainedSolver {
  
public:
	//IpOptSolver(): max_iter(max_iter) {}
	// x0 is the initial guess, x is the result, the return value is the objective value
    virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, const Constraints& constraints, 
            Eigen::VectorXd& x);
    
};

class IpOptVariables : public ifopt::VariableSet {
public:
  IpOptVariables(const Eigen::VectorXd& x_in) : VariableSet(x_in.rows(), "IpOptVariables"), x(x_in) {};
  void SetVariables(const VectorXd& x_in) override {x = x_in;};

  VectorXd GetValues() const override {return x;};

    // Each variable has an upper and lower bound set here
  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    for (int i = 0; i < GetRows(); i++) {bounds.at(i) = ifopt::NoBound;}
    return bounds;
  }
private:
  Eigen::VectorXd x;
};

class IpOptObjective: public ifopt::CostTerm {
public:
  IpOptObjective(Objective& obj_i) : CostTerm("IpOptCost"), obj(obj_i) {}

  double GetCost() const override {
    VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
    return obj.obj(x);
  };

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override {
    if (var_set == "IpOptVariables") {
      VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
      VectorXd g = obj.grad(x);
      for (int i = 0; i < g.rows(); i++) {jac.coeffRef(0,i)=g(i);}
    }
  }
private:
     Objective& obj;
};

class IpOptConstraints : public ifopt::ConstraintSet {
public:
  IpOptConstraints(const Constraints& constraints_i) : ConstraintSet(constraints_i.getConstNum(),"IpOptConstraints"),
                                                     constraints(constraints_i) {}
  // The constraint value minus the constant value "1", moved to bounds.
   // TODO remove the bounds from here!

  VectorXd GetValues() const override{
    VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
    return constraints.Vals(x);
  };
  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override {
    if (var_set == "IpOptVariables") {
      VectorXd x = GetVariables()->GetComponent("IpOptVariables")->GetValues();
      jac_block = constraints.Jacobian(x);
    }
  }

  VecBound GetBounds() const override {
    VecBound bounds(GetRows());
    //TODO: implement get bounds in constraints!!
    for (int i = 0; i < GetRows(); i++) {bounds.at(i) = 0;}
    return bounds;
  }
private:
    const Constraints& constraints;
};

*/