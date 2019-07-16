#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

// This objective is used to make cylindrical-like boundary glueing: 
// 	when we bend a planar piece into a cylinder we don't only want to fit two boundaries but also want to make it smooth along the fitted curve
//	this can be achieved by adding a bending objective along it
class PairedBoundaryVerticesBendingObjective: public Objective {
  
public:
	PairedBoundaryVerticesBendingObjective(const QuadTopology& quadTop, const std::vector<std::pair<int,int>>& pairs, 
			const Eigen::VectorXd& x0_init, const Eigen::Vector3d& axis_direction = Eigen::Vector3d(1,0,0));
	virtual PairedBoundaryVerticesBendingObjective* clone() const {return new PairedBoundaryVerticesBendingObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);
	// Given a boundary index v with only
	int find_neighbour_in_axis_direction(int v, const QuadTopology& quadTop, const Eigen::VectorXd& x, const Eigen::Vector3d& axis_direction);
	
	const QuadTopology& quadTop;
	int vnum;
	std::vector<double> init_edge_lengths;
	//const std::vector<std::pair<int,int>> pairs;
	// Each pair has a boundary bending objective defined on 3 vertices, which we concatenate here
	Eigen::VectorXi obj_vertices;
};