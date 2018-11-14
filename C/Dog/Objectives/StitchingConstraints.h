#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class StitchingConstraints: public Constraints {
public:
	StitchingConstraints(const QuadTopology& quadTop,const DogEdgeStitching& edgeStitching) : quadTop(quadTop), 
																					eS(edgeStitching)
								{const_n= 3*edgeStitching.edge_coordinates.size();
								approx_nnz = 3*4*const_n; /*2 vertices for equality const*/}
	virtual StitchingConstraints* clone() const {return new StitchingConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const;
private:
	const QuadTopology& quadTop;
	const DogEdgeStitching& eS;
};
