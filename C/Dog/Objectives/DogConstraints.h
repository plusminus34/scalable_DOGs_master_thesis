#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

class DogConstraints: public Constraints {
public:
	DogConstraints(const QuadTopology& quadTop, bool offset_planar = true) : quadTop(quadTop), offset_planar(offset_planar)
					{const_n= 3*(quadTop.stars.rows()/5)+1*quadTop.bnd3.rows()/4;}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual std::vector<Eigen::Triplet<double> > JacobianIJV(const Eigen::VectorXd& x) const;

	Eigen::SparseMatrix<double> LambdaHessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) const;
private:
	const bool offset_planar;
	const QuadTopology& quadTop;
};