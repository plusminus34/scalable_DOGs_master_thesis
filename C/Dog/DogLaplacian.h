#include "../QuadMesh/Quad.h"

Eigen::SparseMatrix<double> DOG_laplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr);
Eigen::SparseMatrix<double> uniform_laplacian(const QuadTopology& quadTop);
Eigen::VectorXd DOG_vertex_area(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr);