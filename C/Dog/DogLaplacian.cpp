#include "DogLaplacian.h"

#include "igl/parallel_for.h"
#include "igl/repdiag.h"

using namespace Eigen;

Eigen::SparseMatrix<double> DOG_laplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr) {
	Eigen::MatrixXd squaredL(Fsqr.rows(),4);
	Eigen::MatrixXi edges(4,2);
	edges <<
		0,1,
		1,2,
		2,3,
		3,0;
	igl::parallel_for(
        Fsqr.rows(),
        [&V,&Fsqr,&squaredL](const int i)
        {
          squaredL(i,0) = (V.row(Fsqr(i,0))-V.row(Fsqr(i,1))).squaredNorm();
          squaredL(i,1) = (V.row(Fsqr(i,1))-V.row(Fsqr(i,2))).squaredNorm();
          squaredL(i,2) = (V.row(Fsqr(i,2))-V.row(Fsqr(i,3))).squaredNorm();
          squaredL(i,3) = (V.row(Fsqr(i,3))-V.row(Fsqr(i,0))).squaredNorm();
        },
        1000);
	Eigen::MatrixXd e_length = squaredL.array().sqrt();

	Eigen::SparseMatrix<double> L(V.rows(),V.rows());
	// Inner vertices have 4 neighbours (so 5 per rows including themselves), other vertices have less
	L.reserve(5*V.rows());
	std::vector<Eigen::Triplet<double> > IJV;
	// Loop over squares
	for(int i = 0; i < Fsqr.rows(); i++) {
		// loop over edges of element
		for(int e = 0;e<edges.rows();e++) {
		  int source = Fsqr(i,edges(e,0));
		  int dest = Fsqr(i,edges(e,1));
		  double l_e = e_length(i,e);
		  double l_nb1 = e_length((i+1)%4,e), l_nb2 = e_length((i+3)%4,e);
		  double l_ratio = 1./16* ((l_nb1+l_nb2)/l_e);
		  IJV.push_back(Triplet<double>(source,dest,l_ratio));
		  IJV.push_back(Triplet<double>(dest,source,l_ratio));
		  IJV.push_back(Triplet<double>(source,source,-l_ratio));
		  IJV.push_back(Triplet<double>(dest,dest,-l_ratio));
		}
	}
	L.setFromTriplets(IJV.begin(),IJV.end());

	Eigen::SparseMatrix<double> L3d;
	igl::repdiag(L,3,L3d); // one way
	//cout << "get laplacian matrix with norm " << L.norm() << endl;
	//L = build_uniform_laplacian_mat_3d(quadTop);
	return L3d;
}

Eigen::SparseMatrix<double> uniform_laplacian(const QuadTopology& quadTop) {
  int n = quadTop.v_n;
  int entries_num = quadTop.stars.rows() + quadTop.bnd4.rows() + quadTop.bnd3.rows() + quadTop.bnd2.rows();
  std::vector<Eigen::Triplet<double> > IJV;
  IJV.reserve(entries_num);

  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0 = quadTop.stars(si), p_xf = quadTop.stars(si+1), p_yf = quadTop.stars(si+2), p_xb = quadTop.stars(si+3),p_yb = quadTop.stars(si+4);
    IJV.push_back(Triplet<double>(p_0,p_0,-1));
    IJV.push_back(Triplet<double>(p_0,p_xf,0.25));
    IJV.push_back(Triplet<double>(p_0,p_xb,0.25));
    IJV.push_back(Triplet<double>(p_0,p_yf,0.25));
    IJV.push_back(Triplet<double>(p_0,p_yb,0.25));
  }
  
  for (int si = 0; si < quadTop.bnd4.rows(); si+=5) {
    int p_0 = quadTop.stars(si), p_xf = quadTop.stars(si+1), p_yf = quadTop.stars(si+2), p_xb = quadTop.stars(si+3),p_yb = quadTop.stars(si+4);
    IJV.push_back(Triplet<double>(p_0,p_0,-1));
    IJV.push_back(Triplet<double>(p_0,p_xf,0.25));
    IJV.push_back(Triplet<double>(p_0,p_xb,0.25));
    IJV.push_back(Triplet<double>(p_0,p_yf,0.25));
    IJV.push_back(Triplet<double>(p_0,p_yb,0.25));
  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0 = quadTop.bnd3(si), p_xf = quadTop.bnd3(si+1), p_yf = quadTop.bnd3(si+2), p_yb = quadTop.bnd3(si+3);
    IJV.push_back(Triplet<double>(p_0,p_0,-1));
    IJV.push_back(Triplet<double>(p_0,p_xf,1./3));
    IJV.push_back(Triplet<double>(p_0,p_yf,1./3));
    IJV.push_back(Triplet<double>(p_0,p_yb,1./3));
  }

  for (int si = 0; si < quadTop.bnd2.rows(); si+=3) {
    int p_0 = quadTop.bnd2(si), p_xf = quadTop.bnd2(si+1), p_yf = quadTop.bnd2(si+2);
    IJV.push_back(Triplet<double>(p_0,p_0,-1));
    IJV.push_back(Triplet<double>(p_0,p_xf,0.5));
    IJV.push_back(Triplet<double>(p_0,p_yf,0.5));
  }
  // build matrix
  Eigen::SparseMatrix<double> L; L.resize(n,n);
  L.setFromTriplets(IJV.begin(),IJV.end());

  Eigen::SparseMatrix<double> L_3d;
  igl::repdiag(L,3,L_3d);
  return L_3d;
}

Eigen::VectorXd DOG_vertex_area(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr) {
	Eigen::VectorXd vArea(V.rows()); vArea.setZero();
	Eigen::MatrixXi edges(4,2);
	edges <<
		0,1,
		1,2,
		2,3,
		3,0;
	igl::parallel_for(
        Fsqr.rows(),
        [&V,&Fsqr,&vArea](const int i)
        {

          double l1 = (V.row(Fsqr(i,0))-V.row(Fsqr(i,1))).squaredNorm();
          double l2 = (V.row(Fsqr(i,1))-V.row(Fsqr(i,2))).squaredNorm();
          double l3 = (V.row(Fsqr(i,2))-V.row(Fsqr(i,3))).squaredNorm();
          double l4 = (V.row(Fsqr(i,3))-V.row(Fsqr(i,0))).squaredNorm();
          double sqr_area = 0.25*(l1+l3)*(l2+l4);
          // Add a quarter of the area to each one of the vertices
          for (int j = 0; j < 4; j++) vArea(Fsqr(i,j))+= 0.25*(sqr_area);
        },
        1000);

	// now go through these
	return vArea;
}