#pragma once
#include <Eigen/Dense>
#include "igl/serialize.h"

#include <vector>

struct QuadTopology : public igl::Serializable {
    Eigen::VectorXi stars;
    Eigen::VectorXi bnd4;
    Eigen::VectorXi bnd3; 
    Eigen::VectorXi bnd2;

    std::vector<int> bnd_loop;
    std::vector<bool> is_bnd_v;
    Eigen::MatrixXi E; // quad edges
    std::vector<std::vector<int> > VF; // vertices to faces(quads)
    std::vector<std::vector<int> > A; // quad adjacency list
    std::vector<int> vi_to_star; // inner vertex to the star index in Eigen::VectorXi/ In case it's not an inner vertex, just -1

    int v_n;

    void InitSerialization() {
      Add(stars,std::string("_w_H"));
      Add(bnd4,std::string("_w_FirstFund"));
      Add(bnd3,std::string("_w_GridReg"));
      Add(bnd2,std::string("_w_SecondFund"));
      Add(bnd_loop,std::string("_w_smooth"));
      Add(is_bnd_v,std::string("_w_cylindrical"));
      Add(E,std::string("_w_conical"));
      Add(VF,std::string("_w_Href"));
      Add(A,std::string("_w_sqrGrid"));
      Add(vi_to_star,std::string("_w_vEq"));
    }
};

struct Edge  : public igl::Serializable {

  Edge(): v1(-1),v2(-1){}
  Edge(int v1_, int v2_): v1(v1_),v2(v2_) {
    //if (v1 > v2) { std::swap(v1,v2);}
  }
  Edge(const Edge& edge): v1(edge.v1),v2(edge.v2) {}

  void InitSerialization() {
    Add(v1,std::string("v1"));
    Add(v2,std::string("v2"));
  }
  int v1,v2;

  inline bool operator==(const Edge& rhs) const { /* do actual comparison */ 
    if (v1 == rhs.v1 && v2 == rhs.v2) {
      return true;
    }
    if (v1 == rhs.v2 && v2 == rhs.v1) {
      return true;
    }
    return false;
  }

  bool operator<( const Edge& rhs ) const {
    if (v1 < rhs.v1) return true;
    else if (v1 > rhs.v1) return false;
    else return v2 < rhs.v2;
  }
};

struct EdgeConstI {int const_i=-1; bool edge_1;};
struct FaceEdgeConst{EdgeConstI e_const1, e_const2;};

std::vector<bool> is_border_vertex(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr);

void quad_topology(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Fsqr, QuadTopology& quadTop);

void quad_adjacenecy_list(const Eigen::MatrixXi& Fsqr, std::vector<std::vector<int> >& A);
void vertex_quad_adjacency(const Eigen::MatrixXi& Fsqr, std::vector<std::vector<int> >& VF);
void quad_edges(const Eigen::MatrixXi & F, Eigen::MatrixXi& E);

void get_stars_and_bnd_vertices_nbds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& Fsqr, 
	Eigen::VectorXi& stars, std::vector<int>& vi_to_star, Eigen::VectorXi& bnd4, Eigen::VectorXi& bnd3, Eigen::VectorXi& bnd2);

// F here is triangular faces.. but sqr_f is squre face id
void get_edges_from_face(const Eigen::MatrixXi& F, int sqr_f, std::vector<Edge>& quad_edges);

std::vector<int> get_edge_faces(const QuadTopology& quadTop, const Eigen::RowVectorXi& e);
std::vector<int> get_edge_faces(const QuadTopology& quadTop, Edge& e);

void ivec_to_eigen(const std::vector<int>& vec, Eigen::VectorXi& ei_v);

void mat2_to_vec(const Eigen::MatrixXd& mat, Eigen::VectorXd& vec);
void vec_to_mat2(const Eigen::VectorXd& vec, Eigen::MatrixXd& mat);

void get_planar_square_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int square_h, int square_w);
void get_planar_square_mesh_V(Eigen::MatrixXd& V, int square_h, int square_w);
void get_planar_square_mesh_F(Eigen::MatrixXi& F, int square_h, int square_w);

Eigen::MatrixXi F_to_Fsqr(const Eigen::MatrixXi& F);
Eigen::MatrixXi Fsqr_to_F(const Eigen::MatrixXi& Fsqr);
