#include "Lasso.h"

#include <iostream>
#include <fstream>


#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/point_in_poly.h>
#include <igl/facet_components.h>
#include <igl/barycenter.h>

#include <igl/unproject_in_mesh.h>
#include <igl/Hit.h>

#include <unsupported/Eigen/Splines>

using namespace Eigen;
typedef Spline3d::PointType PointType;
typedef Spline3d::KnotVectorType KnotVectorType;
typedef Spline3d::ControlPointVectorType ControlPointVectorType;

using namespace igl;
using namespace std;

Lasso::Lasso(igl::opengl::glfw::Viewer& v,
             const Eigen::MatrixXd &V_,
             const Eigen::MatrixXi &F_tri):
V(V_),
F(F_tri),
viewer(v)
{
  //ei.init(V.cast<float>(),F);
}

Lasso::~Lasso()
{
}

void Lasso::reinit()
{
//  ei.deinit();
  //ei.init(V.cast<float>(),F);
  strokePoints.clear();
  strokeFacesSqr.clear();
  tri_hits.clear();
  const_edges.clear();
  edge_coordinates.clear();
  F1.clear();
  F2.clear();
  splinePoints.clear();
}

int Lasso::pickVertex(int mouse_x, int mouse_y)
{
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core().viewport(3) - mouse_y;

  Eigen::RowVector3d pt;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;
  int vi = -1;

  std::vector<igl::Hit> hits;
  /*
  igl::unproject_in_mesh(Eigen::Vector2f(x,y),
                         modelview,
                         viewer.core().proj,
                         viewer.core().viewport,
                         ei,pt,hits);
  */

  igl::unproject_in_mesh(Eigen::Vector2f(x,y), viewer.core().view,// * viewer.core().model,
      viewer.core().proj, viewer.core().viewport, V, F, pt,hits);

  if (hits.size()> 0) {
    int fi = hits[0].id;
    Eigen::RowVector3d bc;
    bc << 1.0-hits[0].u-hits[0].v, hits[0].u, hits[0].v;
    bc.maxCoeff(&vi);
    vi = F(fi,vi);
  }
  return vi;
}

Edge Lasso::pickEdge(int mouse_x, int mouse_y) {
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core().viewport(3) - mouse_y;

  Eigen::RowVector3d pt;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;
  int v1 = -1,v2=-1;

  std::vector<igl::Hit> hits;
  /*
  igl::unproject_in_mesh(Eigen::Vector2f(x,y),
                         modelview,
                         viewer.core().proj,
                         viewer.core().viewport,
                         ei,pt,hits);
  */

  igl::unproject_in_mesh(Eigen::Vector2f(x,y), viewer.core().view,// * viewer.core().model,
      viewer.core().proj, viewer.core().viewport, V, F, pt,hits);

  if (hits.size()> 0) {
    int fi = hits[0].id;
    Eigen::RowVector3d bc;
    bc << 1.0-hits[0].u-hits[0].v, hits[0].u, hits[0].v;
    bc.maxCoeff(&v1);
    bc(v1) = 0;
    v1 = F(fi,v1);
    bc.maxCoeff(&v2);
    v2 = F(fi,v2);
  }
  Edge edge;
  edge.v1 = v1; edge.v2 = v2;
  return edge;
}

int Lasso::strokeAdd(int mouse_x,
                    int mouse_y)
{
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core().viewport(3) - mouse_y;

  std::vector<unsigned> pt2D; pt2D.push_back(x); pt2D.push_back(y);
  stroke2DPoints.push_back(pt2D);

  Eigen::RowVector3d pt;
  int fi = -1;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;

  if (d<0)//first time
  {
    std::vector<igl::Hit> hits;
    /*
    igl::unproject_in_mesh(Eigen::Vector2f(x,y),
                           modelview,
                           viewer.core().proj,
                           viewer.core().viewport,
                           ei,pt,hits);
    */
      igl::unproject_in_mesh(Eigen::Vector2f(x,y), viewer.core().view,// * viewer.core().model,
      viewer.core().proj, viewer.core().viewport, V, F, pt,hits);
    if (hits.size()> 0)
    {
      fi = hits[0].id;
      Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core().proj,viewer.core().viewport);
      d = proj[2];
    }
  }

  // This is lazy, it will find more than just the first hit
  Eigen::Vector3f pt2 = igl::unproject(Eigen::Vector3f(x,y,0.95*d), modelview, viewer.core().proj, viewer.core().viewport);
  pt = pt2.transpose().cast<double>();


  strokePoints.push_back(pt);

  return fi;

}


int Lasso::strokeAddCurve(int mouse_x,
                    int mouse_y)
{
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core().viewport(3) - mouse_y;

  Eigen::RowVector3d pt;
  int fi = -1;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;

  std::vector<igl::Hit> hits;
    /*
    igl::unproject_in_mesh(Eigen::Vector2f(x,y),
                           modelview,
                           viewer.core().proj,
                           viewer.core().viewport,
                           ei,pt,hits);
    */
    igl::unproject_in_mesh(Eigen::Vector2f(x,y), viewer.core().view,// * viewer.core().model,
    viewer.core().proj, viewer.core().viewport, V, F, pt,hits);
  if (hits.size()> 0) {
    fi = hits[0].id;
    Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core().proj,viewer.core().viewport);
    d = proj[2];
    strokeFacesSqr.insert(int(fi/2)); // square face, so divide by two from triangle faces
    tri_hits.push_back(hits[0]);
    //selected_sqr_faces_to_tri_hit[strokeFacesSqr.size].push_back(hits[0]); // the index will be square one, but the hit structure will contain triangular one (this is useful)
    //selected_sqr_faces_to_tri_hit[fi/2].push_back(hits[0]); // the index will be square one, but the hit structure will contain triangular one (this is useful)
    strokePoints.push_back(pt);
  }
  return fi;
}

void Lasso::set_curve_manually(std::vector< Eigen::Matrix<double, 1,3>  >& _strokePoints, std::set<int>& _strokeFacesSqr, std::vector<igl::Hit>& _tri_hits) {
  strokePoints = _strokePoints; strokeFacesSqr = _strokeFacesSqr; tri_hits = _tri_hits;
}

void Lasso::set_curve_manually(Eigen::MatrixXd& _strokePoints, Eigen::VectorXi& faces_tri, Eigen::MatrixXd& BarycentricCoords) {
  // set strokepoints
  strokePoints.resize(_strokePoints.rows()); for (int i = 0 ; i < strokePoints.size(); i++) {strokePoints[i] = _strokePoints.row(i);}
  // set strokeFacesSqr
  strokeFacesSqr.clear(); for (int i = 0 ; i < strokePoints.size(); i++) {strokeFacesSqr.insert(faces_tri(i)/2); /*cout << "added faces_tri(i) = " << faces_tri(i) << endl;*/}
  // set tri hits
  tri_hits.clear();
  for (int i = 0; i < BarycentricCoords.rows(); i++) {
    igl::Hit hit; hit.id = faces_tri(i);
    hit.u = BarycentricCoords(i,1); hit.v = BarycentricCoords(i,2);
    tri_hits.push_back(hit);
  }
  cout << "tri_hits.size() = " << tri_hits.size() << endl;
}

void Lasso::get_hit_closest_edge_coordinates(igl::Hit hit, const std::vector<Edge> &edges, Edge& closest_edge, std::pair<double,double>& edge_coords) {
  int tri_f = hit.id;
  Eigen::RowVector3d bc;
  bc << 1.0-hit.u-hit.v, hit.u, hit.v;

  Edge e1(F(tri_f,0),F(tri_f,1)); Edge e2(F(tri_f,0),F(tri_f,2)); Edge e3(F(tri_f,1),F(tri_f,2));
  double e1_d = bc(0)+bc(1); double e2_d = bc(0)+bc(2); double e3_d = bc(1)+bc(2);

  //for (auto edge:edges) {cout << "edge: " << edge.v1 << "," << edge.v2 << endl;}
  //cout << "e1 = " << e1.v1 << "," << e1.v2 << endl;
  //cout << "e2 = " << e2.v1 << "," << e2.v2 << endl;
  //cout << "e3 = " << e3.v1 << "," << e3.v2 << endl;

  //cout << "Before: e1_d = " << e1_d << " e2_d = " << e2_d << " e3_d = " << e3_d << endl;
  // zero out the triangle edge
  if  (std::find(std::begin(edges), std::end(edges), e1) == std::end(edges) ) {e1_d = 0;}
  if  (std::find(std::begin(edges), std::end(edges), e2) == std::end(edges) ) {e2_d = 0;}
  if  (std::find(std::begin(edges), std::end(edges), e3) == std::end(edges) ) {e3_d = 0;}
  //cout << "After: e1_d = " << e1_d << " e2_d = " << e2_d << " e3_d = " << e3_d << endl;

  if ( (e1_d > e2_d) && (e1_d > e3_d) ) {
    //cout << "Winner:e1" << endl;
    closest_edge = e1; edge_coords.first = bc(0); edge_coords.second = bc(1);
  } else if ( (e2_d > e1_d) && (e2_d > e3_d) ) {
    //cout << "Winner:e2" << endl;
    closest_edge = e2; edge_coords.first = bc(0); edge_coords.second = bc(2);
  } else {
    //cout << "Winner:e3" << endl;
    closest_edge = e3; edge_coords.first = bc(1); edge_coords.second = bc(2);
  }
  //cout << "bc = " << bc << endl;
  //cout << "1) edge_coords.first = " << edge_coords.first << " edge_coords.second = " << edge_coords.second << endl;
  // normalize edge_coords
  double mult = 1./(edge_coords.first + edge_coords.second);
  edge_coords.first*= mult; edge_coords.second*= mult;
  //cout << "2) edge_coords.first = " << edge_coords.first << " edge_coords.second = " << edge_coords.second << endl;
}

void Lasso::get_two_face_neighbour_groups() {
  std::vector<int> nbd_quads; std::set<int> nbd_quads_with_selected_faces;
  Eigen::MatrixXi TT,TTi; igl::triangle_triangle_adjacency(F,TT,TTi);
  for (auto f : strokeFacesSqr) {
    int tri1 = 2*f; int tri2 = 2*f+1;
    for (int tti = 0; tti < 3; tti++) {if (TT(tri1,tti) != -1) nbd_quads_with_selected_faces.insert(TT(tri1,tti)/2);}
    for (int tti = 0; tti < 3; tti++) {if (TT(tri2,tti) != -1) nbd_quads_with_selected_faces.insert(TT(tri2,tti)/2);}
  }
  std::set_difference(nbd_quads_with_selected_faces.begin(), nbd_quads_with_selected_faces.end(), strokeFacesSqr.begin(), strokeFacesSqr.end(),
                        std::inserter(nbd_quads, nbd_quads.begin()));
  //for (auto f: nbd_quads) { cout << "nbd quad f = " << f << endl;}

  Eigen::MatrixXi FMinusSelected = F; // A bit of a hack: work with the triangular connectivity just so we could se fact_components
  for (int i = 0; i < F.rows(); i++) {
    if (strokeFacesSqr.find(i/2) != strokeFacesSqr.end()) {
      FMinusSelected.row(i) << -1,-1,-1;
    }
  }
  Eigen::VectorXi cid; igl::facet_components(FMinusSelected,cid);
  int first_comp = cid(2*nbd_quads[0]); //working with triangles..
  for (auto f : nbd_quads) {
    if (cid(2*f) == first_comp) {
      F1.push_back(f);
    } else {
      F2.push_back(f);
    }
  }
}

void Lasso::set_spline_points_from_stroke_points(int spline_pt_number) {
  splinePoints.clear();
  ControlPointVectorType points(3,spline_pt_number);// = ControlPointVectorType::Random(3,100);
  int point_jumps = strokePoints.size()/spline_pt_number;
  for (int i = 0; i < spline_pt_number-1; i++) {
    points.col(i) = strokePoints[point_jumps*i];
  }
  points.col(spline_pt_number-1) = strokePoints[strokePoints.size()-1]; //make sure we interpolate the lst point
  KnotVectorType chord_lengths; // knot parameters
  Eigen::ChordLengths(points, chord_lengths);

  //cout << "chord_lengths.cols() = " << chord_lengths.cols() << endl;

  const Spline3d spline = SplineFitting<Spline3d>::Interpolate(points,2,chord_lengths);
  //cout << "strokePoints.size() = " << strokePoints.size() << endl;
  double range = chord_lengths(chord_lengths.cols()-1)-chord_lengths(0);

  int s_num = 2000;
  splinePoints.push_back(strokePoints[0]);
  for(int i=1;i<s_num;i++) {
      double t = chord_lengths(0) + range*((double)i)/(s_num);
      PointType pt = spline(t);
      splinePoints.push_back(pt);
    }
    splinePoints.push_back(strokePoints[strokePoints.size()-1]);

    //cout << "and points = " << points << endl;

}

void Lasso::updateStrokeAndHitsFromSpline() {
    strokeFacesSqr.clear();tri_hits.clear(); strokePoints.clear();

    Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;
    for (auto pt: splinePoints) {
      Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core().proj,viewer.core().viewport);
      std::vector<igl::Hit> hits;

      igl::unproject_in_mesh(Eigen::Vector2f(proj[0],proj[1]), viewer.core().view,// viewer.core().model,
        viewer.core().proj, viewer.core().viewport, V, F, pt,hits);


      if (hits.size()> 0) {
        int fi = hits[0].id;

        d = proj[2];
        strokeFacesSqr.insert(int(fi/2)); // square face, so divide by two from triangle faces
        tri_hits.push_back(hits[0]);
        strokePoints.push_back(pt);
      }
    }
}

void Lasso::strokeFinishCurve(int spline_pt_number) {
  if (spline_pt_number >= 3) {
    set_spline_points_from_stroke_points(spline_pt_number);
    updateStrokeAndHitsFromSpline();
  }
  std::vector<int> hits_per_face;
  int cur_face = tri_hits[0].id/2; int cur_face_n = 0; // it's zero because we'll iteratre the same hit again..
  for (auto tri_hit: tri_hits) {
    int new_face = tri_hit.id/2;
    if (new_face == cur_face) {
      cur_face_n++;
    } else {
      hits_per_face.push_back(cur_face_n);
      cur_face_n = 1;
      cur_face = new_face;
    }
  }
  hits_per_face.push_back(cur_face_n);
  //for (auto f_hits: hits_per_face) {cout << "face hits: " << f_hits << endl;}

  Edge closest_edge; std::pair<double,double> edge_coords; double distance;

  // first pt should fit the first edge, then use last point of every face
  int sqr_f = tri_hits[0].id/2;
  std::vector<Edge> edges; get_edges_from_face(F,sqr_f,edges);
  get_hit_closest_edge_coordinates(tri_hits[0], edges, closest_edge, edge_coords);
  const_edges.push_back(closest_edge); edge_coordinates.push_back(edge_coords);

  int hits_idx = 0;
  // Go through all hits by faces
  for (int i = 0; i < hits_per_face.size(); i++) {
    int hits_num = hits_per_face[i];
    //cout << "hits_num = " << hits_num << endl;
    int sqr_f = tri_hits[hits_idx].id/2;
    edges.clear();get_edges_from_face(F,sqr_f, edges);

    // Get the last hit in the face
    hits_idx += hits_num-1;
    //cout << "adding point with hits_idx = " << hits_idx << endl;
    // remove previous edge
    edges.erase(std::remove(edges.begin(), edges.end(), const_edges[i]), edges.end());

    // find next closest edge for the last edge
    get_hit_closest_edge_coordinates(tri_hits[hits_idx], edges, closest_edge, edge_coords);
    const_edges.push_back(closest_edge); edge_coordinates.push_back(edge_coords);
    /*
    int cur_hits = hits_idx;
    for (; hits_idx < cur_hits + hits_num; hits_idx++) {
      igl::Hit hit = tri_hits[hits_idx];
      Eigen::RowVector3d bc;
      bc << 1.0-hit.u-hit.v, hit.u, hit.v;

      cout << "sqr_f = " << sqr_f << " hits_idx = " << hits_idx << " bc = " << bc << endl;
    }
    */
    hits_idx++;
  }
  //cout << "hits_per_face.size() = " << hits_per_face.size() << endl;
  //cout << "const_edges.size() = " << const_edges.size() << endl;
  get_two_face_neighbour_groups();

  strokePoints.clear();
}


void Lasso::strokeFinish(Eigen::VectorXi &selected_vertices)
{

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.core().model;

  //marker for selected vertices
  Eigen::VectorXi is_selected; is_selected.setZero(V.rows(),1);

  //project all vertices, check which ones land inside the polyline
  for (int vi =0; vi<V.rows(); ++vi)
  {
    Eigen::Vector3f vertex = V.row(vi).transpose().cast<float>();
    Eigen::Vector3f proj = igl::project(vertex, modelview, viewer.core().proj,viewer.core().viewport);
    if (igl::point_in_poly(stroke2DPoints, proj[0], proj[1]))
      is_selected[vi] = 1;

  }

  //the selection might consist of front facing and back facing facets.
  //we will only select the connected component that is the most frontal

  //first, isolate the faces that have at least one selected vertex
  int nf = 0;
  Eigen::MatrixXi Fsel(F.rows(),3);
  for (int fi = 0; fi<F.rows(); ++fi)
  {
    bool mark = false;
    for (int i = 0; i<3; ++i)
      if (is_selected[F(fi,i)])
      {
        mark = true;
        break;
      }
    if (mark)
      Fsel.row(nf++) = F.row(fi);
  }
  Fsel.conservativeResize(nf, Eigen::NoChange);
  //compute their barycenters
  Eigen::MatrixXd MFsel;
  igl::barycenter(V, Fsel, MFsel);

  //now, find all connected components of selected faces
  Eigen::VectorXi cid;
  igl::facet_components(Fsel, cid);

  //compute centroids of connected components
  int ncomp = cid.maxCoeff()+1;
  Eigen::MatrixXd region_centroids;
  region_centroids.setZero(ncomp,3);
  Eigen::VectorXi total; total.setZero(ncomp,1);
  for (long fi = 0; fi<Fsel.rows(); ++fi)
  {
    int r = cid[fi];
    region_centroids.row(r) += MFsel.row(fi);
    total[r]++;
  }
  for (long i = 0; i<ncomp; ++i)
    region_centroids.row(i) = region_centroids.row(i).array()/total[i];

  //project all centroids and isolate only the most frontal one
  float mind = 1e10;
  int r = -1;
  for (long i = 0; i<ncomp; ++i)
  {
    Eigen::Vector3f t = region_centroids.row(i).transpose().cast<float>();
    Eigen::Vector3f proj = igl::project(t,
                                        modelview,
                                        viewer.core().proj,
                                        viewer.core().viewport);
    float depth = proj[2];
    if (mind > depth)
    {
      r = i;
      mind = depth;
    }
  }

  //all vertices belonging to other components are unmarked
  for (long fi = 0; fi<Fsel.rows(); ++fi)
  {
    if (cid[fi] != r)
    for (int i = 0; i<3; ++i)
      is_selected[Fsel(fi,i)] = 0;
  }

  //return the selected vertices
  int nc = is_selected.sum();
  selected_vertices.resize(nc,1);
  int num = 0;
  for (int vi =0; vi<V.rows(); ++vi)
    if (is_selected[vi])
      selected_vertices[num++] = vi;



  strokePoints.clear();
  stroke2DPoints.clear();
  d = -1;
}
