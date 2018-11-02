#include <igl/opengl/glfw/Viewer.h>

#include <igl/combine.h>

//#include "CreasePatterns/PatternBoundary.h"
/*
#include "CreasePatterns/DogCreasePattern.h"
#include "CreasePatterns/OrthogonalGrid.h"
#include "CreasePatterns/PlanarArrangement.h"
#include "CreasePatterns/SVGReader.h"


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
*/

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>
#include "CreasePatterns/arr_print.h"

#include <CGAL/config.h>
// Instantiate the traits class using a user-defined kernel
// and Segment_traits_2.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_2;
// Identical instantiation can be achieved using the default Kernel:
// typedef CGAL::Arr_polyline_traits_2<>                    Geom_traits_2;
typedef Geom_traits_2::Point_2                            Point_2;
typedef Geom_traits_2::Segment_2                          Segment_2;
typedef Geom_traits_2::Curve_2                            Polyline_2;
typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;

using namespace std;

int main(int argc, char* argv[]) {
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);
  Geom_traits_2::Construct_curve_2 polyline_construct =
    traits.construct_curve_2_object();
  Point_2 points1[5];
  points1[0] = Point_2(261.5,300);
  points1[1] = Point_2(247.5,113.5);
  points1[2] = Point_2(0,67.5);
  
  Polyline_2 pi1 = polyline_construct(&points1[0], &points1[3]);
  insert(arr, pi1);
  print_arrangement(arr);
  cout << "arr.number_of_vertices() = " << arr.number_of_vertices() << endl;

  std::cout << "My CGAL library is " <<  CGAL_VERSION_NR << " (1MMmmb1000)" << std::endl; 
  std::cout << std::endl;
  std::cout << "where MM is the major number release, mm is the minor number release, and "
            << "b is the bug fixing number release." << std::endl;

  return 0;

}