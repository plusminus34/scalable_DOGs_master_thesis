#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_simple_point_location.h>

#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT                                        Number_type;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_2;
// Identical instantiation can be achieved using the default Kernel:
// typedef CGAL::Arr_polyline_traits_2<>                    Geom_traits_2;
typedef Geom_traits_2::Point_2                            Point_2;
typedef Geom_traits_2::Segment_2                          Segment_2;
typedef Geom_traits_2::Curve_2                            Polyline_2;
typedef Geom_traits_2::X_monotone_curve_2         		  Polyline_2_Monotone;

//typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;
typedef CGAL::Arrangement_with_history_2<Geom_traits_2> 	Arrangement_2; // Use arrangement with history to keep track of old polylines
typedef CGAL::Arr_simple_point_location<Arrangement_2>  Point_location;
typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;

typedef CGAL::Snap_rounding_traits_2<Kernel>     Snap_Traits;
typedef std::list<Polyline_2>                    Polyline_list_2;

typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2_without_history;