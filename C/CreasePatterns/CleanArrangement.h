#pragma once

#include <Eigen/Dense>
#include "PlanarArrangement.h"
#include "PatternBoundary.h"


// Filter segments outside of the boundary, and snap vertices that are close to the boundary to it
void snap_to_boundary(const PatternBoundary& patternBoundary, const CGAL::Bbox_2& bbox, PlanarArrangement& inArrangement, PlanarArrangement& outArrangement);