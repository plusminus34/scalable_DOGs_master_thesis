import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *

from graph_planar_embeddings import *


def crease_pattern(border_polygon, lines):
	remove_points_outside_border(border_polygon, lines)
	snap_points_towards_another(border_polygon, lines)
	"""
	G = build_planar_graph(border_polygon, lines)
	f = get_graph_faces(G)
	faces_polygons = get_face_polygons(G,f)
	"""

def remove_points_outside_border(border_polygon, lines):
	pass

def snap_points_towards_another(border_polygon, lines):
	pass

"""

1) Remove all points that are outside the border (call prepared geometry and then contained and/or intersects since its a polygon)
2) Snap points towards one another (n to n-1, than n and n-1 to n-2 and more). Consider defining that epsilon in a smart way (still much smaller than most edge lengths?)
3) Build planar graph. The vertices are all the points that are left. The edges will be defined by the curves.
4) Create the faces. From the faces create the polygons again, and plot all in different colors.
5) Test on various svg's.
6) Create a grid. Subdivide the grid to add internal curves. Then intersect various polygons with the grid.
7) Rasterize various parts and save the bookeeping information

"""
	