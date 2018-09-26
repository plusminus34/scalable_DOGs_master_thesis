import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.prepared import prep

from graph_planar_embeddings import *
from descartes import PolygonPatch
from shapely.ops import cascaded_union,linemerge,split

from drawing import *


def crease_pattern(border_polygon, lines):
	remove_points_outside_border(border_polygon, lines)
	snap_points_towards_another(border_polygon, lines)
	"""
	G = build_planar_graph(border_polygon, lines)
	f = get_graph_faces(G)
	faces_polygons = get_face_polygons(G,f)
	"""

def split_line_to_polygon(border_polygon, line):
	res = split(line,border_polygon)
	print 'res = ', res
	res2 = cascaded_union(res)
	print 'res2 = ', res2
	new_line = linemerge(res2)
	print 'new_line = ', new_line
	line =  linemerge(cascaded_union(split(line, border_polygon)))
	print 'line before = ', line
	return line

# now got a splitted line

# also splits the lines if needed
def remove_points_outside_border(border_polygon, lines):
	prepared_polygon = prep(border_polygon)
	new_lines = []
	for line in lines:
		line = split_line_to_polygon(border_polygon, line)
		#print 'line = ', line
		#print 'here'
		#bla = prepared_polygon.intersection(line)
		#some_points = MultiPoint(line.coords[:])
		#prepared_polygon.intersects(some_points)
		#hits = filter(prepared_polygon.intersects, some_points)
		#line = 
		#print 'hits = ', hits
		new_lines.append(line)
	return new_lines

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

def test_plot_polygon_and_lines(fig_num, border_polygon, lines):
	fig = plt.figure(fig_num, figsize=(5,5), dpi=90)
	ax = fig.add_subplot(111)

	# plot border polygon
	pol_patch = PolygonPatch(border_polygon)
	ax.add_patch(pol_patch)

	# plot lines
	for line in lines:
		plot_line(ax, line)
		plot_coords(ax, line)

def test_crease_pattern():
	eps = 1e-2
	border_polygon = Polygon([(0, 0), (0,1), (1, 1), (1, 0)])
	line1 = LineString([(-eps,-eps), (0.1,0), (0.5,0.5), (2,1)])
	lines = [line1]

	# figure before remove_points_outside_border
	test_plot_polygon_and_lines(1,border_polygon,lines)	


	line1.intersects(border_polygon)
	border_polygon.intersects(line1)
	print 'worked!'
	# figure after
	lines = remove_points_outside_border(border_polygon, [line1])
	test_plot_polygon_and_lines(2,border_polygon,lines)	

	# new figure


	# show all
	plt.show()
	
if __name__ == "__main__":
	test_crease_pattern()