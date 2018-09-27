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


def crease_pattern(border_polygon, polylines):
	# snap_points_towards_another
	snap_points_towards_another(border_polygon, polylines)
	# figure after
	border_polygon,polylines = remove_points_outside_border(border_polygon, polylines)
	face_polygons = build_polygons(border_polygon, polylines)
	return face_polygons

def add_curve_edges_to_graph(G,vertices,coords):
	for idx in range(coords.shape[0]-1):
		pos1, pos2 = coords[idx], coords[idx+1]
		idx1 = np.where((vertices == pos1).all(axis=1))[0][0]
		idx2 = np.where((vertices == pos2).all(axis=1))[0][0]
		print 'adding edge between ', idx1, ' and ', idx2
		G.add_edge(idx1,idx2)

def build_polygons(border_polygon, polylines):
	G = build_planar_graph(border_polygon, polylines)

	faces = get_graph_faces(G)
	polygons = []
	print "nx.get_node_attributes(G,'pos') = ", nx.get_node_attributes(G,'pos')
	positions = nx.get_node_attributes(G,'pos')
	#print 'faces = ', faces
	for f in faces:
		print 'face with ', f
		indices = [pt[0] for pt in f]
		vals = [positions[idx] for idx in indices]
		#print 'vals = ', vals
		new_poly = Polygon(vals)
		if new_poly.area != border_polygon.area:
			polygons.append(new_poly)
	return polygons

def build_planar_graph(border_polygon, polylines):
	G = nx.Graph()

	# build vertices
	vertices = np.array(border_polygon.exterior.coords[:])
	for pol in polylines:
		vertices = np.concatenate((vertices, pol.coords[:]))
	# get unique vertices
	vertices = unique_rows(vertices)
	v_n = vertices.shape[0]
	print 'v_n = ', v_n
	for v in range(v_n):
		G.add_node(v, pos = vertices[v,:])

	# add edges from lines
	#for pol in polylines:
	pol_v = np.array(border_polygon.exterior.coords[:])
	lines_v = [np.array(pol.coords[:]) for pol in polylines]
	add_curve_edges_to_graph(G,vertices,pol_v)
	for coords in lines_v:
		add_curve_edges_to_graph(G,vertices,coords)
	print 'G.edges() = ', G.edges()
	return G

def split_line_to_polygon(polygon, line):
	return linemerge(cascaded_union(split(line, polygon)))

def filter_line_points_outside_polygon(polygon, line):
	points_inside_polygon = filter(polygon.intersects, MultiPoint(line.coords[:]))
	return LineString(points_inside_polygon)

def split_polygon_by_line(polygon, line):
	return cascaded_union(split(polygon,line))

# also splits the lines if needed, and also split the polygon itself
def remove_points_outside_border(border_polygon, polylines):
	#prepared_polygon = prep(border_polygon)
	new_poly_lines = []
	for pol_line in polylines:
		pol_line = split_line_to_polygon(border_polygon, pol_line)
		pol_line = filter_line_points_outside_polygon(border_polygon, pol_line)
		new_poly_lines.append(pol_line)

		border_polygon = split_polygon_by_line(border_polygon, pol_line)
	return border_polygon, new_poly_lines

def snap_points_towards_another(border_polygon, polylines):
	pass

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

"""

1) Remove all points that are outside the border (call prepared geometry and then contained and/or intersects since its a polygon)
2) Snap points towards one another (n to n-1, than n and n-1 to n-2 and more). Consider defining that epsilon in a smart way (still much smaller than most edge lengths?)
3) Build planar graph. The vertices are all the points that are left. The edges will be defined by the curves.
4) Create the faces. From the faces create the polygons again, and plot all in different colors.
5) Test on various svg's.
6) Create a grid. Subdivide the grid to add internal curves. Then intersect various polygons with the grid.
7) Rasterize various parts and save the bookeeping information

"""

def test_plot_polygon_and_lines(fig_num, border_polygon, polylines):
	fig = plt.figure(fig_num, figsize=(5,5), dpi=90)
	ax = fig.add_subplot(111)

	# plot border polygon
	pol_patch = PolygonPatch(border_polygon)
	ax.add_patch(pol_patch)

	# plot lines
	for line in polylines:
		plot_line(ax, line)
		plot_coords(ax, line)

def test_crease_pattern(border_polygon = [], polylines = []):
	if border_polygon == []:
		eps = 1e-2
		#border_polygon = prep(Polygon([(0, 0), (0,1), (1, 1), (1, 0)]))
		border_polygon = Polygon([(0, 0), (0,1), (1, 1), (1, 0)])
		line1 = LineString([(-eps,-eps), (0.1,0), (0.5,0.5), (2,1)])
		line2 = LineString([(0,0.5), (1.2,1.5)])
		line3 = LineString([(0,0.4), (1 +eps,1+eps)])
		polylines = [line1, line2, line3]

	# figure before remove_points_outside_border
	test_plot_polygon_and_lines(1,border_polygon,polylines)	

	# snap_points_towards_another
	snap_points_towards_another(border_polygon, polylines)

	# figure after
	border_polygon,polylines = remove_points_outside_border(border_polygon, polylines)
	test_plot_polygon_and_lines(3,border_polygon,polylines)

	face_polygons = build_polygons(border_polygon, polylines)
	print 'face_polygons = ', face_polygons

	fig = plt.figure(4, figsize=(5,5), dpi=90)
	ax = fig.add_subplot(111)

	pol_colors = get_spaced_colors(len(face_polygons))
	i = 0
	for pol in face_polygons:
		color = np.array(pol_colors[i])/255.
		print 'color = ', color
		print 'pol.area = ', pol.area
		pol_patch = PolygonPatch(pol, facecolor=color)
		print 'pol_patch = ', pol_patch
		ax.add_patch(pol_patch)
		i += 1
	# plot lines
	for line in polylines:
		plot_line(ax, line)
		plot_coords(ax, line)

	# show all
	plt.show()
	
if __name__ == "__main__":
	test_crease_pattern()