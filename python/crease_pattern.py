import networkx as nx
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, wsvg
from svg_utils import *
from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.prepared import prep

from graph_planar_embeddings import *
from shapely.ops import cascaded_union,linemerge,split,snap

from drawing import *
from planar_dog import *

def crease_pattern(border_polygon, polylines):
	border_polygon, polylines = snap_polygons_border_to_another(border_polygon, polylines)
	face_polygons = build_polygons(border_polygon, polylines)
	return face_polygons, polylines

def snap_polygons_border_to_another(border_polygon, polylines):
	# snap_points_towards_another
	snap_curves_towards_another(polylines)
	
	# remove points outside border and snap polygons to another
	border_polygon,polylines = remove_points_outside_border(border_polygon, polylines)
	polylines = split_polylines_to_each_other(polylines)
	# snap_points_towards_another
	snap_curves_towards_another(polylines)
	return border_polygon, polylines

def add_curve_edges_to_graph(G,vertices,coords):
	for idx in range(coords.shape[0]-1):
		pos1, pos2 = coords[idx], coords[idx+1]
		idx1 = np.where((vertices == pos1).all(axis=1))[0][0]
		idx2 = np.where((vertices == pos2).all(axis=1))[0][0]
		#print 'adding edge between ', idx1, ' and ', idx2
		G.add_edge(idx1,idx2)

def build_polygons(border_polygon, polylines):
	G = build_planar_graph(border_polygon, polylines)

	faces = get_graph_faces(G)
	polygons = []
	#print "nx.get_node_attributes(G,'pos') = ", nx.get_node_attributes(G,'pos')
	positions = nx.get_node_attributes(G,'pos')
	#print 'faces = ', faces
	for f in faces:
		#print 'face with ', f
		indices = [pt[0] for pt in f]
		vals = [positions[idx] for idx in indices]
		#print 'vals = ', vals
		new_poly = Polygon(vals)
		if abs(new_poly.area - border_polygon.area) > 1e-4:
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
	#print 'v_n = ', v_n
	for v in range(v_n):
		G.add_node(v, pos = vertices[v,:])

	# add edges from lines
	#for pol in polylines:
	pol_v = np.array(border_polygon.exterior.coords[:])
	lines_v = [np.array(pol.coords[:]) for pol in polylines]
	add_curve_edges_to_graph(G,vertices,pol_v)
	for coords in lines_v:
		add_curve_edges_to_graph(G,vertices,coords)
	#print 'G.edges() = ', G.edges()
	return G

def split_line_to_geometry(polygon, line):
	#union = cascaded_union(split(line, polygon))
	split_res = split(line, polygon)
	# If there was a split, unite the lines
	if len(list(split_res)) > 1:
		return linemerge(cascaded_union(split_res))
	else: # otherwise return the original line
		return line

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
		pol_line = split_line_to_geometry(border_polygon, pol_line)
		pol_line = filter_line_points_outside_polygon(border_polygon, pol_line)
		new_poly_lines.append(pol_line)

		border_polygon = split_polygon_by_line(border_polygon, pol_line)
	return border_polygon, new_poly_lines

def extend_polyline_by_epsilon(polyline, eps):
	b_point, bb_point = np.array(polyline.coords[-1]), np.array(polyline.coords[-2])
	diff_vec = b_point-bb_point
	#print 'diff_vec = ', diff_vec
	direction_vec = diff_vec/np.linalg.norm(diff_vec)
	#print 'b_point = ', b_point
	#print 'diff_vec + eps*direction_vec = ', diff_vec + eps*direction_vec
	new_last_coord = bb_point + diff_vec + eps*direction_vec # extend it by an epsilon
	#print 'new_last_coord = ', new_last_coord

	new_coords = polyline.coords[0:-1]
	#print 'old_coords_minus_last = ', new_coords
	#print '(new_last_coord[0],new_last_coord[1]) = ', (new_last_coord[0],new_last_coord[1])
	new_coords.append((new_last_coord[0],new_last_coord[1]))
	#print 'new coords = ', new_coords
	return LineString(new_coords)

def split_polylines_to_each_other(polylines):
	polylines_new = []
	for pol_line in polylines:
		new_line = pol_line
		#print 'pol_line before = ', new_line
		for pol_line2 in polylines:
			# hack: due to precision errors in line splitting intersections, we need to create a "snapper" that is longer by some epsilon
			snapper = extend_polyline_by_epsilon(pol_line2, 1e-4)
			
			if pol_line != pol_line2 and pol_line.intersects(snapper):
				new_line = split_line_to_geometry(snapper,new_line)
				
		polylines_new.append(new_line)
	return polylines_new

def snap_curves_towards_another(polylines):
	# snap the first to second, then both to third, etc
	for idx1 in range(0, len(polylines)-1):
		for idx2 in range(0, idx1+2):
			if idx1!=idx2:
				polylines[idx2] = snap(polylines[idx2], polylines[idx1], 1e-3)
	return polylines

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def get_default_test_params():
	eps = 1e-2
	border_polygon = Polygon([(0, 0), (0,1), (1, 1), (1, 0)])
	line1 = LineString([(-eps,-eps), (0.1,0), (0.5,0.5), (2,1)])
	line2 = LineString([(0,0.5), (1.2,1.5)])
	line3 = LineString([(0,0.4), (1 +eps,1+eps)])
	polylines = [line1, line2, line3]
	return border_polygon, polylines

def test_crease_pattern(border_polygon = [], polylines = []):
	if border_polygon == []:
		border_polygon, polylines = get_default_test_params()

	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

	# figure before remove_points_outside_border
	plot_border_polygon_and_lines(ax1,border_polygon,polylines, 'SVG Input')	

	border_polygon, polylines = snap_polygons_border_to_another(border_polygon, polylines)

	plot_border_polygon_and_lines(ax2,border_polygon,polylines, 'Snapped and splitted curves')
	
	face_polygons = build_polygons(border_polygon, polylines)
	face_polygons_num = len(face_polygons)
	plot_face_polygons(face_polygons, polylines, ax3, 'Faces decomposition (' + str(face_polygons_num) + ' faces)')
	
	res_x,res_y = 25,25
	grid = grid_from_boundary(border_polygon, res_x,res_y)
	plot_face_polygons(face_polygons, polylines, ax4, 'Faces with grid')
	plot_grid(grid, ax4)

	# show all
	plt.show()


	
if __name__ == "__main__":
	test_crease_pattern()