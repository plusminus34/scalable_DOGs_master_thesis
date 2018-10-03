from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.ops import transform
from graph_planar_embeddings import *
import numpy as np
import networkx as nx
from geometry_utils import *

def grid_from_boundary(border_polygon, res_x = 20, res_y = 20):
	assert is_planar_polygon_rectangle(border_polygon), "Unsupported: Currently only supporting rectangular polygons"
	# For now assume a rectangular border polygon
	minx, miny, maxx, maxy = border_polygon.bounds
	step_x = (maxx-minx)/res_x
	step_y = (maxy-miny)/res_y
	grid_lines = []
	for y_i in range(res_y+1):
		points = [(minx + step_x*x_i,miny+y_i*step_y) for x_i in range(res_x+1)]
		grid_lines.append(LineString(points))

	for x_i in range(res_x+1):
		points = [(minx + step_x*x_i,miny+y_i*step_y) for y_i in range(res_y+1)]
		grid_lines.append(LineString(points))

	return grid_lines

def split_grid_by_intersections(grid, intersections):
	#assert is_planar_polygon_rectangle(border_polygon), "Unsupported: Currently only supporting rectangular polygons"
	# For now assume a rectangular border polygon
	#minx, miny, maxx, maxy = border_polygon.bounds

	for pt in intersections:
		y = pt[1]
		#print 'sub y = ', y
		grid = subdivide_planar_grid_at_x(grid,y)
	return grid

# faces to graph, then sort edges and "unite" them (0,1 and 1,2 to 0,2..)
def polygons_to_grid(polygons):
	grid = []
	return grid

def grid_squares_to_mesh_numpy(grid_squares):
	V = np.empty((0,2))
	for p in grid_squares:
		V = np.concatenate((V, p.exterior.coords[:]))
	V = unique_rows(V)
	F = np.empty((len(grid_squares),4))
	f_i = 0
	for sqr in grid_squares:
		#print 'sqr.exterior.coords[:] = ', sqr.exterior.coords[:]
		pos1, pos2, pos3, pos4= sqr.exterior.coords[0:4]
		idx1, idx2 = vertices_pos_to_index(V,pos1),vertices_pos_to_index(V,pos2)
		idx3, idx4 = vertices_pos_to_index(V,pos3),vertices_pos_to_index(V,pos4)
		F[f_i,:] = idx1,idx2,idx3,idx4
		f_i += 1
	return V,F

# build a graph and find faces
def grid_to_polygons(grid):
	G = nx.Graph()
	vertices = np.empty((0,2))
	vertices = add_curve_vertices_to_graph(G, grid, vertices)

	# add edges
	lines_v = [np.array(pol.coords[:]) for pol in grid]

	for coords in lines_v:
		add_curve_edges_to_graph(G,vertices,coords)
	return graph_to_polygons(G)

def subdivide_planar_grid_at_x(grid,sub_y):
	#print 'grid len before = ', len(grid)
	idx = 0
	for line in grid:
		x,y = np.array(line.xy[0]),np.array(line.xy[1])
		if np.all(y==y[0]):
			# This is a curve with constant 'y'
			#print 'cur_y = ', y[0]
			if y[0] == sub_y:
				return grid # do nothing
			if y[0] > sub_y:
				new_y = line
				new_y = transform(lambda x,y: [x,sub_y], new_y)
				for i in range(len(grid)):
					grid[i] = split_line_to_geometry(extend_polyline_by_epsilon(new_y,1e-4),grid[i])
				grid = grid[:idx]+[new_y]+grid[idx:]
				return grid
		idx +=1

	return grid

def is_planar_polygon_rectangle(poly):
	minx, miny, maxx, maxy = poly.bounds
	rectangle = Polygon([(minx,miny),(maxx, miny),(maxx,maxy),(minx,maxy)])
	return poly == rectangle