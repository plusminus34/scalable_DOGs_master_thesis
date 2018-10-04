import numpy as np
import os, sys
sys.path.insert(0, os.getcwd() + "/../../libigl/python")
sys.path.insert(0, os.getcwd() + "/../../libigl/external/nanogui/build/python")
import math
import pyigl as igl
from scipy.spatial import Delaunay

from polygons_to_orthogonal_grids import *

class Edge:
	def __init__(self, v1, v2):
		self.v1 = v1
		self.v2 = v2

class FoldConst:
	def __init__(self, edge1, w1, edge2, w2):
		self.edge1, self.w1, self.edge2, self.w2 = edge1,w1,edge2,w2

def dog_builder(V_list,F_list, polylines):
	# return global V,F, then V and F for rendering, and constraints 
	# The global F just concatenates the vertices and faces list (after adding indices to faces)
	V,F,V_list_pos_dict = get_global_VF(V_list, F_list)

	#
	# The constraints data is defined per vertex on the polylines.
	# A polyline vertex is defined uniquely by its 2D coordinates (so there won't be doubles). 
	# No need to save it but this is a good way to filter multiples of those.
	# Each one of those is either exactly a vertex, or a point on an edge.
	# In any case it can be represented as a point of an edge. So for every mesh that does contain this point:
	#	We save two indices of vertices and an appropriate weight
	#   If we have 2 pairs of this we need one constraint, if we have 3 pairs 2, and in general for n we need n-1 (1 to 2, 2 to 3, etc)
	#   So at the end we just need to save a list of constraints of the form:
	#				2 points on edges being equal: (v_i1,v_i2) 0<=w_1<=1, (v_j1,v_j2) 0<=w_2<=1

	# So the function should find all unique points lying on polylines, then find all meshes containing those,
	# 	then find these as points on (some) edges for each one (n such edges), and save those
	#										 (the optimization code could use them as n-1 constraints)
	#
	polylines_v = get_polylines_v(polylines)
	fold_consts, unique_vertex_on_edges = get_fold_consts(V_list, F_list, V_list_pos_dict, polylines_v)
	print 'len(fold_consts) = ', len(fold_consts)
	assert len(unique_vertex_on_edges) == polylines_v.shape[0], 'Lengths should match'
	
	#
	# Rendering requires culling the faces. 
	# The triangulation can be an arbitrary one (triangulate) if you take the unique vertices locations together with the polylines.
	# If we have a unique order for the polylines vertices (given by another list of edges such that no vertex is there twice)
	#   then we just need the list of edges.
	# So basically: Save a list of vertices indices in the global mesh (for the inner vertices positions, meaning avoiding duplications)
	# 				A list of edge and weights (vertices on edges) for polylines
	#				And an F for triangulation that assumes this order of vertices (so for initialization build it and triangulate)
	#
	# In case the triangulation looks bad, do a manual one

	# TODO1: First find unique indices of V and save them, then delaunay that
	V_render = np.concatenate((np.copy(V),polylines_v))
	F = Delaunay(V_render).simplices
	# use the constraints to get a list of edge and weights vertices
	print 'type(F) = ', type(F)
	print 'F = ', F.shape
	# TODO2: Get the parameters of the edges (edge vertices for each constraint)
	# This could be returned also from get_fold_consts..

	plt.triplot(V_render[:,0], V_render[:,1], F.copy())
	# TODO3: Return the data and save it in a file
	# TODO4: Change your code to read that data, constraints, and rendering

def get_polylines_v(polylines):
	# Get polylines vertices positions (filter for uniqueness)
	vertices = np.empty((0,2))
	for p in polylines:
		vertices = np.concatenate((vertices, np.array(p.coords[:])))
	vertices = unique_rows(vertices)
	return vertices

def get_fold_consts(V_list, F_list, V_list_pos_dict, polylines_v):
	fold_consts = []
	unique_vertex_on_edges = []

	# Find their intersections with the meshes (build polylines from them)
	# Locate the intersections vertices edges and weights
	for v in polylines_v:
		v_edges = []
		edge_weights = []
		first_submesh = True
		for mesh_i in range(len(V_list)):
			V,F = V_list[mesh_i], F_list[mesh_i]
			lines = get_mesh_non_unique_edges_positions(V,F)
			# Check if the point is in the polygon (this could be faster I guess)
			v1_pos,v2_pos,t = find_v_on_lines(v,lines)
			#print 'v = ', v
			if v1_pos:
				#print 'v1_pos = ', v1_pos, ' v2_pos = ', v2_pos
				v1,v2 = V_list_pos_dict[mesh_i][v1_pos], V_list_pos_dict[mesh_i][v2_pos]
				#print 'v1 = ',v1, ' v2 = ', v2
				v_edges.append(Edge(v1,v2))
				edge_weights.append(t)
				if first_submesh:
					unique_vertex_on_edges.append((v1,v2,t))
					first_submesh = False
		# now go through
		const_num = len(v_edges) - 1
		print 'constraint numbers = ', const_num, ' with v = ', v
		#if const_num > 1:
			#print 'constraint numbers = ', const_num, ' with v = ', v
		assert const_num >= 1, 'There should be at least one constraint per fold vertex'
		for c_i in range(const_num):
			fold_consts.append(FoldConst(v_edges[c_i],edge_weights[c_i],v_edges[c_i+1],edge_weights[c_i+1]))

	return fold_consts, unique_vertex_on_edges

def get_mesh_non_unique_edges_positions(V,F):
	lines = []
	for f in F:
		#print 'f = ', f		
		lines.append(LineString([V[f[0]],V[f[1]]]))
		lines.append(LineString([V[f[1]],V[f[2]]]))
		lines.append(LineString([V[f[2]],V[f[3]]]))
		lines.append(LineString([V[f[3]],V[f[0]]]))
	return lines

def find_v_on_lines(v, lines):
	for l in lines:
		if Point(v).intersects(l):
			# solve t*v1 + (1-t)*v2 = v so t*v1-t*v2=v-v2, so t*(v1-v2) = v-v2, so t = (v-v2)/(v1-v2)
			v1,v2 = np.array(l.coords[0]), np.array(l.coords[1])
			#print '(v[0]-v2[0]) = ', (v[0]-v2[0])
			#print '(v1[0]-v2[0]) = ', (v1[0]-v2[0])
			if (v1[0]!=v2[0]):
				t = (v[0]-v2[0])/(v1[0]-v2[0])
			else:
				t = (v[1]-v2[1])/(v1[1]-v2[1])
			#print 't = ', t
			return tuple(v1),tuple(v2),t # return vertices as tuples as we need them as hashable dictionary keys
	return None,None,None

def get_global_VF(V_list, F_list):
	V_list_pos_dict = []
	V,F = V_list[0], F_list[0]
	max_f = F.max()
	for mesh_i in range(1, len(V_list)):
		V = np.concatenate((V, V_list[mesh_i]))
		F = np.concatenate((F, F_list[mesh_i] + F.max() + 1))

	# For every submesh creates a dictionary between vertex positions to the index of the
	# 	 vertex in the global mesh
	base_v_offset = 0
	for submesh_i in range(len(V_list)):
		V_sub_dict = {}
		for row_i in range(V_list[submesh_i].shape[0]):
			pos = V_list[submesh_i][row_i]
			V_sub_dict[tuple(pos)] = row_i + base_v_offset
		base_v_offset += V_list[submesh_i].shape[0]
		V_list_pos_dict.append(V_sub_dict)
	#print 'V_list_pos_dict = ', V_list_pos_dict
	return V,F, V_list_pos_dict

def test_dog_builder(svg_file):
	f, ((ax1, ax2),(ax3,ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
	border_poly,polylines = svg_creases_to_polygonal_data(svg_file)
	face_polygons, polylines = crease_pattern(border_poly, polylines)

	res_x, res_y = 5,5
	V_list, F_list, grid, grid_polylines = polygons_to_orthogonal_grids(face_polygons, border_poly, polylines, res_x, res_y)
	face_polygons_num = len(V_list)
	plot_face_polygons(face_polygons, grid_polylines, ax1, 'Faces and grid (' + str(face_polygons_num) + ' faces)')
	#pol_patch = PolygonPatch(border_poly)
	#ax3.add_patch(pol_patch)
	plot_grid(grid, ax1, 1.5, '#ffffff')
	for line in grid_polylines:
		plot_coords(ax1, line, '#444444')
		plot_line(ax1, line, 0.5, '#dddddd') # line width = 1
	plot_line(ax1,LineString(border_poly.exterior.coords),1,'#000000')
	for f in F_list:
		print 'f.rows() = ', f.shape[0]

	dog_builder(V_list, F_list, grid_polylines)

	#lines = get_mesh_non_unique_edges_positions(V_list[2],F_list[2])
	#plot_border_polygon_and_lines(ax2,border_poly, lines)

	#lines = get_mesh_non_unique_edges_positions(V_list[3],F_list[3])
	#plot_border_polygon_and_lines(ax3,border_poly, lines)
	plt.show()

if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_builder(sys.argv[1])
	else:
		test_dog_builder("../crease_patterns/1_curve.svg")
