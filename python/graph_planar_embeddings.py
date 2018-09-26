import networkx as nx
import matplotlib.pyplot as plt
import math

def get_graph_faces(G):
	comb_emb = comb_embedding_from_graph(G)
	return get_faces(G.edges(), comb_emb)

def comb_embedding_from_graph(G):
	# vertices should be zero based
	A=[(n, nbrdict) for n, nbrdict in G.adjacency()]
	pos=nx.get_node_attributes(G,'pos')

	comb_emb = {}
	for node in G.nodes():
		origin = pos[node]
		v_A = list(A[node][1])
		nb_pos = [[pos[nb][0],pos[nb][1]] for nb in v_A]
		#print 'pos[nb] = ', list[pos[nb]]
		print '--- node = ', node, ' at position = ', origin
		print 'nb_pos before = ', nb_pos
		print 'v_A before = ', v_A
		#nb_pos = sorted(nb_pos, key=lambda nb_pos: clockwiseangle_and_distance(nb_pos,origin))
		v_A = sorted(v_A, key=lambda node: clockwiseangle_and_distance(pos[node],origin))
		#print 'nb_pos after = ', nb_pos
		#sorted(v_A, key=lambda v_A: clockwiseangle_and_distance(v_A,nb_pos,origin))
		#v_A = [x for _,x in sorted(zip(nb_pos,v_A),key=lambda nb_pos: clockwiseangle_and_distance(nb_pos,origin))]
		print 'v_A after = ', v_A
		print '------\n'
		comb_emb[node] = v_A
	print 'comb_emb = ', comb_emb
	return comb_emb

def get_faces(edges,embedding):
	"""
	edges: is an undirected graph as a set of undirected edges
	embedding: is a combinatorial embedding dictionary. Format: v1:[v2,v3], v2:[v1], v3:[v1] clockwise ordering of neighbors at each vertex.)

	"""

	# Establish set of possible edges
	edgeset = set()
	for edge in edges: # edges is an undirected graph as a set of undirected edges
		edge = list(edge)
		edgeset |= set([(edge[0],edge[1]),(edge[1],edge[0])])

	# Storage for face paths
	faces = []
	path  = []
	for edge in edgeset:
		path.append(edge)
		edgeset -= set([edge])
		break  # (Only one iteration)

	# Trace faces
	while (len(edgeset) > 0):
		neighbors = embedding[path[-1][-1]]
		next_node = neighbors[(neighbors.index(path[-1][-2])+1)%(len(neighbors))]
		tup = (path[-1][-1],next_node)
		if tup == path[0]:
			faces.append(path)
			path = []
			for edge in edgeset:
				path.append(edge)
				edgeset -= set([edge])
				break  # (Only one iteration)
		else:
			path.append(tup)
			edgeset -= set([tup])
	if (len(path) != 0): faces.append(path)
	return iter(faces)

#def clockwiseangle_and_distance(vA, point, origin):
def clockwiseangle_and_distance(point, origin):
	refvec = [0, 1]
	# Vector between point and the origin: v = p - o
	vector = [point[0]-origin[0], point[1]-origin[1]]
	# Length of vector: ||v||
	lenvector = math.hypot(vector[0], vector[1])
	# If length is zero there is no angle
	if lenvector == 0:
		return -math.pi, 0
	# Normalize vector: v/||v||
	normalized = [vector[0]/lenvector, vector[1]/lenvector]
	dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
	diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
	angle = math.atan2(diffprod, dotprod)
	# Negative angles represent counter-clockwise angles so we need to subtract them 
	# from 2*pi (360 degrees)
	if angle < 0:
		return 2*math.pi+angle, lenvector
	# I return first the angle because that's the primary sorting criterium
	# but if two vectors have the same angle then the shorter distance should come first.
	#return angle, lenvector
	return angle

def comb_embedding_from_graph():
	G = nx.Graph()
	G.add_node(0,pos=(0,0))
	G.add_node(1,pos=(1,0))
	G.add_node(2,pos=(0,1))
	G.add_node(3,pos=(1,1))
	G.add_node(4,pos=(2,1))


	G.add_edge(0,1)
	G.add_edge(0,2)
	G.add_edge(1,3)
	G.add_edge(2,3)
	G.add_edge(0,3)
	G.add_edge(4,3)
	G.add_edge(4,1)

	pos=nx.get_node_attributes(G,'pos')
	nx.draw(G,pos)

	comb_emb = comb_embedding_from_graph(G)
	print 'G.edges() = ', G.edges()
	faces = get_faces(G.edges(), comb_emb)
	for f in faces:
		print 'found face = ', f
	#print 'faces = ', faces
	#pts = [[1,4],[2,4],[3,4],[1,3],[2,3],[3,3],[1,2],[2,2],[3,2]]
	#print 'pts before = ', pts
	#origin = [2,3]
	#sorted(pts, key=lambda point: clockwiseangle_and_distance(point,origin))
	#pts = sorted(pts, key=lambda point: clockwiseangle_and_distance(point,origin))
	#print 'pts after = ', pts
	plt.show()


if __name__ == "__main__":
	comb_embedding_from_graph()