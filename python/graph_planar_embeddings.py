def comb_embedding_from_graph(G):
    # vertices should be zero based
    A=[(n, nbrdict) for n, nbrdict in G.adjacency()]
    pos=nx.get_node_attributes(G,'pos')

    comb_emb = []
    for node in G.nodes()
        origin = pos[node]
        v_A = list[A[node]]
        print 'vA before = ', vA
        nb_pos = [pos[nb] for nb in v_A]
        sorted(v_A, key=lambda sort_clockwise_center: clockwiseangle_and_distance(nb_pos, origin))
        print 'vA after = ', vA

def Faces(edges,embedding)
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
        neighbors = self.embedding[path[-1][-1]]
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
    return angle, lenvector

"""
In [1]: import networkx as nx

In [2]: G=nx.Graph()

In [3]: G.add_node(1,pos=(1,1))

In [4]: G.add_node(2,pos=(2,2))

In [5]: G.add_edge(1,2)

In [6]: pos=nx.get_node_attributes(G,'pos')

In [7]: pos
Out[7]: {1: (1, 1), 2: (2, 2)}

In [8]: nx.draw(G,pos)
"""