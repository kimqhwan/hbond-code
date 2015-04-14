import networkx as nx
import numpy as np

n=1000
k=20
p=[]

o_arr = []

for i in range (0, 20+1):
	p.append((2)**(i-20))
for i in range (0, len(p)):
	g_sw = nx.connected_watts_strogatz_graph(n,k,p[i],tries=500)
	average_short_path_length = nx.average_shortest_path_length(g_sw)
	clustering = nx.average_clustering(g_sw)
	o_arr.append([p[i], average_short_path_length, clustering])

short_path_length_0 = o_arr[0][1]
clustering_0 = o_arr[0][2]
for i in range (0, len(p)):
	o_arr[i][1] /= short_path_length_0
	o_arr[i][2] /= clustering_0


np.savetxt('sw.xvg', o_arr, fmt='%5f')

