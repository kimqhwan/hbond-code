import networkx as nx
import numpy as np
import sys
import argparse
import time

def parsecmd():
	parser=argparse.ArgumentParser(description='description')
    
	parser.add_argument('-f', type=argparse.FileType('r'), dest='ifile')
	parser.add_argument('-o', type=argparse.FileType('w'), dest='ofile')
	parser.add_argument('--name',nargs='*',type=str, dest='name')
	parser.add_argument('--number', nargs='*', type=int, dest='number')
	options = parser.parse_args()

	return options


start = time.time()

options=parsecmd()

g = nx.read_graphml(options.ifile, node_type=int)

total_mol_num = np.sum(options.number)
mol_num_arr = options.number

o_arr = []

## Check availability of graph and set subgraph
sub_node_arr = []
for i in range (0, len(options.name)):
	sub_node_arr.append([])

attribute_dic = nx.get_node_attributes(g, 'name')
for node in attribute_dic:
	attribute = attribute_dic[node]
	for i in range (0, len(options.name)): 
		if options.name[i] == attribute:
			sub_node_arr[i].append(node)

for i in range (0, len(options.number)): 
	if options.number[i] != len(sub_node_arr[i]):
		print 'Your input number != your graph nodes number'


## Calculate graph property
length = nx.all_pairs_shortest_path_length(g)
buf_spl = 0.
for i in range (1, total_mol_num+1):
	buf_spl_per_mol = 0.
	for node in length[i]:
		buf_spl_per_mol += length[i][node]
	buf_spl_per_mol /= (total_mol_num-1)
	buf_spl += buf_spl_per_mol
buf_spl /= total_mol_num

o_arr.append(buf_spl)


## Set subgraph and calculate its property
if len(options.number) >= 2:
	for i in range (0, len(options.number)):
		sub_g = g.subgraph(sub_node_arr[i])
		length = nx.all_pairs_shortest_path_length(sub_g)

		buf_spl = 0.
		for j in sub_node_arr[i]:
			buf_spl_per_mol = 0.
			for node in length[j]:
				buf_spl_per_mol += length[j][node]
			buf_spl_per_mol /= (options.number[i]-1)
			buf_spl += buf_spl_per_mol
		buf_spl /= options.number[i]

		o_arr.append(buf_spl)

## Save data
np.savetxt(options.ofile, o_arr, fmt='%5f')

end = time.time() - start
print 'Total time is '+str(end)


