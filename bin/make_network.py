#!/usr/local/bin/python3.4

import numpy as np
import sys

ifilename = sys.argv[1]
nodefilename = sys.argv[2]
ofilename = sys.argv[3]

ifile = open(ifilename, 'r')
ofile = open(ofilename, 'w')
nodenum = np.loadtxt(nodefilename, dtype=int)

oarr=[]

for i in range (0, nodenum):
	oarr.append([i+1])


for line in ifile:
	pair=line.split()
	node1=int(pair[0])
	node2=int(pair[1])	
	oarr[node1-1].append(node2)
	oarr[node2-1].append(node1)

### Save graphml file fo
## 1. Header file.
ofile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
ofile.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns"\n')
ofile.write('xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n')
ofile.write('xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns\n')
ofile.write('http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\n">')
## 2. Graph information
ofile.write('<graph id="G" edgedefault="undirected">\n')
## 3. Node information
for i in range (0, nodenum):
	ofile.write('<node id="'+str(i+1)+'"/>\n')
## 4. Edge information
for i in range (0, nodenum):
	for j in range (1, len(oarr[i])):
		source=oarr[i][0]; target=oarr[i][j]
		if source < target:
			ofile.write('<edge source="'+str(source)+'" target="'+str(target)+'"/>\n')
## 5. End
ofile.write('</graph>\n')
ofile.write('</graphml>\n')

ifile.close()
ofile.close()
	
