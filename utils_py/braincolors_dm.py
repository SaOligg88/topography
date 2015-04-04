import sys
import xlrd
import re
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools
sys.path.insert(0,'./python-colormath/colormath')
from colormath.color_objects import LCHuvColor, LabColor, sRGBColor
from colormath.color_conversions import convert_color
from colormath.color_diff import delta_e_cie2000
from elementtree import ElementTree as et
import h5py 
from random import shuffle

##################
# SET PARAMETERS #
##################

# Output color text file
save_colors = 1
verbose = 1

# Choose plotting procedure(s) and save plots 
save_plots = 1
plot_colormap = 0  # Plot colormap
plot_graph = 1  # Plot whole graph
plot_subgraphs = 0  # Plot each individual colored subgraph

# Paths
in_dir = '../data/' #'../../results/'  # 'input/input_BrainCOLOR'
out_dir = '../results/colormaps/'
out_images = out_dir
vers = '3'

# Region adjacency table
in_table = in_dir + 'clus.mat' #'region_adjacency_matrix.xls'
#col_ID = 0  # column with region abbreviations (not used by the program)
#col_group = 0  # column with region group numbers
#col_abbr = 0  # column with region abbreviations
#col_name = 0  # column with full region names (not used by the program)
col_start_data = 0  # first column with data
row_start_data = 0  # first row with data
everyother = 1  # use <everyother> alternate row(s);
                # set to 2 for redundant labels across, e.g. brain hemispheres

# Color parameters
init_angle = 22.5
chroma = 70  # color "saturation" level
Lumas_init = np.array([30, 60, 15, 45])#, 50, 70])  # vary luminance values for adjacent colors
repeat_hues = 1  # repeat each hue for the different Lumas

# Graph layout
graph_node_size = 1000
graph_edge_width = 2
graph_font_size = 10
subgraph_node_size = 3000
subgraph_edge_width = 5
subgraph_font_size = 18
axis_buffer = 10

# Use weights in input adjacency matrix
use_input_weights = 1
# Plot whole graph with colored subgraphs
plot_graph_color = 1  
# Necessary for generating permutations
run_permutations = 1
n_perm = 100000

#########
# BEGIN #
#########

# Convert weighted connection matrix to weighted graph
A = np.array(h5py.File(in_table,'r').get('clus')['netScore'])
A = A/np.max(A)  # normalize weights
G = nx.from_numpy_matrix(A)
Ntotal = G.number_of_nodes()

roi_numbers = np.array(xrange(0,Ntotal)).astype(int)
roi_abbrs = [str(s).strip() for s in list(roi_numbers)]

for inode in range(Ntotal):
    G.node[inode]['abbr'] = roi_abbrs[inode] 
    G.node[inode]['code'] = roi_numbers[inode]
# Secondary parameters
color_angle = 360.0 / Ntotal
if repeat_hues:
    nLumas = len(Lumas_init)
    nangles = np.ceil(Ntotal / np.float(nLumas))
    color_angle = nLumas * color_angle
else: 
    nangles = Ntotal

# Plot graph
if plot_graph:
    plt.figure(Ntotal+2)
    labels={}
    for i in range(Ntotal):
        labels[i] = G.node[i]['abbr']
    pos = nx.graphviz_layout(G,prog="neato")
    nx.draw(G,pos,node_color='cyan',node_size=graph_node_size,width=graph_edge_width,with_labels=False)
    nx.draw_networkx_labels(G, pos, labels, font_size=graph_font_size, font_color='black')
    plt.axis('off')

out_colors = out_dir + 'region_colors_' + str(Ntotal) + '_' + vers + '.csv'
f = open(out_colors,'w')
N = Ntotal #len(glist)        
g = G #G.subgraph(glist)

# Define colormap as uniformly distributed colors in CIELch color space
Lumas = Lumas_init.copy()
while len(Lumas) < N: 
    Lumas = np.hstack((Lumas,Lumas_init))            
if repeat_hues:
    nangles_g = np.ceil(N / np.float(nLumas))
else: 
    nangles_g = N
hues = np.arange(init_angle, init_angle + nangles_g*color_angle, color_angle)
if repeat_hues:
    hues = np.hstack((hues * np.ones((nLumas,1))).transpose())
init_angle += nangles_g*color_angle

# Compute the differences between every pair of colors in the colormap
permutation_max = np.zeros(N)
NxN_matrix = np.zeros((N,N))

print "beginning perms"
if run_permutations:
    # Convert subgraph into an adjacency matrix (1 for adjacent pair of regions)
    neighbor_matrix = np.array(nx.to_numpy_matrix(g))
    if use_input_weights:
        pass
    else:
        neighbor_matrix = (neighbor_matrix > 0).astype(np.uint8)
    
    # Compute permutations of colors and color pair differences
    DEmax = 0
    permutations = []
    for s in xrange(0,n_perm): 
	permutations.append(np.array(np.random.permutation(range(0,N)))) #itertools.permutations(range(0,N),N)]
    print(" ".join([str(N),'labels,',str(len(permutations)),'permutations:']) )
    count = 0

    for permutation in permutations:
        delta_matrix = NxN_matrix.copy()
	
        for i1 in range(N):
            for i2 in range(N):
                if (i2 > i1) and (neighbor_matrix[i1,i2] > 0):
                    delta_matrix[i1,i2] = delta_e_cie2000(LabColor(Lumas[permutation[i1]],chroma,hues[permutation[i1]]),
							  LabColor(Lumas[permutation[i2]],chroma,hues[permutation[i2]]) )
        if use_input_weights:
             DE = np.sum((delta_matrix * neighbor_matrix))
        else:
             DE = np.sum(delta_matrix)

        # Store the color permutation with the maximum adjacency cost
        if DE > DEmax:
            DEmax = DE
            permutation_max = permutation
	
	if verbose:
	    print count
	count += 1

# Color subgraphs
if plot_graph_color:
    plt.figure(Ntotal+2)    
    for iN in range(N):
        ic = np.int(permutation_max[iN])
        lch = LCHuvColor(Lumas[ic],chroma,hues[ic])          
        rgb = convert_color(lch,sRGBColor)                                 
	coloring = np.array([rgb.rgb_r, rgb.rgb_g, rgb.rgb_b])
	coloring[coloring > 1] = 1
        nx.draw_networkx_nodes(g,pos,node_size=graph_node_size,nodelist=[g.node.keys()[iN]],node_color=coloring)
        f.write(str(coloring[0]) + ', ' + str(coloring[1]) + ', ' + str(coloring[2]) + '\n')  
                    
plt.savefig(out_images + "graph_" + vers + ".png")
f.close()
    
