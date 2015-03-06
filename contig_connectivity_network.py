#!/usr/bin/env python
###############################################################################
#
# __contig_connectivity_network__.py - description!
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2015"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################

# system imports
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
import networkx as nx
import matplotlib.pyplot as plt
#import os
#import errno
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NetworkPlot(object):
    def __init__(self, coverage_links_file):
        self.CLD = TFP.CoverageLinksData(coverage_links_file)
        self.contig_connectivity    = {}
        self.TG_connectivity        = {}
        self.node_colour            = {}
        self.node_colour_values     = {}
    
    def Wrapper(self, outfile, outfmt, label_font_size, draw_network_labels):
        # get contig connectivity data
        self.getData()
        
        # plot network
        self.Plot(outfile, outfmt, label_font_size, draw_network_labels)
        
    def getData(self):
        # Loop through contigs
        for i in range(len(self.CLD.TGs)):
            TG1 = self.CLD.TGs[i]
            
            for j in range(i+1, len(self.CLD.TGs)):
                TG2 = self.CLD.TGs[j]
                
                # Loop through TGs
                for contig1 in self.CLD.TGs_contigs[TG1]:
                    
                    for contig2 in self.CLD.TGs_contigs[TG2]:
                        
                        if contig1 == contig2:
                            self.addLink(TG1, TG2)
                            
    def addLink(self, TG1, TG2):
        try:
            self.TG_connectivity[TG1][TG2] = 1
        except KeyError:
            self.TG_connectivity[TG1] = {TG2:1}
        try:
            self.TG_connectivity[TG2][TG1] = 1
        except KeyError:
            self.TG_connectivity[TG2] = {TG1:1}
    
    def buildNetworkData(self, G):
        # get nodes
        nodes = self.CLD.TGs
        
        # add nodes to graph 
        G.add_nodes_from(nodes)
        
        # build node properties
        for i in self.CLD.TGs:
            TG = i
            self.colourNodesByEcoli(TG)
        
        # build edge properties
        for TG1 in self.TG_connectivity.keys():
            for TG2 in self.TG_connectivity[TG1]:
                G.add_edge(TG1,
                           TG2)
        
        self.node_colour_values = [self.node_colour.get(node) for node in G.nodes()]
        
    def colourNodesByEcoli(self, TG):
        # set red as default
        self.node_colour[TG] = 'green'
        for contig in self.CLD.TGs_contigs[TG]:
            if self.CLD.contig_to_genus[contig] == 'Escherichia':
                # change to blue if E.coli in TG
                self.node_colour[TG] = 'red'
    
    def colourNodesByStatus(self, TG):
        if self.CLD.TG_status[TG] == 'bogus':
            self.node_colour[TG] = 'red'
        elif  self.CLD.TG_status[TG] == 'real':
            self.node_colour[TG] = 'green'
        elif self.CLD.TG_status[TG] == 'unresolved':
            self.node_colour[TG] = 'orange'
        
    def Plot(self, outfile, outfmt, label_font_size, draw_network_labels):
        # initialise graph
        G = nx.Graph()
        
        # build network data
        self.buildNetworkData(G)
        
        # Build network plot
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        #pos= nx.spring_layout(G,iterations=500)
        pos= nx.fruchterman_reingold_layout(G, k=0.03)
        #pos= nx.circular_layout(G)
        
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               alpha=0.7,
                               node_color = self.node_colour_values,
                               with_labels=True,
                               node_size=150)
        
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = "#C3C3C3", #E1E1E1",
                               alpha=0.5)
        
        # draw network labels
        if draw_network_labels:
            nx.draw_networkx_labels(G,
                                    pos,
                                    font_size=label_font_size,
                                    font_color= 'black')
        
        # set axis properties
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')
        
        # write to file
        plt.savefig("%s" % (outfile),format="%s" % (outfmt))
        #plt.show()
            
            
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    print cmd
    args = shlex.split(cmd)
    p = subprocess.Popen(args) # shell=bash is not recommended. Only use when '>' must be in cmd. 
    return p.communicate()
    #p = Popen(cmd.split(' '), stdout=PIPE)
    #return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    NP = NetworkPlot(args.coverage_links_file)
    NP.Wrapper(args.outfile,
               args.outfmt,
               args.label_font_size,
               args.draw_network_labels)
                


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('coverage_links_file', help="")
    parser.add_argument('-o','--outfile', default='network.png',help="")
    parser.add_argument('-of','--outfmt', default='png',help="")
    parser.add_argument('-lfs','--label_font_size', type=int, default=12,help="")
    parser.add_argument('-dnl','--draw_network_labels', default=False,help="Draw network labels: True or False")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
