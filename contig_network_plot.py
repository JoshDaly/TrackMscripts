#!/usr/bin/env python
###############################################################################
#
# __contig_network_plot__.py - description!
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
import math as math
#import os
#import errno
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq

# local imports
from cb2cols import Cb2Cols as CB2
from trackm_file_parser import BamMData

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NetworkPlot(object):
    def __init__(self, bamm_links_file, bamm_cov_file):
        self.BD     = BamMData(bamm_links_file,
                               bamm_cov_file)
        self.node_size                  = {}
        self.node_size_values           = []
        self.node_colour                = {}
        self.node_colour_values         = []
        self.edge_width_links           = []
        self.max_cov                    = 0
        # Colour Brewer
        cb2                             = CB2()
        col_set                         = "qualSet1"
        col_set_gradient                = "seqReds"
        self.ColBrewColours             = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient     = cb2.maps[col_set_gradient].values()[0:10]
        
    def Plot(self, outfile, outfmt, label_font_size):
        # initialise graph
        G = nx.Graph()
        
        # build network data
        self.buildNetworkData(G)
        
        # Build network plot
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        #pos= nx.spring_layout(G,k=1,iterations=500)
        pos= nx.circular_layout(G)
        
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               node_size= self.node_size_values,
                               node_color = self.node_colour_values, 
                               alpha=1,
                               with_labels=True)
        
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = "#E1E1E1",
                               alpha=0.5,
                               width=self.edge_width_links)
        
        # draw network labels
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
    
    def buildNetworkData(self, G):
        # get nodes
        nodes = self.BD.bamm_cov_data.keys()
        
        # add nodes to graph 
        G.add_nodes_from(nodes)
        
        # set node colours using coverage
        self.getMaxCoverage()
        
        # build node properties
        for contig in self.BD.bamm_cov_data.keys():
            
            # Size nodes by contig len
            self.node_size[contig] = self.BD.contig_len[contig] 
        
            # Colour nodes coverage
            self.colourNodesByCoverage(contig)
            
        self.node_size_values   = [float(math.sqrt(int(self.node_size.get(node)))) for node in G.nodes()]  
        self.node_colour_values = [self.node_colour.get(node) for node in G.nodes()]
        
        # build edge properties
        for contig1 in self.BD.bamm_links_data.keys():
            for contig2 in self.BD.bamm_links_data[contig1]:
                G.add_edge(contig1,
                           contig2,
                           capacity = self.BD.bamm_links_data[contig1][contig2])
        
        # Edgewidth is # of transfer groups
        for (u,v,d) in G.edges(data=True):
            self.edge_width_links.append(math.sqrt(int(str(d).split(':')[-1].split()[0].split('}')[0])))
            
    def getMaxCoverage(self):
        covs = []
        for contig in self.BD.bamm_cov_data.keys():
            for cov in self.BD.bamm_cov_data[contig]:
                covs.append(cov)
            
        # get max coverage
        self.max_cov = max(covs)
        
    def colourNodesByCoverage(self, contig):
        cov_mean = self.averageArray(self.BD.bamm_cov_data[contig])
        adjusted_coverage = float(cov_mean)/self.max_cov
        
        print adjusted_coverage
            
        if adjusted_coverage >= 0 and adjusted_coverage < 0.2:
            self.node_colour[contig] = self.colBrewColoursGradient[4]
        elif adjusted_coverage >= 0.2 and adjusted_coverage < 0.4:
            self.node_colour[contig] = self.colBrewColoursGradient[5]
        elif adjusted_coverage >= 0.4 and adjusted_coverage < 0.6: 
            self.node_colour[contig] = self.colBrewColoursGradient[6]
        elif adjusted_coverage >= 0.6 and adjusted_coverage < 0.8: 
            self.node_colour[contig] = self.colBrewColoursGradient[7]
        elif adjusted_coverage >= 0.8 and adjusted_coverage <= 1: 
            self.node_colour[contig] = self.colBrewColoursGradient[8]
        
        
    def averageArray(self, array):
        a = np.array(array)
        return np.mean(a)
            
        

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
    NP = NetworkPlot(args.bamm_links_file,
                    args.bamm_cov_file)
    NP.Plot(args.outfile,
            args.outfmt,
            args.label_font_size)
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('bamm_links_file', help="")
    parser.add_argument('bamm_cov_file', help="")
    parser.add_argument('-o','--outfile', default='network.png',help="")
    parser.add_argument('-of','--outfmt', default='png',help="")
    parser.add_argument('-lfs','--label_font_size', type=int, default=12,help="")
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
