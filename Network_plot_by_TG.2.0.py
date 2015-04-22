#!/usr/bin/env python
###############################################################################
#
# __Network_plot_by_TG.2.0__.py - description!
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
__email__ = ""
__status__ = "Development"

###############################################################################

# system imports
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
from collections import Counter
import math

# modules to be loaded
import numpy as np
np.seterr(all='raise')
import networkx as nx
import matplotlib.pyplot as plt

# local imports
from cb2cols import Cb2Cols as CB2
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class View(object):
    def __init__(self, hit_data, transfer_group_file, contam_pidsqids_file):
        self.TG                             = TFP.GroupData(transfer_group_file)
        self.HD                             = TFP.HitData(hit_data)
        self.clean_analysis                 = contam_pidsqids_file
        if self.clean_analysis:
            self.CP                         = TFP.ContaminatedPidsqids(contam_pidsqids_file)
            
class ScatterPlot(object):
    def __init__(self):
        pass

class NetworkPlot(object):
    def __init__(self, hit_data, transfer_group_file, contam_pidsqids_file):
        #self.AD                             = TFP.AnnotationData(annotation_file)
        #self.TD                             = TFP.TaxonomyData(taxon_file)
        #self.KD                             = TFP.KmerData(kmer_file)
        self.TG                             = TFP.GroupData(transfer_group_file)
        self.HD                             = TFP.HitData(hit_data)
        self.clean_analysis                 = contam_pidsqids_file
        if self.clean_analysis:
            self.CP                         = TFP.ContaminatedPidsqids(contam_pidsqids_file)
        self.nodes                          = {}
        self.edges                          = {} 
        self.node_size                      = {}
        self.node_size_values               = []
        self.edgewidth                      = []
        self.nodeCols                       = {}
        self.nodeDict                       = {}
        self.node_colour_values             = []
        # ColorBrewer colours
        cb2                                 = CB2()
        col_set                             = "qualSet1"
        col_set_gradient                    = "seqReds"
        self.ColBrewColours                 = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient         = cb2.maps[col_set_gradient].values()[0:10]
        
        
    def wrapper(self, outfile, outfmt, data_type, centrality_type, node_repulsion):
        # ready hit data for network plot
        self.prepareHitData(data_type)
        
        # make network plot
        self.networkPlot(outfile, outfmt, centrality_type, data_type, node_repulsion)
        
    def prepareHitData(self, data_type):
        # loop through pidsqids
        for pidsqid in self.HD.hit_data.keys():
            
            # check if pidsqid passes criteria
            if self.shallYouPass(pidsqid, data_type):
                
                # build network data
                self.buildHitData(pidsqid, data_type)
            
    def shallYouPass(self, pidsqid, data_type):
        # shall you pass khazad dum, Balrog? 
        status = True
        if data_type:
            if data_type == 'both' or data_type == 'habitat':
                if self.clean_analysis:
                    # remove contaminated pidsqids
                    status = self.checkIfPidsqidContaminated(pidsqid)
                return status
            elif data_type == 'intra':
                if self.HD.intra_or_inter[pidsqid] == 'intra':
                    if self.clean_analysis:
                        # remove contaminated pidsqids
                        status = self.checkIfPidsqidContaminated(pidsqid)
                    return status
            elif data_type == 'inter':
                if self.HD.intra_or_inter[pidsqid] == 'inter':
                    if self.clean_analysis:
                        # remove contaminated pidsqids
                        status = self.checkIfPidsqidContaminated(pidsqid)
                    return status
            else:
                print "####################"
                print "Please select either inter, intra, both or habitat ie. -dt intra "
                print "####################"
                sys.exit()
        else:
            if self.clean_analysis:
                # remove contaminated pidsqids
                status = self.checkIfPidsqidContaminated(pidsqid)
            return status
        

    def checkIfPidsqidContaminated(self, pidsqid):
        if pidsqid in self.CP.contam_pidsqids:
            return False
        else:
            return True
        
    def buildHitData(self, pidsqid, data_type):
        # get pidsqid data
        try:
            transfer_group  = self.TG.group_membership[pidsqid]
            if data_type == 'habitat':
                # nodes = habitats!
                node1       = self.HD.habitat[pidsqid][0] 
                node2       = self.HD.habitat[pidsqid][1]
            else:
                # nodes = genus!
                node1       = self.HD.genus[pidsqid][0]
                node2       = self.HD.genus[pidsqid][1]
            gid1            = self.HD.pidsqid_to_gid[pidsqid][0]
            gid2            = self.HD.pidsqid_to_gid[pidsqid][1]
            
            # add node
            self.addNode(node1, gid1)
            self.addNode(node2, gid2)
            
            # add edge, edges weighted by number of TGs
            self.addEdge(node1, node2, transfer_group)
            
        except KeyError:
            pass
        
    def colourNodesWrapper(self, node, data_type):
        if data_type == 'habitat':
            self.colourNodesByHabitat(node)
        else:
            self.colourNodesByPhylum(node)
    
    def colourNodesByHabitat(self, habitat):
        if habitat not in self.nodeCols:
            if habitat == 'Blood':
                self.nodeDict[habitat] = self.ColBrewColours[0]
                self.nodeCols[habitat] = self.ColBrewColours[0]
                
            elif habitat == 'Urogenital tract':
                self.nodeDict[habitat] = self.ColBrewColours[1]
                self.nodeCols[habitat] = self.ColBrewColours[1]
            
            elif habitat == 'plant':
                self.nodeDict[habitat] = self.ColBrewColours[2]
                self.nodeCols[habitat] = self.ColBrewColours[2]
                
            elif habitat == 'Gastrointestinal tract':
                self.nodeDict[habitat] = self.ColBrewColours[3]
                self.nodeCols[habitat] = self.ColBrewColours[3]
            
            elif habitat == 'Oral':
                self.nodeDict[habitat] = self.ColBrewColours[4]
                self.nodeCols[habitat] = self.ColBrewColours[4]
                
            elif habitat == 'Airways':
                self.nodeDict[habitat] = self.ColBrewColours[5]
                self.nodeCols[habitat] = self.ColBrewColours[5]
                
            elif habitat == 'skin':
                self.nodeDict[habitat] = self.ColBrewColours[6]
                self.nodeCols[habitat] = self.ColBrewColours[6]
                
            elif habitat == 'internal_organs':
                self.nodeDict[habitat] = self.ColBrewColours[7]
                self.nodeCols[habitat] = self.ColBrewColours[7]
                
            elif habitat == 'Ear':
                self.nodeDict[habitat] = self.ColBrewColours[8]
                self.nodeCols[habitat] = self.ColBrewColours[8]
                
            elif habitat == 'Eye':
                self.nodeDict[habitat] = "#01B5B5"
                self.nodeCols[habitat] = "#01B5B5" 
    
    def colourNodesByPhylum(self, genus):
        if genus not in self.nodeCols:
            if self.HD.genus_to_phylum[genus] == "p__Bacteroidetes":
                self.nodeDict["p__Bacteroidetes"] =  self.ColBrewColours[0] #"#ea00ff"
                self.nodeCols[genus] = self.ColBrewColours[0] #"#ea00ff"

            elif self.HD.genus_to_phylum[genus] == "p__Lentisphaerae":
                self.nodeDict["p__Lentisphaerae"] = self.ColBrewColours[1] #"#ffb700"
                self.nodeCols[genus] = self.ColBrewColours[1] #"#ffb700"
                
            elif self.HD.genus_to_phylum[genus] == "p__Firmicutes":
                self.nodeDict["p__Firmicutes"] = self.ColBrewColours[2] #"#0047ff"
                self.nodeCols[genus] = self.ColBrewColours[2] #"#0047ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Spirochaetes":
                self.nodeDict["p__Spirochaetes"] = self.ColBrewColours[3] #"#14ff00"
                self.nodeCols[genus] = self.ColBrewColours[3] #"#14ff00"   
                
            elif self.HD.genus_to_phylum[genus] == "p__Synergistetes":
                self.nodeDict["p__Synergistetes"] = self.ColBrewColours[4] #"#6600CC"
                self.nodeCols[genus] = self.ColBrewColours[4] #"#6600CC"
                
            elif self.HD.genus_to_phylum[genus] == "p__Actinobacteria":
                self.nodeDict["p__Actinobacteria"] = self.ColBrewColours[5] #"#ffff00"
                self.nodeCols[genus] = self.ColBrewColours[5] #"#ffff00"
                
            elif self.HD.genus_to_phylum[genus] == "p__Tenericutes":
                self.nodeDict["p__Tenericutes"] = self.ColBrewColours[6] #"#006600"
                self.nodeCols[genus] = self.ColBrewColours[6] #"#0080ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Fusobacteria":
                self.nodeDict["p__Fusobacteria"] = self.ColBrewColours[7] #"#00e0ff"
                self.nodeCols[genus] = self.ColBrewColours[7] #"#00e0ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Epsilonmicrobia": 
                self.nodeDict["p__Epsilonmicrobia"] = self.ColBrewColours[8]
                self.nodeCols[genus] = self.ColBrewColours[8] 
            
            elif self.HD.genus_to_phylum[genus] == "p__Proteobacteria":
                self.nodeDict["p__Proteobacteria"] = "#01B5B5" # self.ColBrewColours[8] #"#ff1e00"
                self.nodeCols[genus] = "#01B5B5" #self.ColBrewColours[8] #"#ff1e00"
            else:
                print "%s not in a recognised phylum" % genus
    
    def addNode(self, node, gid):
        try:
            self.nodes[node][gid] = 1 
        except KeyError:
            self.nodes[node]= {gid:1}
            
    def addEdge(self, node1, node2, transfer_group):
        # make it a one side pyramid of data
        if node1 > node2: 
            try:
                self.edges[node1][node2][transfer_group] = 1
            except KeyError:
                try:
                    self.edges[node1][node2] = {transfer_group:1}
                except KeyError:
                    self.edges[node1] = {node2:{transfer_group:1}}
        else:
            try:
                self.edges[node2][node1][transfer_group] = 1
            except KeyError:
                try:
                    self.edges[node2][node1] = {transfer_group:1}
                except KeyError:
                    self.edges[node2] = {node1:{transfer_group:1}}
    
    def networkPlot(self, outfile, outfmt, centrality_type, data_type, node_repulsion):
        # initialise network plot
        G = nx.Graph()
        
        # build network data
        self.buildNetworkData(G, data_type)
        
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        pos= nx.spring_layout(G,k=node_repulsion,iterations=500)
        
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               node_size= self.node_size_values,
                               node_color = self.node_colour_values,
                               alpha=0.7,
                               with_labels=True)
        # draw network labels
        nx.draw_networkx_labels(G,
                                pos,
                                font_size=12,
                                font_color= 'black')
        # draw network edges
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = '#E5E5E5',
                               width=self.edgewidth)
        
        # set axis properties
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')
        
        
        # calculate centrality
        self.calculateCentrality(G, centrality_type, outfile)
        
        # adjacency matrix
        self.createAdjacencyMatrix(G, outfile)
        
        # write to file
        plt.savefig("%s.%s" % (outfile, outfmt),format="%s" % (outfmt), dpi = 300)
        
    def createAdjacencyMatrix(self, G, outfile):
        
        matrix = []
        
        # produce adjacency matrix
        A = nx.to_dict_of_dicts(G)
        
        # nodes 
        nodes = A.keys()
        
        # populate matrix with zeros
        self.zeroArray(matrix, nodes)
        
        # loop through data and populate matrix
        for i in range(len(nodes)):
            for j in range(i, len(nodes)):
                
                # get nodes names 
                node1 = nodes[i]
                node2 = nodes[j]
                
                # calculate hit
                try:
                    hits = str(A[node1][node2]['capacity'])
                except KeyError:
                    hits = '0'
                
                # add to matrix
                matrix[i][j] = hits
        
        # insert node names to top of array
        nodes.insert(0, 'node')
        matrix.insert(0, nodes)
        
        # insert nodes name to start of each array
        for i in range(len(nodes)):
            
            # dont consider first element
            if i != 0:
                
                # insert node name to array
                matrix[i].insert(0, nodes[i])
        
        # matrix outfile
        matrix_outfile = open("%s.adjacency_matrix.txt" % outfile, 'w')
        
        # write matrix to file
        for i in matrix:
            string_to_write = "\t".join(i)
            matrix_outfile.write("%s\n" % string_to_write)

    def zeroArray(self, matrix, array):
        for i in range(len(array)):
            zeros = []
            for v in range(len(array)):
                zeros.append('0')
            matrix.append(zeros)
    
    def calculateCentrality(self, G, centrality_type, outfile):
        if centrality_type == 'degree':
            
            # calculate degree centrality
            centrality = nx.degree_centrality(G)
            
            # create output file
            centrality_outfile = open("%s.centrality.csv" % outfile, 'w')
            
            # loop through nodes
            for node in centrality.keys():
                
                # determine centrality
                degree_centrality = centrality[node]
                
                # write to file
                centrality_outfile.write("%s\t%s\n" % (node, degree_centrality))
            
            # close file
            centrality_outfile.close()
            
        elif centrality_type == 'eigenvector':
            
            # calculate eigenvector centrality
            eigenvector_centrality = nx.eigenvector_centrality(G)
            
            # create output file
            eigenvector_centrality_outfile = open("%s.eigenvector_centrality.csv" % outfile, 'w')
            
            # loop through nodes
            for node in centrality.keys():
                
                # determine centrality
                eigen_centrality  = eigenvector_centrality[node]
                
                # write to file
                eigenvector_centrality_outfile.write("%s\t%s\n" % (node, eigen_centrality))
            
            # close file
            eigenvector_centrality_outfile.close()
                
    def buildNetworkData(self, G, data_type):
        """build data to be used to create network plot"""

        # get list of genus nodes
        nodes = self.nodes.keys()
        
        # add nodes 
        G.add_nodes_from(nodes)
        
        # build node properties
        for node in nodes:
            
            # Size nodes by no. of genomes
            self.node_size[node] = len(self.nodes[node])*100
            
            # colour nodes by phylum
            self.colourNodesWrapper(node, data_type)
            
        # build node data
        self.node_size_values = [self.node_size.get(node) for node in G.nodes()]
        self.node_colour_values   = [self.nodeCols.get(node) for node in G.nodes()]
        
        # build edge properties
        for node1 in self.edges.keys():
            for node2 in self.edges[node1]:
                
                G.add_edge(node1,
                           node2, 
                           capacity = len(self.edges[node1][node2]), # aka edge width 
                           weight = len(self.edges[node1][node2]))
                
        # Edgewidth is # of transfer groups
        for (u,v,d) in G.edges(data=True):
            self.edgewidth.append(math.sqrt(int(str(d).split(':')[-1].split()[0].split('}')[0])))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    NP = NetworkPlot(args.hit_data,
                     args.transfer_group_file,
                     args.contam_pidsqids_file)
    NP.wrapper(args.outfile,
               args.outfmt,
               args.data_type,
               args.centrality_type,
               args.node_repulsion)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hit_data', help="File containing hit data")
    parser.add_argument('transfer_group_file', help="File containing transfer groups")
    parser.add_argument('-cpf','--contam_pidsqids_file', default=False, help="File containing a list of the contaminated pidsqids")
    parser.add_argument('-o','--outfile', default='network_plot', help="Output file name prefix.")
    parser.add_argument('-of','--outfmt', default='png', help="format of network plot. Default = png.")
    parser.add_argument('-k','--node_repulsion', default=0.15, type=float, help="Set the node repulsion value 0-1. Default = 0.15.")
    parser.add_argument('-dt','--data_type', default='all', help="Limit transfers to: inter, intra or both[default].")
    parser.add_argument('-ct','--centrality_type', default='degree', help="Determine node centrality and write to file. degree (default) or eigenvector.")
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
