#!/usr/bin/env python
###############################################################################
#
# __analyse_crispr_spacer_hits__.py - description!
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
__copyright__ = "Copyright 2014"
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
import networkx as nx
import matplotlib.pyplot as plt
#import os
#import errno
#import numpy as np
#np.seterr(all='raise')

# local imports
import trackm_file_parser as TFP
from cb2cols import Cb2Cols as CB2

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CrisprData(object):
    def __init__(self, blast_file, group_file, hit_file, taxonomy_file):
        self.BD                 = TFP.BlastFileParser(blast_file)
        self.GD                 = TFP.GroupData(group_file)
        self.HD                 = TFP.HitData(hit_file, taxonomy_file)
        self.transfer_groups    = {}
        self.crispr_data        = {}
        self.spacer_to_genus    = {}
        
    def wrapper(self):
        # loop through blast data
        for spacer in self.BD.blast_data.keys():
            
            # get spacer taxon data
            #gid     = spacer.split("_")[0]
            #genus   = self.HD.gid_to_genus[gid]
            #self.spacer_to_genus[spacer] = genus
            
            for pidsqid in self.BD.blast_data[spacer]:
                transfer_group = self.whichTransferGroup(pidsqid)
                
                # build transfer group data
                self.addTransferGroup(spacer, transfer_group) 
                
                # build crispr data
                genus1 = self.HD.genus[pidsqid][0]
                genus2 = self.HD.genus[pidsqid][1]
                
                self.addGenus(spacer, genus1)
                self.addGenus(spacer, genus2)
            
        self.printOutData()
    
    def addGenus(self, spacer, genus):
        try:
            self.crispr_data[spacer][genus] = 1 
        except KeyError:
            self.crispr_data[spacer] = {genus:1}
                
    def printOutData(self):
        for spacer in self.crispr_data.keys():
            genuses = ''
            for genus in self.crispr_data[spacer]:
                genuses += "%s " % genus
            print "%s\t%s" % (spacer, genuses)
        
        
        #for spacer in self.transfer_groups.keys():
            
        #    for TG in self.transfer_groups[spacer]:
        #        print "%s\t%s" % (spacer,
        #                          TG)
                
    def addTransferGroup(self, spacer, transfer_group):
        try:
            self.transfer_groups[spacer][transfer_group] += 1 
        except KeyError:
            try:
                self.transfer_groups[spacer][transfer_group] = 1
            except KeyError:
                self.transfer_groups[spacer] = {transfer_group:1}  
    
    def whichTransferGroup(self, pidsqid):
        return self.GD.group_membership[pidsqid]

class Plot(object):
    def __init__(self, crispr_data):
        self.CD                         = crispr_data
        self.phylumDict                 = {}
        self.phylumCols                 = {}
        self.habitatDict                = {}
        self.habitatCols                = {}
        self.nodes_info                 = {}
        self.node_size                  = {}
        self.node_size_values           = []
        self.phylum_colour_values       = []
        self.habitat_colour_values      = []
        # ColorBrewer colours
        cb2                         = CB2()
        col_set                     = "qualSet1"
        col_set_gradient            = "seqReds"
        self.ColBrewColours         = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient = cb2.maps[col_set_gradient].values()[0:10]
    
    def networkXPlot(self, colour_node_by, outfile, outfmt, label_font_size):
        # initialise graph
        G = nx.Graph()
        
        # build graph data
        self.buildNetworkData(G)
        
        # Build network plot
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        pos= nx.spring_layout(G,iterations=500)
        
        # draw network nodes
        self.drawNetworkNodes(colour_node_by,
                              pos,
                              G)
        
        # draw network edges 
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = "#E1E1E1")
        
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
        outfile = "%s.%s" % (outfile,
                             outfmt)
        plt.savefig("%s" % (outfile),format="%s" % (outfmt))
                               
        
    def drawNetworkNodes(self, colour_node_by, pos, G):
        if colour_node_by.lower() == "phylum": 
            nx.draw_networkx_nodes(G,
                                   pos,
                                   linewidths=0,
                                   node_size= self.node_size_values,
                                   node_color = self.phylum_colour_values, 
                                   alpha=0.7,
                                   with_labels=True)
            
        elif colour_node_by.lower() == "habitat":
            nx.draw_networkx_nodes(G,
                                   pos,
                                   linewidths=0,
                                   node_size= self.node_size_values,
                                   node_color = self.habitat_colour_values, 
                                   alpha=0.7,
                                   with_labels=True)
            
        else:
            print "Please select node colour scheme: phylum or habitat"
            sys.exit()
        
    def buildNetworkData(self, G):
        # create list of nodes
        nodes = self.returnNodes()
        
        # add nodes
        G.add_nodes_from(nodes)
        
        # build node properties
        for node in nodes:
        
            # Size nodes by no. of genomes
            self.node_size[node] = len(self.nodes_info[node])*100
        
            # Colour nodes by phylum
            self.colourNodesByPhylum(node)
            
            # colour nodes by habitat
            #self.colourNodesByHabitat(node)
            
        # build node data
        self.node_size_values = [self.node_size.get(node) for node in G.nodes()]
        self.phylum_colour_values   = [self.phylumCols.get(node) for node in G.nodes()]
        self.habitat_colour_values   = [self.habitatCols.get(node) for node in G.nodes()]
        
        # build edge properties
        for node1 in self.CD.BD.blast_data.keys():
            for pidsqid in self.CD.BD.blast_data[node1]:
                node2 = self.CD.HD.genus[pidsqid][0]
                node3 = self.CD.HD.genus[pidsqid][1]
                G.add_edge(node1,
                           node2)
                G.add_edge(node1,
                           node3)
        
            
    def colourNodesByHabitat(self,node):
        try:
            self.habitatCols[node] = self.habitatDict[self.CD.HD.habitat_homog_data[node]]
        except KeyError:
            if self.CD.HD.habitat_homog_data[node] == "Gastrointestinal tract":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] =  self.ColBrewColours[0] #"#ea00ff"
                self.habitatCols[node] = self.ColBrewColours[0] #"#ea00ff"

            elif self.CD.HD.habitat_homog_data[node] == "Eye":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[1] #"#ffb700"
                self.habitatCols[node] = self.ColBrewColours[1] #"#ffb700"
                
            elif self.CD.HD.habitat_homog_data[node] == "Oral":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[2] #"#0047ff"
                self.habitatCols[node] = self.ColBrewColours[2] #"#0047ff"

            elif self.CD.HD.habitat_homog_data[node] == "Airways":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[3] #"#14ff00"
                self.habitatCols[node] = self.ColBrewColours[3] #"#14ff00"

            elif self.CD.HD.habitat_homog_data[node] == "internal_organs":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[4] #"#6600CC"
                self.habitatCols[node] = self.ColBrewColours[4] #"#6600CC"

            elif self.CD.HD.habitat_homog_data[node] == "Blood":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[5] #"#ffff00"
                self.habitatCols[node] = self.ColBrewColours[5] #"#ffff00"

            elif self.CD.HD.habitat_homog_data[node] == "skin":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[6] #"#006600"
                self.habitatCols[node] = self.ColBrewColours[6] #"#0080ff"

            elif self.CD.HD.habitat_homog_data[node] == "Urogenital tract":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[7] #"#00e0ff"
                self.habitatCols[node] = self.ColBrewColours[7] #"#00e0ff"
                
            elif self.CD.HD.habitat_homog_data[node] == "Plant":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = self.ColBrewColours[8] #"#00e0ff"
                self.habitatCols[node] = self.ColBrewColours[8] #"#00e0ff"

            elif self.CD.HD.habitat_homog_data[node] == "Ear":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = '#00CCCC' #self.ColBrewColours[8] #"#ff1e00"
                self.habitatCols[node] = '#00CCCC' #self.ColBrewColours[8] #"#ff1e00"
            
            elif self.CD.HD.habitat_homog_data[node] == "mixed":
                self.habitatDict[self.CD.HD.habitat_homog_data[node]] = '#C9C9C9' #self.ColBrewColours[8] #"#ff1e00"
                self.habitatCols[node] = '#C9C9C9' #self.ColBrewColours[8] #"#ff1e00
            
    def colourNodesByPhylum(self, node):
        try:
            self.phylumCols[node] = self.phylumDict[self.CD.HD.genus_to_phylum[node]]
        except KeyError:
            try:
                if self.CD.HD.genus_to_phylum[node] == "p__Bacteroidetes":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] =  self.ColBrewColours[0] #"#ea00ff"
                    self.phylumCols[node] = self.ColBrewColours[0] #"#ea00ff"
    
                elif self.CD.HD.genus_to_phylum[node] == "p__Lentisphaerae":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[1] #"#ffb700"
                    self.phylumCols[node] = self.ColBrewColours[1] #"#ffb700"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Firmicutes":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[2] #"#0047ff"
                    self.phylumCols[node] = self.ColBrewColours[2] #"#0047ff"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Spirochaetes":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[3] #"#14ff00"
                    self.phylumCols[node] = self.ColBrewColours[3] #"#14ff00"   
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Synergistetes":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[4] #"#6600CC"
                    self.phylumCols[node] = self.ColBrewColours[4] #"#6600CC"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Actinobacteria":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[5] #"#ffff00"
                    self.phylumCols[node] = self.ColBrewColours[5] #"#ffff00"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Tenericutes":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[6] #"#006600"
                    self.phylumCols[node] = self.ColBrewColours[6] #"#0080ff"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Fusobacteria":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[7] #"#00e0ff"
                    self.phylumCols[node] = self.ColBrewColours[7] #"#00e0ff"
                    
                elif self.CD.HD.genus_to_phylum[node] == "p__Epsilonmicrobia": 
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = self.ColBrewColours[8]
                    self.phylumCols[node] = self.ColBrewColours[8] 
                
                elif self.CD.HD.genus_to_phylum[node] == "p__Proteobacteria":
                    self.phylumDict[self.CD.HD.genus_to_phylum[node]] = "#01B5B5" # self.ColBrewColours[8] #"#ff1e00"
                    self.phylumCols[node] = "#01B5B5" #self.ColBrewColours[8] #"#ff1e00"
            except KeyError:
                self.phylumCols[node] = "#5E5E5E"
                
    def returnNodes(self):
        self.nodes_info = {}
        
        # loop through blast data
        for spacer in self.CD.BD.blast_data.keys():
            try:
                self.nodes_info[spacer]['a'] = 1
            except KeyError:
                self.nodes_info[spacer] = {'a':1}
                
            for pidsqid in self.CD.BD.blast_data[spacer]:
                genus1 = self.CD.HD.genus[pidsqid][0] 
                genus2 = self.CD.HD.genus[pidsqid][1]
                gid1   = self.CD.HD.pidsqid_to_gid[pidsqid][0] 
                gid2   = self.CD.HD.pidsqid_to_gid[pidsqid][0]
                
                try:
                    self.nodes_info[genus1][gid1] = 1
                except KeyError:
                    self.nodes_info[genus1] = {gid1:1}
                try:
                    self.nodes_info[genus2][gid2]
                except KeyError:
                    self.nodes_info[genus2] = {gid2:1}
        
        return self.nodes_info.keys()
                    
        
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

    CD = CrisprData(args.blast_file,
                    args.group_file,
                    args.hit_file,
                    args.taxonomy_file 
                    )
    CD.wrapper()
    
    P = Plot(CD)
    P.networkXPlot(args.colour_node_by,
                   args.out_file,
                   args.outfmt,
                   args.label_font_size
                   )


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file', help="File containing CRISPR Blast output")
    parser.add_argument('group_file', help="File containing transfer groups")
    parser.add_argument('hit_file', help="File containing hit data")
    parser.add_argument('taxonomy_file', help="File containing Phils updated taxonomy")
    parser.add_argument('-colour_node_by','--colour_node_by', default='phylum', help="Colour node by phylum or habitat")
    parser.add_argument('-outfmt','--outfmt', default='png',help="Set the output file format. Default: png.")
    parser.add_argument('-o','--out_file', default='Network_plot',help="Set output file, and directory.")
    parser.add_argument('-lfs','--label_font_size', default=8, type=int, help="Node label font size.")
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
