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

# modules to be loaded
import numpy as np
np.seterr(all='raise')
import networkx as nx
import matplotlib.pyplot as plt

# local imports
from cb2cols import Cb2Cols as CB2
import trackm.fileparser as TFP

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
        self.genus_node_size                = {}
        self.genus_node_size_values         = []
        self.edgewidth                      = []
        self.phylumCols                     = {}
        self.phylumDict                     = {}
        self.phylum_colour_values           = []
        # ColorBrewer colours
        cb2                                 = CB2()
        col_set                             = "qualSet1"
        col_set_gradient                    = "seqReds"
        self.ColBrewColours                 = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient         = cb2.maps[col_set_gradient].values()[0:10]
        
        
    def wrapper(self, outfile, outfmt):
        # ready hit data for network plot
        self.prepareHitData()
        
        # make network plot
        self.networkPlot(outfile, outfmt)
        
    def prepareHitData(self):
        # loop through pidsqids
        for pidsqid in self.HD.hit_data.keys():
            
            # check if pidsqid passes criteria
            if self.shallYouPass(pidsqid):
                
                # build network data
                self.buildHitData(pidsqid)
            
    def shallYouPass(self, pidsqid):
        # shall you pass khazad dum, Balrog? 
        status = True
        if self.clean_analysis:
            # remove contaminated pidsqids
            status = self.checkIfPidsqidContaminated()
        
        return status

    def checkIfPidsqidContaminated(self, pidsqid):
        if pidsqid in self.CP.contam_pidsqids:
            return False
        else:
            return True
        
    def buildHitData(self, pidsqid):
        # get pidsqid data
        try:
            transfer_group  = self.TG.group_membership[pidsqid]
            genus1          = self.HD.genus[pidsqid][0]
            genus2          = self.HD.genus[pidsqid][1]
            gid1            = self.HD.pidsqid_to_gid[pidsqid][0]
            gid2            = self.HD.pidsqid_to_gid[pidsqid][1]
            
            # add node
            self.addGenusNode(genus1, gid1)
            self.addGenusNode(genus2, gid2)
            
            # add edge, edges weighted by number of TGs
            self.addEdge(genus1, genus2, transfer_group)
            
        except KeyError:
            pass
        
    def colourNodesByPhylum(self, genus):
        if genus not in self.phylumCols:
            if self.HD.genus_to_phylum[genus] == "p__Bacteroidetes":
                self.phylumDict["p__Bacteroidetes"] =  self.ColBrewColours[0] #"#ea00ff"
                self.phylumCols[genus] = self.ColBrewColours[0] #"#ea00ff"

            elif self.HD.genus_to_phylum[genus] == "p__Lentisphaerae":
                self.phylumDict["p__Lentisphaerae"] = self.ColBrewColours[1] #"#ffb700"
                self.phylumCols[genus] = self.ColBrewColours[1] #"#ffb700"
                
            elif self.HD.genus_to_phylum[genus] == "p__Firmicutes":
                self.phylumDict["p__Firmicutes"] = self.ColBrewColours[2] #"#0047ff"
                self.phylumCols[genus] = self.ColBrewColours[2] #"#0047ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Spirochaetes":
                self.phylumDict["p__Spirochaetes"] = self.ColBrewColours[3] #"#14ff00"
                self.phylumCols[genus] = self.ColBrewColours[3] #"#14ff00"   
                
            elif self.HD.genus_to_phylum[genus] == "p__Synergistetes":
                self.phylumDict["p__Synergistetes"] = self.ColBrewColours[4] #"#6600CC"
                self.phylumCols[genus] = self.ColBrewColours[4] #"#6600CC"
                
            elif self.HD.genus_to_phylum[genus] == "p__Actinobacteria":
                self.phylumDict["p__Actinobacteria"] = self.ColBrewColours[5] #"#ffff00"
                self.phylumCols[genus] = self.ColBrewColours[5] #"#ffff00"
                
            elif self.HD.genus_to_phylum[genus] == "p__Tenericutes":
                self.phylumDict["p__Tenericutes"] = self.ColBrewColours[6] #"#006600"
                self.phylumCols[genus] = self.ColBrewColours[6] #"#0080ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Fusobacteria":
                self.phylumDict["p__Fusobacteria"] = self.ColBrewColours[7] #"#00e0ff"
                self.phylumCols[genus] = self.ColBrewColours[7] #"#00e0ff"
                
            elif self.HD.genus_to_phylum[genus] == "p__Epsilonmicrobia": 
                self.phylumDict["p__Epsilonmicrobia"] = self.ColBrewColours[8]
                self.phylumCols[genus] = self.ColBrewColours[8] 
            
            elif self.HD.genus_to_phylum[genus] == "p__Proteobacteria":
                self.phylumDict["p__Proteobacteria"] = "#01B5B5" # self.ColBrewColours[8] #"#ff1e00"
                self.phylumCols[genus] = "#01B5B5" #self.ColBrewColours[8] #"#ff1e00"
            else:
                print "%s not in a recognised phylum" % genus
                
    
    
    def addGenusNode(self, genus, gid):
        try:
            self.nodes[genus][gid] = 1 
        except KeyError:
            self.nodes[genus]= {gid:1}
            
    def addEdge(self, genus1, genus2, transfer_group):
        # make it a one side pyramid of data
        if genus1 > genus2: 
            try:
                self.edges[genus1][genus2][transfer_group] = 1
            except KeyError:
                try:
                    self.edges[genus1][genus2] = {transfer_group:1}
                except KeyError:
                    self.edges[genus1] = {genus2:{transfer_group:1}}
        else:
            try:
                self.edges[genus2][genus1][transfer_group] = 1
            except KeyError:
                try:
                    self.edges[genus2][genus1] = {transfer_group:1}
                except KeyError:
                    self.edges[genus2] = {genus1:{transfer_group:1}}
    
    def networkPlot(self, outfile, outfmt):
        # initialise network plot
        G = nx.Graph()
        
        # build network data
        self.buildNetworkData(G)
        
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        pos= nx.spring_layout(G,k=0.15,iterations=500)
        
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               node_size= self.genus_node_size_values,
                               node_color = self.phylum_colour_values,
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
        
        # write to file
        plt.savefig("%s" % (outfile),format="%s" % (outfmt), dpi = 300)
        
    def buildNetworkData(self, G):
        """build data to be used to create network plot"""

        # get list of genus nodes
        genus_nodes = self.nodes.keys()
        
        # add nodes 
        G.add_nodes_from(genus_nodes)
        
        # build node properties
        for genus in genus_nodes:
            
            # Size nodes by no. of genomes
            self.genus_node_size[genus] = len(self.nodes[genus])*100
            
            # colour nodes by phylum
            self.colourNodesByPhylum(genus)
            
        # build node data
        self.genus_node_size_values = [self.genus_node_size.get(node) for node in G.nodes()]
        self.phylum_colour_values   = [self.phylumCols.get(node) for node in G.nodes()]
        
        # build edge properties
        for genus1 in self.edges.keys():
            for genus2 in self.edges[genus1]:
                
                G.add_edge(genus2,
                           genus1, 
                           capacity = len(self.edges[genus1][genus2]))
                           #weight = len(self.edges[genus1][genus2]))
                
        # Edgewidth is # of transfer groups
        for (u,v,d) in G.edges(data=True):
            self.edgewidth.append(int(str(d).split(':')[-1].split()[0].split('}')[0]))

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
               args.outfmt)
    
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
