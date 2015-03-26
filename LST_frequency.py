#!/usr/bin/env python
###############################################################################
#
# __LST_frequency__.py - description!
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
#import os
#import errno
#import glob
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq
#import matplotlib.pyplot as plt
#import networkx as nx

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Data(object):
    def __init__(self, name):
        self.__name  = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

class LSTFrequency(object):
    def __init__(self, hit_data_file, taxonomy_file, pairs_file):
        self.HD             = TFP.HitData(hit_data_file)
        self.TD             = TFP.TaxonomyData(taxonomy_file)
        self.PD             = TFP.PairsFileData(pairs_file)
        self.ani_data       = {}
        self.gid_ani_group  = {}
        self.taxon_data     = {}
        self.pidsqid_count  = {}
        self.collated_data  = {}
    
    def wrapper(self):
        # gather taxonomy data
        self.gatherTaxonomyData()
        
        # make ANI groups
        self.calculateComparisons()
        
        # are there any times when the genus == genus for the hitdata
        self.getPidsqidCounts()
        
        # combine pidsqid counts with no. of comparisons
        self.collateDataWrapper()
        
        for rank in self.collated_data.keys():
            print rank
            for rank_value in self.collated_data[rank]:
                comparisons = int(self.collated_data[rank][rank_value][1])
                pidsqids    = int(self.collated_data[rank][rank_value][0])
                pidsqid_per_100 = float(pidsqids)/comparisons * 100
                print rank_value, str(pidsqids), str(comparisons), str(pidsqid_per_100)
        
    def collateDataWrapper(self):
        for gid in self.TD.gids.keys():
            try:
                pidsqid_count   = len(self.pidsqid_count[gid])
            except KeyError:
                pidsqid_count   = 0
            comparisons     = len(self.PD.pairs_data[gid])
            self.collateData(gid, pidsqid_count, comparisons)
            
    def collateData(self, gid, pidsqid_count, comparisons):
        for rank in self.taxon_data[gid]:
            
            rank_data = self.taxon_data[gid][rank]
            try:
                self.collated_data[rank][rank_data][0] += pidsqid_count
                self.collated_data[rank][rank_data][1] += comparisons
            except KeyError:
                try:
                    self.collated_data[rank][rank_data] = [pidsqid_count, comparisons]
                except KeyError:
                    self.collated_data[rank] = {rank_data:[pidsqid_count, comparisons]}
                    
    
    def getPidsqidCounts(self):
        for pidsqid in self.HD.pidsqid_to_gid.keys():
            gid1 = self.HD.pidsqid_to_gid[pidsqid][0]
            gid2 = self.HD.pidsqid_to_gid[pidsqid][1]
            self.addPidsqid(pidsqid, gid1)
            self.addPidsqid(pidsqid, gid2)
            # check genus
            #genus1 = self.getGenus(gid1)
            #genus2 = self.getGenus(gid2)
    
    def addPidsqid(self, pidsqid, gid):
        try:
            self.pidsqid_count[gid][pidsqid] = 1 
        except KeyError: 
            self.pidsqid_count[gid] = {pidsqid:1}
            
            
    def getGenus(self, gid):
        if self.TD.taxon_genus[gid] == '':
            return self.TD.taxon_organism[gid]
        else:
            return self.TD.taxon_genus[gid]
    
    def gatherTaxonomyData(self):
        for gid in self.TD.gids.keys():
            _genus  = self.TD.taxon_genus[gid]
            if _genus == '':
                _genus = self.TD.taxon_organism[gid]
            _family = self.TD.taxon_family[gid]
            _order  = self.TD.taxon_order[gid]
            _class  = self.TD.taxon_class[gid]
            _phylum = self.TD.taxon_phylum[gid]
            self.addTaxonData(gid, 'genus', _genus)
            self.addTaxonData(gid, 'family', _family)
            self.addTaxonData(gid, 'order', _order)
            self.addTaxonData(gid, 'class', _class)
            self.addTaxonData(gid, 'phylum', _phylum)
            
    def addTaxonData(self, gid, rank, rank_value):
        # assign rank value to rank i.e genus -> escherichia

        try:
            self.taxon_data[gid][rank] = rank_value
        except KeyError:
            self.taxon_data[gid] = {rank:rank_value}
    
    def calculateComparisons(self):
        gids = self.TD.gids.keys()
        for i in range(len(gids)):
            
            gid1 = gids[i]
            
            # initialise gid
            self.ani_data[gid1] = {}
            
            for j in range(i+1, len(gids)):
                
                gid2 = gids[j]
                ani  = float(self.PD.pairs_data[gid1][gid2])
                
                if ani >= 95:
                    # within ANI cutoff
                    try:
                        self.ani_data[gid1][gid2] = 1 
                    except KeyError:
                        self.ani_data[gid1] = {gid2:1} 
        
        #self.getConnectedNodes()

    def getConnectedNodes(self):
        node_dict = {}
        
        nodes = set()

        # build trees from ANI data

        regions_to_nodes= {}
        
        for region in self.ani_data.keys():
            if region not in regions_to_nodes:
                D = Data(region)
                regions_to_nodes[region] = D
                nodes |= {D}
            for linked_region in self.ani_data[region].keys():
                if linked_region not in regions_to_nodes:
                    D = Data(linked_region)
                    regions_to_nodes[linked_region] = D
                    nodes |= {D}

        for region1 in self.ani_data.keys():
            for region2 in self.ani_data[region1].keys():
                regions_to_nodes[region1].add_link(regions_to_nodes[region2])
                
        # Find all the connected components
        number = 1

        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            
            names = ",".join(names)
            print "ani_group_%i\t%d\t%s" % (number,names_count, names)
            number += 1
                
                
    
    def connected_components(self,nodes):
        # List of connected components found. The order is random.
        result = []

        # Make a copy of the set, so we can modify it.
        nodes = set(nodes)

        nn = float(len(nodes))
        remains = nn

        # Iterate while we still have nodes to process.
        while nodes:

            # Get a random node and remove it from the global set.
            n = nodes.pop()

            # This set will contain the next group of nodes connected to each other.
            group = {n}

            # Build a queue with this node in it.
            queue = [n]

            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:

                # Consume the next item from the queue.
                n = queue.pop(0)

                # Fetch the neighbors.
                neighbors = n.links

                # Remove the neighbors we already visited.
                neighbors.difference_update(group)

                # Remove the remaining nodes from the global set.
                remains -= float(len(neighbors))
                #print "%f completed" % (1 - remains/nn)
                nodes.difference_update(neighbors)

                # Add them to the group of connected nodes.
                group.update(neighbors)

                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(group)

        # Return the list of groups.
        return result
                    

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
    LSTF = LSTFrequency(args.hit_data_file,
                        args.taxonomy_file,
                        args.pairs_file)
    LSTF.wrapper()
                
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hit_data_file', help="")
    parser.add_argument('taxonomy_file', help="")
    parser.add_argument('pairs_file', help="")
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
