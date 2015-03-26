#!/usr/bin/env python
###############################################################################
#
# __make_transfer_groups__.py - description!
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
from Bio import SeqIO
from Bio.Seq import Seq
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

class TransferGroupsNucmer(object):
    def __init__(self, hitdata):
        self.HD             = TFP.HitData(hitdata)
        self.NP             = TFP.NucMerParser()
        self.nucmer_data    = {}
        self.transfer_len   = {}
        
    def wrapper(self, nucmer_file, identity_thresh, len_overlap, fasta_file):
        # grab transfer length information from multifasta file
        self.parseFastaFile(fasta_file)
        
        # parse through nucmer file
        self.parseNucmerFile(nucmer_file, identity_thresh, len_overlap)
        
        # use nucmer data to form transfer groups
        self.makeTransferGroups()
        
        
    def parseFastaFile(self, fasta_file):
        for accession,sequence in SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta")).items():
            self.transfer_len[accession] = len(sequence.seq)
    
    def parseNucmerFile(self, nucmer_file, identity_thresh, len_overlap):
        with open(nucmer_file) as fh:
            for hit in self.NP.readNuc(fh):
                
                # check sequence len and identity 
                self.checkThresholds(hit[self.NP._ID_1],
                                     hit[self.NP._LEN_1], 
                                     hit[self.NP._ID_2],
                                     hit[self.NP._LEN_2], 
                                     hit[self.NP._IDENTITY], 
                                     identity_thresh, 
                                     len_overlap)
                    
                    

    def checkThresholds(self, id1, len1, id2, len2, identity, identity_thresh, len_overlap):
        # make percentage
        len_overlap = len_overlap/float(100)
        
        # check identity
        if identity >= identity_thresh:
            
            # check length overlap
            len_overlap_1 = len1/float(self.transfer_len[id1])
            len_overlap_2 = len2/float(self.transfer_len[id2])
            
            # overlap must exceed threshold for at least one of the sequences
            if len_overlap_1 >= len_overlap or len_overlap_2 >= len_overlap:
                
                self.addConnection(id1, id2)
            
        else:
            return False

    def addConnection(self, id1, id2):
        try:
            self.nucmer_data[id1][id2] = 1 
        except KeyError:
            self.nucmer_data[id1] = {id2:1}
        try:
            self.nucmer_data[id2][id1] = 1 
        except KeyError:
            self.nucmer_data[id2] = {id1:1}
    
    def makeTransferGroups(self):
        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_to_nodes= {}

        for region in self.nucmer_data.keys():
            if region not in regions_to_nodes:
                D = Data(region)
                regions_to_nodes[region] = D
                nodes |= {D}
                for linked_region in self.nucmer_data[region].keys():
                    if linked_region not in regions_to_nodes:
                        D = Data(linked_region)
                        regions_to_nodes[linked_region] = D
                        nodes |= {D}

        for region1 in self.nucmer_data.keys():
            for region2 in self.nucmer_data[region1].keys():
                regions_to_nodes[region1].add_link(regions_to_nodes[region2])

        # Find all the connected components
        number = 1

        # initialise outfile
        #outfile = open(outfile, 'w')
        
        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            names = ",".join(names)
            #outfile.write("Group %i:\t%d\t%s\n" % (number,names_count, names))
            print "Group %i:\t%d\t%s" % (number,names_count, names)
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

class TransferGroupsBlast(object):
    def __init__(self, blast_file, hitdata):
        BD              = TFP.BlastFileParser(blast_file)
        self.HD         = TFP.HitData(hitdata)
        self.blastData  = BD.blast_data_tg 
    
    def wrapper(self, group_number, evalue_cutoff, remove_ecoli, outfile):
        self.getConnectedNodes(group_number, evalue_cutoff, remove_ecoli, outfile)
        
    def getConnectedNodes(self,group_number, evalue_cutoff, remove_ecoli, outfile):
        """Connect dem nodes!"""

        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_2_nodes= {}

        for region in self.blastData.keys():
            if self.removeEcoli(remove_ecoli, region): 
                if region not in regions_2_nodes:
                    D = Data(region)
                    regions_2_nodes[region] = D
                    nodes |= {D}
                for linked_region in self.blastData[region].keys():
                    if self.removeEcoli(remove_ecoli, linked_region):
                        if linked_region not in regions_2_nodes:
                            D = Data(linked_region)
                            regions_2_nodes[linked_region] = D
                            nodes |= {D}

        for region1 in self.blastData.keys():
            if self.removeEcoli(remove_ecoli, region1):
                for region2 in self.blastData[region1].keys():
                    if self.removeEcoli(remove_ecoli, region2):
                        evalue = float(self.blastData[region][linked_region][0])
                        if evalue < evalue_cutoff:
                            regions_2_nodes[region1].add_link(regions_2_nodes[region2])

        # Find all the connected components
        number = 1
        # drop to interpretor!
        #code.interact(local=locals())

        # initialise outfile
        outfile = open(outfile, 'w')
        
        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            
            if group_number >0:
                if names_count == group_number:
                    for name in names:
                        print name
                    break
                else:
                    number += 1
            else:
                names = ",".join(names)
                outfile.write("Group %i:\t%d\t%s\n" % (number,names_count, names))
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

    def removeEcoli(self, remove_ecoli, pidsqid):
        if remove_ecoli:
            if self.HD.genus[pidsqid][0] == 'Escherichia' or self.HD.genus[pidsqid][1] == 'Escherichia':
                return False
            else:
                return True
        else:
            return True
    
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
    if args.blast_file:
        TGB = TransferGroupsBlast(args.blast_file,
                                  args.hitdata)
        TGB.wrapper(args.group,
                    args.evalue_cutoff,
                    args.remove_ecoli,
                    args.outfile)
    elif args.nucmer_file:
        TGN = TransferGroupsNucmer(args.hitdata
                                   )
        TGN.wrapper(args.nucmer_file,
                    args.identity,
                    args.length_overlap,
                    args.contigs_file)
    else:
        print 'Please provide either a blast or nucmer file.'

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata', help="")
    parser.add_argument('-b','--blast_file', default=False, help="")
    parser.add_argument('-n','--nucmer_file', default=False, help="")
    parser.add_argument('-c','--contigs_file', help="")
    parser.add_argument('-g','--group', type=int,default=0,help="Print out group-specific files. Leave as default to print all.")
    parser.add_argument('-o','--outfile', default='transfer_groups.csv',help="Output file. Default=transfer_groups.csv")
    parser.add_argument('-evc','--evalue_cutoff', type=int,default=1,help="Set Evalue score cutoff")
    parser.add_argument('-re','--remove_ecoli', default=True,help="Remove E.coli True (default) or False.")
    parser.add_argument('-i','--identity', type=int, default=99,help="Set the identity percentage threshold. e.g. 99 = 99%. Default=99%.")
    parser.add_argument('-lo','--length_overlap', type=int, default=50,help="Set the minimum length overlap percentage. e.g. 99 = 99%. Default=50%.")
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
