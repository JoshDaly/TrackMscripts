#!/usr/bin/env python
###############################################################################
#
# __cog_centric_table__.py - description!
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
import os
import errno
import numpy as np
np.seterr(all='raise')
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

class COGData(object):
    def __init__(self, hitdata_file, taxon_file, group_file, annotation_file):
        self.HD                 = TFP.HitData(hitdata_file, taxon_file)
        self.GD                 = TFP.GroupData(group_file)
        self.AD                 = TFP.AnnotationData(annotation_file)
        self.cogs               = {}
        self.cogs_phylum        = {}
        self.cog_genera         = {}
        self.cog_habitat        = {}
        self.cog_tgs            = {}
        self.cogs_avg_transfer  = {}
        self.cogs_avg_contigs   = {}
        self.transfer_group     = {}
        
    def wrapper(self, outdirectory):
        
        # loop through COGs
        for cog in self.AD.cog_to_pidsqid.keys():
            for pidsqid in self.AD.cog_to_pidsqid[cog]:
                if pidsqid in self.GD.group_membership:
                    
                    # call hitdata with pidsqid
                    phylum1         = self.HD.phylum[pidsqid][0]
                    phylum2         = self.HD.phylum[pidsqid][1]
                    genus1          = self.HD.genus[pidsqid][0]
                    genus2          = self.HD.genus[pidsqid][1]
                    habitat1        = self.HD.habitat[pidsqid][0]
                    habitat2        = self.HD.habitat[pidsqid][1]
                    transferSize    = self.HD.transfer_size[pidsqid]
                    contigSize1     = self.HD.contig_size_gid[pidsqid][0]
                    contigSize2     = self.HD.contig_size_gid[pidsqid][1]
                    transferGroup   = self.GD.group_membership[pidsqid]
                    
                    # check if intra or inter
                    if self.HD.intra_or_inter[pidsqid] == 'inter':
                        # phylum 
                        self.addCOGData(self.cogs_phylum, cog, phylum1, 'inter')
                        self.addCOGData(self.cogs_phylum, cog, phylum2, 'inter')
                        # genus
                        self.addCOGData(self.cog_genera, cog, genus1, 'inter')
                        self.addCOGData(self.cog_genera, cog, genus2, 'inter')
                        # habitat
                        self.addCOGData(self.cog_habitat, cog, habitat1, 'inter')
                        self.addCOGData(self.cog_habitat, cog, habitat2, 'inter')
                        # transfer length
                        self.addCOGLengthData(self.cogs_avg_transfer, cog, transferSize, 'inter')
                        # contig length
                        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize1, 'inter')
                        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize2, 'inter')
                        # transfer group
                        self.addCOGData(self.transfer_group, cog, transferGroup, 'inter')
                    
        self.printData(outdirectory)
        """                    
        else:
        
        # phylum 
        self.addCOGData(self.cogs_phylum, cog, phylum1, 'intra')
        self.addCOGData(self.cogs_phylum, cog, phylum2, 'intra')
        # genus
        self.addCOGData(self.cog_genera, cog, genus1, 'intra')
        self.addCOGData(self.cog_genera, cog, genus2, 'intra')
        # habitat
        self.addCOGData(self.cog_habitat, cog, habitat1, 'intra')
        self.addCOGData(self.cog_habitat, cog, habitat2, 'intra')
        # transfer length
        self.addCOGLengthData(self.cogs_avg_transfer, cog, transferSize, 'intra')
        # contig length
        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize1, 'intra')
        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize2, 'intra')
        # transfer group
        self.addCOGData(self.transfer_group, cog, transferGroup, 'intra')
        """     
        
    
    def printData(self, outdirectory):
        # output files
        outfile_inter = os.path.join(outdirectory, "cog_centric_inter.all_TGs.csv")
        #outfile_intra = os.path.join(outdirectory, "cog_centric_intra.all_TGs.csv")
        
        f1 = open(outfile_inter, 'w')
        #f2 = open(outfile_intra, 'w') 
        
        # write header
        f1.write("cog\tdescription\tnum_genus\tnum_phyla\tnum_habitats\ttransfers_avg\ttransfer_std\tcontig_avg\tcontig_std\ttransfer_groups\n")
        #f2.write("cog\tdescription\tnum_genus\tnum_phyla\tnum_habitats\ttransfers_avg\ttransfer_std\tcontig_avg\tcontig_std\ttransfer_groups\n")
        
        # print inter phylum transfers
        for cog in self.cog_genera['inter']:
            num_genus                   = len(self.cog_genera['inter'][cog])
            num_phyla                   = len(self.cogs_phylum['inter'][cog])
            num_habitats                = len(self.cog_habitat['inter'][cog])
            transfer_len_array          = self.cogs_avg_transfer['inter'][cog]
            transfer_avg,transfer_std   = self.averageLength(transfer_len_array)
            contig_len_array            = self.cogs_avg_contigs['inter'][cog]
            contig_avg,contig_std       = self.averageLength(contig_len_array)
            num_transfer_groups         = len(self.transfer_group['inter'][cog])
            description                 = self.AD.cog_description[cog]
            
            f1.write("%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\n" % (cog,
                                                                   description,
                                                                   num_genus,
                                                                   num_phyla,
                                                                   num_habitats,
                                                                   transfer_avg,
                                                                   transfer_std,
                                                                   contig_avg,
                                                                   contig_std,
                                                                   num_transfer_groups
                                                                   ))
        """            
        for cog in self.cog_genera['intra']:
        num_genus                   = len(self.cog_genera['intra'][cog])
        num_phyla                   = len(self.cogs_phylum['intra'][cog])
        num_habitats                = len(self.cog_habitat['intra'][cog])
        transfer_len_array          = self.cogs_avg_transfer['intra'][cog]
        transfer_avg,transfer_std   = self.averageLength(transfer_len_array)
        contig_len_array            = self.cogs_avg_contigs['intra'][cog]
        contig_avg,contig_std       = self.averageLength(contig_len_array)
        num_transfer_groups         = len(self.transfer_group['intra'][cog])
        description                 = self.AD.cog_description[cog]
        
        f2.write("%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\n" % (cog,
        description,
        num_genus,
        num_phyla,
        num_habitats,
        transfer_avg,
        transfer_std,
        contig_avg,
        contig_std,
        num_transfer_groups
        ))
        """
    
    def averageLength(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)
    
    def addCOGLengthData(self, dict, cog, key, type):
        try:
            dict[type][cog].append(key)
        except KeyError:
            try:
                dict[type][cog] = [key]
            except KeyError:
                dict[type] = {cog:[key]}
                
    def addCOGData(self, dict, cog, key, type):
        try:
            dict[type][cog][key] = 1 
        except KeyError:
            try:
                dict[type][cog] = {key:1}
            except KeyError:
                dict[type] = {cog:{key:1}}
    

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
    CD = COGData(args.hitdata_file,
                 args.taxon_file,
                 args.group_file,
                 args.annotation_file)

    CD.wrapper(args.outdirectory)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata_file', help="")
    parser.add_argument('taxon_file', help="")
    parser.add_argument('group_file', help="")
    parser.add_argument('annotation_file', help="")
    parser.add_argument('outdirectory', help="")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
