#!/usr/bin/env python
###############################################################################
#
# __hit_data_stats__.py - description!
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
import math as math
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

class HitData(object):
    def __init__(self, hitdata_file, taxon_file):
        self.HD                     = TFP.HitData(hitdata_file, taxon_file)
        self.contig_size_stats      = {}
        self.genome_status_stats    = {} 
        self.total_inter_pidsqid    = 0
    
    def wrapper(self, stat_type):
        try:
            if stat_type.lower() == 'contig_size':
                self.getContigSizeStats()
            elif stat_type.lower() == 'genome_status_of_pidsqids':
                self.getGenomeStatusStats()
        except  AttributeError:
            print "***ERROR***"
            print ">>>>>>Please select either -st --stat_type as contig_size or genome_status_of_pidsqids"
            os.sys.exit()
            
    def getContigSizeStats(self):
        for pidsqid in self.HD.hit_data.keys():
            if self.HD.intra_or_inter[pidsqid] == 'inter':
                
                # add 1 to inter phyla pidsqid count 
                self.total_inter_pidsqid += 1
                
                # A = <700,000
                # B = <700,000 & >700,000
                # C = >700,000 & >700,000
                contig_size1 = self.HD.contig_size_gid[pidsqid][0]
                contig_size2 = self.HD.contig_size_gid[pidsqid][1]
                contig_size_cutoff = 700000
                
                # gather contig length data
                self.checkContigSize(contig_size1, contig_size2, contig_size_cutoff)
            
        # print data
        self.returnContigSizeStats()
            
    def checkContigSize(self, contig_size1, contig_size2, contig_size_cutoff):
        if contig_size1 <= contig_size_cutoff and contig_size2 <= contig_size_cutoff:
            self.addCountToDict(self.contig_size_stats, 'A')
        elif contig_size1 <= contig_size_cutoff and contig_size2 > contig_size_cutoff:
            self.addCountToDict(self.contig_size_stats, 'B')
        elif contig_size1 > contig_size_cutoff and contig_size2 <= contig_size_cutoff:
            self.addCountToDict(self.contig_size_stats, 'B')
        elif contig_size1 > contig_size_cutoff and contig_size2 > contig_size_cutoff:
            self.addCountToDict(self.contig_size_stats, 'C')
            
    def addCountToDict(self, dict, type):
        try:
            dict[type] += 1 
        except KeyError:
            dict[type] = 1 
            
    def returnContigSizeStats(self):
        for type in self.contig_size_stats.keys():
            percentage = str(np.round(self.contig_size_stats[type]/float(self.total_inter_pidsqid), 2))
            count = self.contig_size_stats[type]
            print "%s: %s %d" % (type,
                                 percentage,
                                 count)
        
    def getGenomeStatusStats(self):
        for pidsqid in self.HD.hit_data.keys():
            if self.HD.intra_or_inter[pidsqid] == 'inter':
                
                # add 1 to inter phyla pidsqid count 
                self.total_inter_pidsqid += 1
                
                genome_status1 = self.HD.genome_status[pidsqid][0]
                genome_status2 = self.HD.genome_status[pidsqid][1]

                # gather data
                self.checkGenomeStatus(genome_status1, genome_status2)
        
        # print data
        self.returnGenomeStatusStats()
                
    def checkGenomeStatus(self, genome_status1, genome_status2):
        if "finished" in genome_status1.lower() and 'finished' in genome_status2.lower():
            self.addCountToDict(self.genome_status_stats, 'finished')
        elif "draft" in genome_status1.lower() and 'finished' in genome_status2.lower():
            self.addCountToDict(self.genome_status_stats, 'draft+finished')
        elif "finished" in genome_status1.lower() and 'draft' in genome_status2.lower():
            self.addCountToDict(self.genome_status_stats, 'draft+finished')
        elif "draft" in genome_status1.lower() and 'draft' in genome_status2.lower():
            self.addCountToDict(self.genome_status_stats, 'draft')
            
    def returnGenomeStatusStats(self):
        for type in self.genome_status_stats.keys():
            percentage = str(np.round(self.genome_status_stats[type]/float(self.total_inter_pidsqid), 2))
            count = self.genome_status_stats[type]
            print "%s: %s %d" % (type,
                                 percentage,
                                 self.genome_status_stats[type])
            
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
    HD = HitData(args.hitdata_file,
                 args.taxon_file)
    HD.wrapper(args.stat_type)
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata_file', help="")
    parser.add_argument('taxon_file', help="")
    parser.add_argument('-st','--stat_type', help="Indicate which statistic you would like to print: contig_size or genome_status_of_pidsqids")
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
