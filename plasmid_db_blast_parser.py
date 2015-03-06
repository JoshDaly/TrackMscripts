#!/usr/bin/env python
###############################################################################
#
# __plasmid_db_blast_parser__.py - Description!
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
import glob
from multiprocessing import Pool
from subprocess import Popen, PIPE
#from Bio import SeqIO
#from Bio.Seq import Seq
#import os
#import errno
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BlastData(object):
    def __init__(self, group_file, taxon_file, hit_file):
        self.GD                         = TFP.GroupData(group_file)
        self.HD                         = TFP.HitData(hit_file, taxon_file)
        self.transfer_groups_plasmid    = {} 
        self.blast_data                 = {}
        self.plasmid_pidsqids           = {}
    
    def parseBLASTFile(self, blastfile):
        with open(blastfile) as fh:
            for l in fh:
                if l[0:6] == "Query=":
                    pidsqid = l.rstrip().split(" ")[-1]

                    for l in fh:
                        if l[0:4] == 'lcl|':
                            whitespace      = l.rstrip().split()
                            straightline    = whitespace[0].rstrip().split('|')
                            plasmid_id      = straightline[1]
                            score = float(whitespace[-2])
                            evalue = float(whitespace[-1])
                            description = " ".join(whitespace[1:-2])
                            
                            # add to pidsqid list
                            self.plasmid_pidsqids[pidsqid] = 1 
                            
                            # collect top hit only
                            self.blast_data[pidsqid] = [plasmid_id, score, evalue, description]
                            
                            # add pidsqid to TG
                            transfer_group = self.GD.group_membership[pidsqid]
                            try:
                                self.transfer_groups_plasmid[transfer_group] += [pidsqid]
                            except KeyError:
                                self.transfer_groups_plasmid[transfer_group] = [pidsqid]
                            
                            break
                        elif l[0:5] == '*****':
                            # no hits! 
                            break
                            
        
    def plasmidHits(self, blast_score_cutoff, evalue_cutoff):
        # loop through pidsqid that hit plasmid DB
        for pidsqid in self.blast_data.keys():
            pidsqid_data    = self.blast_data[pidsqid]
            plasmid_id      = pidsqid_data[0]
            score           = pidsqid_data[1]
            evalue          = pidsqid_data[2]
            description     = pidsqid_data[3]
            if evalue <= evalue_cutoff:
                print '%s %s %f %f %s' % (pidsqid, plasmid_id, score, evalue, description)
            else: 
                pass
    
    def plasmidTransferGroupHits(self, blast_score_cutoff, evalue_cutoff):
         
        for TG in self.GD.group_data.keys():
            # initial variables
            TG_total    = 0
            TG_plasmid  = 0
            
            # determine variables
            TG_total   = len(self.GD.group_data[TG])
            try: 
                TG_plasmid = len(self.transfer_groups_plasmid[TG])
            except KeyError:
                pass
            
            # print
            print "Transfer group: %s\t%d\t%d" % (TG,
                                                  TG_total,
                                                  TG_plasmid)
            
    def getInterVsIntra(self, blast_score_cutoff, evalue_cutoff, type):
        
        pidsqid_total_inter   = 0
        pidsqid_total_intra   = 0
        pidsqid_plasmid_inter = 0
        pidsqid_plasmid_intra = 0
        
        # pidsqids
        for TG in self.GD.group_data.keys():
            for pidsqid in self.GD.group_data[TG]:
                gid1 = self.HD.pidsqid_to_gid[pidsqid][0]
                gid2 = self.HD.pidsqid_to_gid[pidsqid][1]
                phylum1 = self.HD.TD.taxon_phylum[gid1]
                phylum2 = self.HD.TD.taxon_phylum[gid2]
                
                try:
                    pidsqid_data    = self.blast_data[pidsqid]
                    plasmid_id      = pidsqid_data[0]
                    score           = pidsqid_data[1]
                    evalue          = pidsqid_data[2]
                    
                    if phylum1 == phylum2:
                        # intra phyla
                        pidsqid_total_intra += 1 
                        try:
                            lenny = self.plasmid_pidsqids[pidsqid]
                            if evalue <= evalue_cutoff and score >= blast_score_cutoff:
                                pidsqid_plasmid_intra += 1
                        except KeyError:
                            pass
                        
                    else:
                        # inter phyla
                        pidsqid_total_inter += 1 
                        try:
                            carl = self.plasmid_pidsqids[pidsqid]
                            if evalue <= evalue_cutoff and score >= blast_score_cutoff:
                                pidsqid_plasmid_inter += 1
                        except KeyError:
                            pass
                except KeyError:
                    if phylum1 == phylum2:
                        # intra phyla
                        pidsqid_total_intra += 1 
                        try:
                            lenny = self.plasmid_pidsqids[pidsqid]
                            if evalue <= evalue_cutoff and score >= blast_score_cutoff:
                                pidsqid_plasmid_intra += 1
                        except KeyError:
                            pass
                        
                    else:
                        # inter phyla
                        pidsqid_total_inter += 1 
                        try:
                            carl = self.plasmid_pidsqids[pidsqid]
                            if evalue <= evalue_cutoff and score >= blast_score_cutoff:
                                pidsqid_plasmid_inter += 1
                        except KeyError:
                            pass
        
        if type == "intra":            
            print "%f\t%d\t%d\t%f" % (blast_score_cutoff,
                                      pidsqid_total_intra,
                                      pidsqid_plasmid_intra,
                                      pidsqid_plasmid_intra/float(pidsqid_total_intra))
            #print 'intra phyla total: %d plasmid: %d percentage: %f' % (pidsqid_total_intra,
            #                                                            pidsqid_plasmid_intra,
            #                                                            pidsqid_plasmid_intra/float(pidsqid_total_intra))
        elif type == "inter":
            print "%f\t%d\t%d\t%f" % (blast_score_cutoff,
                                      pidsqid_total_inter,
                                      pidsqid_plasmid_inter,
                                      pidsqid_plasmid_inter/float(pidsqid_total_inter))
            #print "inter phyla total: %d plasmid: %d percentage: %f" % (pidsqid_total_inter,
            #                                                            pidsqid_plasmid_inter,
            #                                                            pidsqid_plasmid_inter/float(pidsqid_total_inter))
    
    def BlastDataWrapper(self, blastfile, blast_score_cutoff, evalue_cutoff, type):
        # parse blast data, collected top hits
        self.parseBLASTFile(blastfile)
        
        # print out 
        #self.plasmidHits(blast_score_cutoff,
        #                 evalue_cutoff)      
        
        #self.plasmidTransferGroupHits(blast_score_cutoff,
        
        #                              evalue_cutoff)  
        #for blast_score in range(0, 10100, 100):
        self.getInterVsIntra(blast_score_cutoff,
                             evalue_cutoff,
                             type)
        
        #range = np.arange(0,0.1,0.000001)
        #for evalue in range:
        #self.getInterVsIntra(blast_score_cutoff,
        #                     evalue_cutoff,
        #                     type)

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
    BD = BlastData(args.group_file,
                   args.taxon_file,
                   args.hit_file)
    BD.BlastDataWrapper(args.blast_file,
                        args.blast_score_cutoff,
                        args.evalue_cutoff,
                        args.type)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file', help="Blast output file in outfmt 4.")
    parser.add_argument('group_file', help="File containing transfer groups.")
    parser.add_argument('taxon_file', help="File containing Phil's updated taxonomy.")
    parser.add_argument('hit_file', help="File containing kmer data.")
    parser.add_argument('-bsc','--blast_score_cutoff', type=int, default = 0, help="Set BLAST score cutoff.")
    parser.add_argument('-esc','--evalue_cutoff', type=float, default = 0, help="Set evalue cutoff.")
    parser.add_argument('--type', default = 'inter', help="Set inter or intra.")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    #parser.add_argument('-X', '--optional_X', action="store_true", type=int,default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
