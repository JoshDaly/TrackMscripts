#!/usr/bin/env python
###############################################################################
#
# __hitdata_stats__.py - description!
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
import os
import errno
import matplotlib.pyplot as plt
#import numpy as np
#np.seterr(all='raise')

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HitDataStats(object):
    def __init__(self, hitdata, transfer_groups, contam_pidsqids):
        self.HD         = TFP.HitData(hitdata)
        self.TG         = TFP.GroupData(transfer_groups)
        if contam_pidsqids:
            self.CP         = TFP.ContaminatedPidsqids(contam_pidsqids) 
        
    def wrapper(self, type):
        if type == 'phylum_interactions':
            self.createPhylumInteractionMatrix()
        elif type == 'genus_interactions':
            pass
        elif type == 'hitdata_fix':
            self.fixTransferGroupsHitData()
    
    def createPhylumInteractionMatrix(self):
        phylum_interactions = {}
        phylums = {}
        matrix = []
        
        for pidsqid in self.HD.hit_data.keys():
            if self.checkPidsqid(pidsqid):
                phylum1 = self.HD.phylum[pidsqid][0]
                phylum2 = self.HD.phylum[pidsqid][1]
                # add phylum
                phylums[phylum1] = 1
                phylums[phylum2] = 1
                # add phylum interaction data
                self.addPhylum(phylum1, phylum2, pidsqid, phylum_interactions)
                self.addPhylum(phylum2, phylum1, pidsqid, phylum_interactions)
        
        phylum_array = phylums.keys()
        
        # add zeros to matrix
        self.zeroArray(matrix, phylum_array)
        
        for i in range(len(phylum_array)):
            hits_to_append = []
            for v in range(i, len(phylum_array)):
                phylum1 = phylum_array[i]
                phylum2 = phylum_array[v]
                try:
                    hits = len(phylum_interactions[phylum1][phylum2])
                except KeyError:
                    hits = 0
                matrix[i][v] = hits
        
        # make scatter plot
        xs = []
        ys = []
        
        for i,v in enumerate(matrix):
            for x in matrix[i]:
                print i,x
                #xs.append(i)
                #ys.append(v)
                
        #set the number of tick labels
        for i in range(len(phylum_array)):
            ticks.append(i)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        
        # set labels as phyla
        labels = [item.get_text() for item in ax.get_xticklabels()]
        # rotate x labels 90 degrees
        plt.xticks(rotation=90)
        
        plt.tight_layout()
        
        plt.scatter(xs, ys)
        
        plt.show()
        
        """
        # print header
        header = " " # empty top left corner
        for i in phylum_array:
            header += "\t%s" % i
        print header
            
        for i in range(len(phylum_array)):
            string_to_print = "%s" % phylum_array[i]
            for v in range(len(phylum_array)):
                string_to_print += "\t%d" % matrix[i][v]
            print string_to_print
        """
        
    def checkPidsqid(self, pidsqid):
        try: 
            if pidsqid in self.CP.contam_pidsqids:
                return False
            else: 
                return True
        except AttributeError:
            return True
    
    def zeroArray(self, array, phylum_array):
        for i in range(len(phylum_array)):
            zeros = []
            for v in range(len(phylum_array)):
                zeros.append(0)
            array.append(zeros)
            
    def addPhylum(self, phylum1, phylum2, pidsqid, dict):
        try:
            dict[phylum1][phylum2][pidsqid] =1 
        except KeyError:
            try:
                dict[phylum1][phylum2] = {pidsqid:1} 
            except KeyError:
                dict[phylum1] = {phylum2:{pidsqid:1}}
            
    def fixTransferGroupsHitData(self):
        # print header 
        print '\t'.join(['pidsqid','transfer_group','hid','pid','ani_comp','ident','gid_1','habitat_1','phylum_1','genus_1','status_1','sequencingMethod_1','sequencingCentre_1','horizontalTransferred_1','genomeSize_1','scaffoldCount_1','len_1','strand_1','cid_1','contig_1','contigLength_1','sqid_1','start_1','gid_2','habitat_2','phylum_2','genus_2','status_2','sequencingMethod_2','sequencingCentre_2','horizontalTransferred_2','genomeSize_2','scaffoldCount_2','len_2','strand_2','cid_2','contig_2','contigLength_2','sqid_2','start_2'])
        
        for pidsqid in self.HD.hit_data.keys():
            string_to_print = ''
            string_to_print += '%s' % pidsqid
            try:
                transfer_group = self.TG.group_membership[pidsqid]
            except KeyError:
                transfer_group = 'NA'
            for i,v in enumerate(self.HD.hit_data[pidsqid]):
                if i == 0:
                    string_to_print += '\t%s' % transfer_group
                else:
                    string_to_print += '\t%s' % v
            print string_to_print
                
        
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
    HDS = HitDataStats(args.hitdata,
                       args.transfer_groups_file,
                       args.contaminted_pidsqids)
    HDS.wrapper(args.type)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata', help="")
    parser.add_argument('-tg','--transfer_groups_file', help="")
    parser.add_argument('-t','--type', help="Type of summary table to create. Phylum_interactions, genus_interactions.")
    parser.add_argument('-cp','--contaminted_pidsqids', default=False, help="Remove contaminated pidsqids file. Default=False")
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
