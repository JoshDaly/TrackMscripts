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
#import numpy as np
#np.seterr(all='raise')

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HitDataStats(object):
    def __init__(self, hitdata, transfer_groups):
        self.HD         = TFP.HitData(hitdata)
        self.TG         = TFP.GroupData(transfer_groups)
        
    def wrapper(self, type):
        if type == 'Phylum_interactions':
            self.createPhylumInteractionMatrix()
        elif type == 'genus_interactions':
            pass
        elif type == 'hitdata_fix':
            self.fixTransferGroupsHitData()
    
    def createPhylumInteractionMatrix(self):
        for pidsqid in self.HD.hit_data.keys():
            pass
            
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
                       args.transfer_groups_file)
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
