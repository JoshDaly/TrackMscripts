#!/usr/bin/env python
###############################################################################
#
# __combine_annotation_file_BLAST_groupings__.py - description!
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

import argparse
import sys

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GroupFileParser(object):
    def __init__(self,l):
        self.readGroupFile(l)
        
    def readGroupFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num = tabs[0].split(':')[0]
        self.pid_sqids = tabs[2].split(',')
        
class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnoFile(l)
        
    def readAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.col1  = tabs[0]
        
        
class AnnotationData(object):
    def __init__(self):
        self.groupData = {}
    
    def getGroupData(self,group_data):
        with open(group_data,'r') as fh:
            for l in fh:
                GFP = GroupFileParser(l)
                GFP.readGroupFile(l)
                for pid_sqid in GFP.pid_sqids:
                    try:
                        self.groupData[pid_sqid] += GFP.group_num
                    except KeyError:
                        self.groupData[pid_sqid] = GFP.group_num
    
    def combineData(self,anno_file):
        with open(anno_file,'r') as fh:
            for l in fh: 
                AFP = AnnotationFileParser(l)
                AFP.readAnnoFile(l)
                pid_sqid = "_".join(AFP.col1.split('_')[0:2])
                print "\t".join([self.groupData[pid_sqid],
                                 l.rstrip()])
                
                
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
    # call class
    AD = AnnotationData()
                
    # capture group data
    AD.getGroupData(args.group_file)
    
    # combine with annotation data
    AD.combineData(args.anno_file)


            
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('group_file', help="gut_oral_genomes_total_hits_length_clean.csv")
    parser.add_argument('anno_file', help="gut_oral_genomes_total_hits_length_clean.csv")
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
