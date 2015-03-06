#!/usr/bin/env python
###############################################################################
#
# __symbolic_links_for_sras__.py - description!
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
import subprocess
import os
import errno
#import numpy as np
#np.seterr(all='raise')
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

class SymbolicLinks(object):
    def __init__(self, paths_file):
        self.Path   = TFP.PathsFileData(paths_file)
        self.gids   = {}
        self.cmds   = []
        
    def wrapper(self, gids_file):
        # open and store gids file
        self.getGids(gids_file)
        
        # make symbolic links 
        self.makeSymbolicLinks()
        
        # run commands
        for cmd in self.cmds:
            runCommand(cmd)
        
    def makeSymbolicLinks(self):
        for gid in self.gids.keys():
            path = self.Path.gid_to_file[gid]
            outdirectory = '/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/transfer_groups/sra/genomes'
            self.cmds.append('ln -s %s %s' % (path,
                                              outdirectory))
        
    def getGids(self, gids_file):
        with open(gids_file) as fh:
            for l in fh:
                gid = l.rstrip()
                self.gids[gid] = 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
        """Run a command and take care of stdout
    
        expects 'cmd' to be a string like "foo -b ar"
    
        returns (stdout, stderr)
        
        Must be outside class object!!!!!!!
        """
        print cmd
        p = subprocess.Popen(cmd, shell=True)
        return p.communicate()  

def doWork( args ):
    """ Main wrapper"""
    SL = SymbolicLinks(args.paths_file)
    SL.wrapper(args.gids_file)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('paths_file', help="")
    parser.add_argument('gids_file', help="")
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
