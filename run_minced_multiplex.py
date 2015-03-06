#!/usr/bin/env python
###############################################################################
#
# __run_minced__.py - description!
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
import subprocess
import os
import errno

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Batch7GIDs(object):
    def __init__(self, gids_file):
        self.gids = {}
        self.readGIDsFile(gids_file)
        
    def readGIDsFile(self, gids_file):
        with open(gids_file) as fh:
            for l in fh:
                gid = l.rstrip()
                self.gids[gid] = 1 
    
class CRISPRData(object):
    def __init__(self, paths_file, gids_file, outdir, num_threads):
        self.Paths          = TFP.PathsFileData(paths_file)
        self.Batch7         = Batch7GIDs(gids_file)
        self.out_dir        = outdir
        self.num_threads    = num_threads
        self.minced_cmds    = []
        self.mkdir_cmds     = []
    
    def wrapper(self):
        # read through genome files
        for gid in self.Paths.gid_to_file.keys():
            
            # check gid from batch 7
            if self.checkGid(gid):
                
                # add genome mkdir to cmds
                self.makeDirectory(gid, self.out_dir)
                
                # add minced to cmds
                path_to_genome = self.Paths.gid_to_file[gid]
                self.makeMincedCommand(gid, path_to_genome, self.out_dir)
                
        # run mkdir commands
        for cmd in self.mkdir_cmds:
            runCommand(cmd)
        
        # run minced commands
        pool = Pool(self.num_threads)
        print pool.map(runCommand, self.minced_cmds)
        
    def makeDirectory(self, gid, output_dir):
        """
        mkdir /PATH/genome_spacers/gid/
        """
        output = os.path.join(output_dir, gid)
        self.mkdir_cmds.append("mkdir %s" % (output))
    
    def makeMincedCommand(self, gid, path_to_genome, output_dir):
        """
        minced -spacers fasta_file output_file
        """
        output = os.path.join(output_dir, gid, gid)
        self.minced_cmds.append("minced -spacers %s %s" % (path_to_genome,
                                                           output
                                                           ))
    
    def checkGid(self, gid):
        try:
            lenny = self.Batch7.gids[gid]
            return True
        except KeyError:
            return False 

def runCommand(cmd):
        """Run a command and take care of stdout
    
        expects 'cmd' to be a string like "foo -b ar"
    
        returns (stdout, stderr)
        
        Must be outside class object!!!!!!!
        """
        print cmd
        p = subprocess.Popen(cmd, shell=True)
        return p.communicate()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""

    CD = CRISPRData(args.paths_file,
                    args.gids_file,
                    args.out_dir,
                    args.num_threads)
    CD.wrapper()


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('paths_file', help="File containing paths to genomes")
    parser.add_argument('gids_file', help="File containing batch7 gids")
    parser.add_argument('out_dir', help="Specify output directory")
    parser.add_argument('-t','--num_threads', type=int, default = 10, help="Specify the number of threads")
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
