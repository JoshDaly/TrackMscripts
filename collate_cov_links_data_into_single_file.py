#!/usr/bin/env python
###############################################################################
#
# __collate_cov_links_data_into_single_file__.py - description!
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
import glob
#import os
#import errno
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq
#import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CollateFiles(object):
    def __init__(self):
        self.bamm_cov_links_files   = {}
        self.bamm_data              = {}
    
    def wrapper(self, draft_cov_link_dir, finished_cov_link_dir): #links_outfile, cov_outfile):
        # grab bamm stats files
        self.getCovLinkFiles(draft_cov_link_dir)
        self.getCovLinkFiles(finished_cov_link_dir)
        
        # grab bamm cov links data
        self.getCovLinkData()
        
        # print out coverage data
        f = open('cov_data_2.csv','w')
        self.printCovData(f)
        
    def getCovLinkFiles(self, cov_link_dir):
        if cov_link_dir[-1] == '/':
            # remove last element
            cov_link_dir = cov_link_dir[0:-1]
            
        cov_link_files = glob.glob("%s/*/*cov_links_stats.csv" % cov_link_dir)
        for file in cov_link_files:
            gid = file.split("/")[-2]
            self.bamm_cov_links_files[gid] = file
    
    def getCovLinkData(self):
        for gid in self.bamm_cov_links_files.keys():
            cov_link_file = self.bamm_cov_links_files[gid]
            CLD = TFP.CovLinksData(cov_link_file)
            for contig in CLD.cov_link_data[gid]:
                cov_data  = CLD.cov_link_data[gid][contig][1]
                contig_len= CLD.cov_link_data[gid][contig][0]
                try:
                    link_data = CLD.cov_link_data[gid][contig][3]
                    try:
                        self.bamm_data[gid].append("%s:%s:%s:%s"% (contig,cov_data, link_data, contig_len))
                    except KeyError:
                        self.bamm_data[gid] = ["%s:%s:%s:%s"% (contig,cov_data, link_data, contig_len)]
                except (IndexError, KeyError):
                    link_data = '0'
                    try:
                        self.bamm_data[gid].append("%s:%s:%s:%s"% (contig,cov_data, link_data, contig_len))
                    except KeyError:
                        self.bamm_data[gid] = ["%s:%s:%s:%s"% (contig,cov_data, link_data, contig_len)]
    
    def printCovData(self, outfile):
        for gid in self.bamm_data.keys():
            string_to_print = ''
            for i in self.bamm_data[gid]:
                string_to_print += "\t%s" % i
            outfile.write("%s%s\n" % (gid,string_to_print)) 

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
    CF = CollateFiles()
    CF.wrapper(args.draft_cov_link_dir, args.finished_cov_link_dir) #args.links_outfile, args.cov_outfile)
                
            
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('draft_cov_link_dir', help="")
    parser.add_argument('finished_cov_link_dir', help="")
    #parser.add_argument('links_outfile', help="")
    #parser.add_argument('cov_outfile', help="")
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
