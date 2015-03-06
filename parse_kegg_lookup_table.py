#!/usr/bin/env python
###############################################################################
#
# __parse_kegg_lookup_table__.py - description!
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

class KEGGData(object):
    def __init__(self, kegg_data_file, kegg_mapping_file, kegg_lookup_file):
        self.KMD            = TFP.KEGGMappingData(kegg_mapping_file)
        self.KSTD           = TFP.KEGGSummaryTableData(kegg_data_file)
        self.KLT            = TFP.KEGGLookupTable(kegg_lookup_file)
        self.active_keggs   = {}
        
    def wrapper(self, outfile):
        
        f = open(outfile, 'w')
        
        # print header
        f.write("\t".join(["kegg_id",
                           "KO",
                           "associated_kegg_pathways\n"]))
        
        for kegg_id in self.KSTD.kegg_list.keys():
            try:
                KO = self.KLT.KEGG_lookup[kegg_id]
            except KeyError:
                KO = 'NA'
            pathways = ''
            
            try:
                for pathway in self.KMD.kegg_mapping[KO]:
                    pathways += "%s;" % pathway
                f.write("%s\t%s\t%s\n" % (kegg_id,
                                          KO,
                                          pathways))
            except KeyError:
                # not present in a pathway! 
                f.write("%s\t%s\t%s\n" % (kegg_id,
                                          KO,
                                          'NA'))
        
        

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
    KD = KEGGData(args.kegg_data_file,
                  args.kegg_mapping_file,
                  args.kegg_lookup_file)
    KD.wrapper(args.outfile)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kegg_data_file', help="")
    parser.add_argument('kegg_mapping_file', help="")
    parser.add_argument('kegg_lookup_file', help="")
    parser.add_argument('outfile', help="")
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
