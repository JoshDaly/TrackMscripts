#!/usr/bin/env python
###############################################################################
#
# __parser_blast_file_minced.trackm__.py - Parse BLAST output for minced output
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

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BlastData(object):
    def __init__(self, blast_file):
        self.BFP = TFP.BlastFileParser(blast_file)
        
    def writeBlastDataToFile(self, outdir, outfile):
        
        # set output
        output = os.path.join(outdir, outfile)
        f = open(output, 'w')
        
        # print header
        f.write("crispr_spacer_id\tcount\tpidsqid_hits\n")
        
        # loop through query ids
        for crispr_spacer in self.BFP.blast_data.keys():
            # reset hits
            hits = ''
            count = 0
            
            for hit in self.BFP.blast_data[crispr_spacer]:
                count += 1
                hits += "%s," % hit 
                
            # write spacer -> hits
            f.write("%s\t%d\t%s\n" % (crispr_spacer,
                                      count,
                                      hits))
            
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""
    BD = BlastData(args.blast_file)
    BD.writeBlastDataToFile(args.out_dir, args.out_file)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file', help="Input file in blast outfmt 4")
    parser.add_argument('out_dir', help="Output directory")
    parser.add_argument('-o','--out_file', default='blast_data_out.txt' ,help="Output filename")
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
