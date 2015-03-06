#!/usr/bin/env python
###############################################################################
#
# __varun_file_format_KEGG__.py - description!
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

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

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
    # objects
    kegg_lookup = {}
    KOs         = {}

    # read in file, look up table
    with open(args.lookup_file) as fh:
        for l in fh:
            tabs    = l.rstrip().split("\t")
            KO      = tabs[0] 
            kegg_id = tabs[1]
            kegg_lookup[kegg_id]  = KO
    
    # read in inter phyla transfers kegg annotated file
    with open(args.data_file) as fh:
        for l in fh:
            tabs    = l.rstrip().split("\t")
            kegg_id = tabs[1]
            KO = kegg_lookup[kegg_id]
            
            # add to KOs library
            KOs[KO] = 1
    # print out KOs
    for KO in KOs.keys():
        print KO
    
            
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file1', help="gut_oral_genomes_total_hits_length_clean.csv")
    parser.add_argument('input_file1', help="gut_oral_genomes_total_hits_length_clean.csv")
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
