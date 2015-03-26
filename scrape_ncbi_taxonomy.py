#!/usr/bin/env python
###############################################################################
#
# __scrape_ncbi_taxonomy__.py - description!
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
import urllib
import re
import os
import errno
import shlex
import subprocess

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NCBITaxonomy(object):
    def __init__(self):
        pass
    
    def wrapper(self):
        self.searchNCBITaxonomy()
    
    def searchNCBITaxonomy(self):
        search  = "http://www.ncbi.nlm.nih.gov/taxonomy/?term=txid2%5BSubtree%5D"
        page    = urllib.urlopen(search).read()
        print page

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
    NCBIT = NCBITaxonomy()
    NCBIT.wrapper()
            

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('input_file1', help="")
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
