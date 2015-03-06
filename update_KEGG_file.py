#!/usr/bin/env python
###############################################################################
#
# __updatE_KEGG_file__.py - description!
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

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KEGGSummmaryTableParser(object):
    def __init__(self, l):
        self.readKEGGSummaryTable(l)
        
    def readKEGGSummaryTable(self, l):
        tabs = l.rstrip().split("\t")
        self.kegg               = tabs[0]
        self.description        = tabs[1]
        self.num_genus          = tabs[2]
        self.num_phyla          = tabs[3]
        self.num_habitats       = tabs[4]
        self.transfers_avg      = tabs[5]
        self.transfer_std       = tabs[6]
        self.contig_avg         = tabs[7]
        self.contig_std         = tabs[8]
        self.transfer_groups    = tabs[9]

class UpdateKEGG(object):
    def __init__(self,kegg_file):
        self.KD                 = TFP.KEGGData(kegg_file)
        self.kegg_description   = {}
    
    def wrapper(self, kegg_summary_file):
        # read in KEGG summary table data
        self.getKEGGSummaryData(kegg_summary_file)
        
        # write descriptions to file
        self.printData()
        
    def getKEGGSummaryData(self, kegg_summary_file):
        with open(kegg_summary_file) as fh:
            for l in fh:
                KSTP = KEGGSummmaryTableParser(l)
                self.kegg_description[KSTP.kegg] = KSTP.description
       
    def getKEGGDescription(self, kegg_id):
        description = ''
        try:
            description = self.kegg_description[kegg_id]
            return description
        except KeyError:
            return description
                
    def printData(self):
        for pidsqid_gene in self.KD.all_kegg_data.keys():
            kegg_id     = self.KD.all_kegg_data[pidsqid_gene][1]
            description = self.getKEGGDescription(kegg_id)
            print "\t".join([self.KD.all_kegg_data[pidsqid_gene][0],
                             self.KD.all_kegg_data[pidsqid_gene][1],
                             description,
                             self.KD.all_kegg_data[pidsqid_gene][2],
                             self.KD.all_kegg_data[pidsqid_gene][3],
                             self.KD.all_kegg_data[pidsqid_gene][4],
                             self.KD.all_kegg_data[pidsqid_gene][5],
                             self.KD.all_kegg_data[pidsqid_gene][6],
                             self.KD.all_kegg_data[pidsqid_gene][7],
                             self.KD.all_kegg_data[pidsqid_gene][8],
                             self.KD.all_kegg_data[pidsqid_gene][9],
                             self.KD.all_kegg_data[pidsqid_gene][10],
                             self.KD.all_kegg_data[pidsqid_gene][11]
                             ])

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
    UK = UpdateKEGG(args.kegg_file)
    UK.wrapper(args.kegg_summary_file)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kegg_file', help="")
    parser.add_argument('kegg_summary_file', help="")
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
