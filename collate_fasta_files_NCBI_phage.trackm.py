#!/usr/bin/env python
###############################################################################
#
# __collate_fasta_file_NCBI_phage.trackm__.py - description!
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
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os
import errno
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NCBIPhageDB(object):
    def __init__(self):
        self.genome_description = {}
    
    def wrapper(self, directory, outdirectory):
        
        # list of viral genomes
        viral_genomes = glob.glob('%s/*/*.fna' % directory)
        
        # initialise output files
        outfile_fna = os.path.join(outdirectory,"ncbi_phage_db.fna")
        outfile_txt = os.path.join(outdirectory,"ncbi_phage_db.txt")
        f1 = open(outfile_fna,'w')
        f2 = open(outfile_txt, 'w')
        
        # write header to file
        f2.write("ncbi_id\tgenome_name\tlength\n")
        
        for fasta_file in viral_genomes:
            
            # get genome description from path
            genome_description  = fasta_file.split("/")[-2]
            ncbi_id             = fasta_file.split("/")[-1].split(".")[0]
            
            # add to dict
            self.genome_description[ncbi_id] = genome_description
            
            if 'phage' in genome_description.lower():
                fasta_sequences = SeqIO.parse(open(fasta_file),"fasta")
                
                # just in case it is in parts, though should be one piece
                for fasta in fasta_sequences:
                    f1.write(">%s\n" % fasta.description)
                    f1.write("%s\n" % fasta.seq)
                    
                # write to 
                f2.write("%s\t%s\t%d\n" % (ncbi_id,
                                                  genome_description,
                                                  len(fasta.seq)))

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

    PHAGE = NCBIPhageDB()
    PHAGE.wrapper(args.genome_directory,
                  args.output_directory)            
            
        

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_directory', help="")
    parser.add_argument('output_directory', help="")
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
