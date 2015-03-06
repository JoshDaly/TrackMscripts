#!/usr/bin/env python
###############################################################################
#
# __make_plasmid_derep_db__.py - Make plasmid db from dereplicated IMG genomes!
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
import re
import glob
from multiprocessing import Pool
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.Seq import Seq
#import os
#import errno
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class PlasmidDB(object):
    def __init__(self, path_file, taxon_file):
        self.derep_genomes  = {} # list 
        self.plasmid_genome = {}
        self.PD             = TFP.PathsFileData(path_file)
        self.TD             = TFP.TaxonomyData(taxon_file)
    
    def parseMetaData(self, derep_metadata):
        with open(derep_metadata) as fh:
            for l in fh:
                if l[0] != 't':
                    tabs = l.rstrip().split("\t")
                    self.derep_genomes[tabs[0]] = 1
    
    def parseIMGGenomes(self, directory, outfile):
        IMG_genomes = glob.glob('%s/*/*.fna' % directory)
        
        # initialise output file
        if outfile:
            out_file = open(outfile,'w')
        
        #### check that it is a genome, and that it is contained within the dereplicated list, then check if it has a plasmid.
        for fasta_file in IMG_genomes:
            
            # check if genome 
            lenny = re.search('[0-9]+.fna',fasta_file)
            if lenny:
            
                # check if genome in dereplicated DB
                img_genome = fasta_file.rstrip().split('/')[-1][:-4]
                try:
                    carl = self.derep_genomes[img_genome]
                    
                    fasta_sequences = SeqIO.parse(open(fasta_file),"fasta")
                    for fasta in fasta_sequences:
                        if 'plasmid' in fasta.description.lower():
                            self.plasmid_genome[fasta.description] = fasta.seq
                            
                            gid = self.PD.img_to_gid[img_genome]
                    
                            if not self.checkIfUpdatedTaxonomy(gid):
                                print gid
                            
                            if outfile:
                                out_file.write('>%s\n' % fasta.description)
                                out_file.write('%s\n' % fasta.seq)
                            
                except KeyError:
                    # replicated genome
                    pass
    
    def checkIfUpdatedTaxonomy(self, gid):
        try:
            lenny = self.TD.taxon_genus[gid]
            return True
        except KeyError:
            return False
    
    def averagePlasmidGenomeLength(self):
        cumulative_length   = []
        
        # calculate the total length of plasmid genomes
        for plasmid in self.plasmid_genome.keys():
            cumulative_length.append(len(self.plasmid_genome[plasmid]))
        
        cumulative_length_np = np.array(cumulative_length)
        
        print "Average dereplicated plasmid genome length %d" % np.average(cumulative_length_np)
        print "Standard deviation %d" % np.std(cumulative_length_np)
    
    def buildPlasmidDB(self, directory, derep_metadata, outfile):
        # build dict of derep genomes
        self.parseMetaData(derep_metadata)
        
        # parse through directory, only grabbing plasmid 
        # genomes of dereplicated IMG genomes! 
        self.parseIMGGenomes(directory, outfile)
        
        # print to screen
        self.averagePlasmidGenomeLength()

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
    PDB = PlasmidDB(args.path_file,
                    args.taxon_file)
    PDB.buildPlasmidDB(args.genome_directory,
                       args.derep_metadata,
                       args.outfile)
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('genome_directory', help="Directory containing dereplicated IMG genomes.")
    parser.add_argument('derep_metadata', help="Dereplicated IMG metadata.")
    parser.add_argument('path_file', help="File containing paths to img files.")
    parser.add_argument('taxon_file', help="File containing updated taxon information.")
    parser.add_argument('-o','--outfile', default = False, help="Output file.")
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
