#!/usr/bin/env python
###############################################################################
#
# __run_blastn_multiplexed__.py - description!
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
import os
import errno
import glob
import shlex
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Blast(object):
    def __init__(self, path_to_genomes, hitdata):
        self.PATH           = TFP.PathsFileData(path_to_genomes)
        self.HD             = TFP.HitData(hitdata)
        self.contigs        = {}
        self.fasta_files    = {}
        self.blast_cmds     = []
        
    def wrapper(self, gid_contig_file, blastdir, num_threads):
        # parse gid contig file
        self.parseGidContigFile(gid_contig_file)
        
        # grab contigs.fna
        self.generateContigFastas(blastdir)
        
        # read in fasta files
        self.getContigFastaFiles(blastdir)
        
        # run blastn
        self.runBlastn(blastdir, num_threads)
        
    def generateContigFastas(self, blastdir):
        # loop through data, grabbing fasta info from files
        for gid in self.contigs.keys():
            path_to_genome = self.PATH.gid_to_file[gid]
            
            # loop through genome using BioPython
            for accession,sequence in SeqIO.to_dict(SeqIO.parse(path_to_genome,"fasta")).items():
                self.printContigToFile(gid, accession, sequence, blastdir)
                    
    def doesDirectoryExist(self, dir):
        if os.path.exists(dir):
            pass
        else:
            os.mkdir(dir)
    
    def doesFastaFileExist(self, outfile):
        if os.path.exists(outfile):
            return False
        else:
            return True
    
    def printContigToFile(self, gid, accession, seq, blastdir):
        # check that the outdirectory exists
        outdir = os.path.join(blastdir, gid)
        self.doesDirectoryExist(outdir)
        
        # check if the fasta seq exists
        outfile = os.path.join(blastdir, "%s/%s_%s.fna" % (gid, gid, accession))
        if self.doesFastaFileExist(outfile):
            try:
                contig  = self.contigs[gid][accession]
                outfile = os.path.join(blastdir, "%s/%s_%s.fna" % (gid, gid, accession))
                f       = open(outfile, 'w')
                f.write(">%s_%s\n" % (gid,accession))
                f.write("%s\n" % (seq.seq))
                
            except KeyError:
                pass
        
    def parseGidContigFile(self, gid_contig_file):
        with open(gid_contig_file) as fh:
            for l in fh:
                tabs    = l.rstrip().split("\t")
                gid     = tabs[0]
                contig  = tabs[1]
                try:
                    self.contigs[gid][contig] = contig
                except KeyError:
                    self.contigs[gid] = {contig:contig}
    
    def getContigFastaFiles(self, blastdir):
        blastdir = os.path.join(blastdir, "*/*.fna")
        fasta_files = glob.glob(blastdir)
        
        # loop through fasta files
        for fasta_file in fasta_files:
            gid     =  fasta_file.split("/")[-1].split("_")[0]
            contig  = fasta_file.split("/")[-1].split("_")[1:]
            contig  = "_".join(contig).split(".")[0]
            try:
                self.fasta_files[gid][contig] = fasta_file
            except KeyError:
                self.fasta_files[gid] = {contig:fasta_file}
            
    def runBlastn(self, blastdir, num_threads):
        count = 0 
        # loop through fasta files
        for gid in self.fasta_files.keys():
            for contig in self.fasta_files[gid]:
                # check if in gid_contig_file
                if self.inGidContigFile(gid, contig):
                    # check if finished genome
                    if self.checkGenomeStatus(gid):
                        # check if blast output exists
                        if self.doesBLASTOutputExist(blastdir, gid, contig):
                            count +=1 
                            fasta_file      = self.fasta_files[gid][contig]
                            outfile         = os.path.join(blastdir, gid, "%s_%s_vs_nr.blast20151802" % (gid,contig))
                            self.blast_cmds.append("blastn -db /srv/db/ncbi/20151802/nt/nt -query %s -outfmt 1 -out %s" % (fasta_file,
                                                                                                                   outfile))
        # run cmds in multiple threads
        #print str(len(self.blast_cmds))
        
        pool = Pool(num_threads)
        print pool.map(runCommand, self.blast_cmds)
    
    def inGidContigFile(self, gid, contig):
        try:
            lenny = self.contigs[gid][contig]
            return True
        except KeyError:
            return False
    
    def checkGenomeStatus(self, gid):
        try:
            status = self.HD.gid_status[gid]
            # remove Finished and pseudo Finished genomes
            if 'finished' in status.lower():
                return False
            else:
                return True
        except KeyError:
            return False
    
    def doesBLASTOutputExist(self, blastdir, gid, contig):
        outfile = os.path.join(blastdir, gid, "%s_%s_vs_nr.blast20151802" % (gid,contig))
        if os.path.exists(outfile):
            return False
        else:
            return True
               
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

def doWork( args ):
    """ Main wrapper"""
    B = Blast(args.paths_file,
              args.hitdata)
    B.wrapper(args.gid_contig_file,
              args.blastdir,
              args.num_threads)
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('gid_contig_file', help="")
    parser.add_argument('paths_file', help="")
    parser.add_argument('blastdir', help="")
    parser.add_argument('hitdata', help="")
    parser.add_argument('-t','--num_threads', type=int, default=1, help="")
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
