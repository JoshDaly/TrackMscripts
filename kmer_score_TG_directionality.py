#!/usr/bin/env python
###############################################################################
#
# __kmer_score_TG_directionality__.py - Script to calculate the kmer score or directionality of each Transfer group member. 
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

import argparse
import sys
import glob

from multiprocessing import Pool
from subprocess import Popen, PIPE

from Bio import SeqIO
from Bio.Seq import Seq

from scipy.spatial.distance import pdist

#import os
#import errno

import numpy as np
np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class PathsFileParser(object):
    def __init__(self,l):
        self.readPathsFile(l)
    
    def readPathsFile(self,l):
        tabs = l.rstrip().split("\t")
        self.gid             = tabs[0]
        self.path_to_file    = tabs[1]
        self.img             = tabs[1].rstrip().split("/")[-2]
        
class PathsFileData(object):
    def __init__(self, paths_file):
        self.gid_to_file = {}
        self.img_to_gid  = {}
        self.buildPathsData(paths_file)
    
    def buildPathsData(self, paths_file):
        with open(paths_file) as fh:
            for l in fh:
                if l[0] != "#":
                    PFP = PathsFileParser(l)
                    PFP.readPathsFile(l)
                    self.gid_to_file[PFP.gid] = PFP.path_to_file
                    self.img_to_gid[PFP.img] = PFP.gid

class LGTInfoStore( object ):
    """Helper class for storing information about directionality of LGT events"""
    # __init__ links to object, and must be present in every Class! 
    def __init__( self ):
        self.lgtTmer = {}
        self.genomeTmers = {}
        self.lgtGenomes = {}
        self.what = None
        self.Dist_dict = {} # store distances
        self.lgtScores = {}
            
    def addLGT( self, pidsqid, gid1, gid2 ):
        """populate lgtgenomes dictionary"""
        self.lgtScores[pidsqid] = None
        self.lgtGenomes[pidsqid] = [gid1, gid2] 

    def addLGTTmer( self, pidsqid, tmer ):
        self.lgtTmer[pidsqid] = tmer 
        
    def addGenomeTmer( self, gid, tmer ):
        self.genomeTmers[gid] = tmer

    def getClosestGID( self , rounded):
        """Calculate the kmer score
        
        returns (score, (closestGID, dist), (furthestGID, dist))
        """
        LGTs = self.lgtGenomes.keys()
        dgs = []
        for pidsqid in LGTs:
            
            dg1 = pdist([self.genomeTmers[self.lgtGenomes[pidsqid][0]], self.lgtTmer[pidsqid] ])
            dg2 = pdist([self.genomeTmers[self.lgtGenomes[pidsqid][1]], self.lgtTmer[pidsqid] ])
            dg1_str = ''.join(map(str,dg1))
            dg2_str = ''.join(map(str,dg2)) 
            rounded_score = float(np.round(dg1/(dg1+dg2),decimals=2))
            score = float(dg1/(dg1+dg2))
            
            self.lgtScores[pidsqid] = [score,float(np.mean([dg1,dg2])),dg1_str,dg2_str]
            
            if rounded:
                try:
                    self.Dist_dict[rounded_score]+=1
                except KeyError:
                    self.Dist_dict[rounded_score]=1
            else:
                self.Dist_dict[score]=[float(np.mean([dg1,dg2])),dg1_str,dg2_str]
                
    def printInfoHeader(self, out_file):
        out_file.write("pidsqid\tgid_1\tgid_2\tkmer_score\tmean_dg\tdg1\tdg2\n")
    
    def printInfo(self, out_file):
        for pidsqid in self.lgtScores.keys():
            kmer_score  = str(self.lgtScores[pidsqid][0])
            mean_dg     = str(self.lgtScores[pidsqid][1])
            dg1         = str(self.lgtScores[pidsqid][2])
            dg2         = str(self.lgtScores[pidsqid][3])
            g1          = self.lgtGenomes[pidsqid][0]
            g2          = self.lgtGenomes[pidsqid][1]
            out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (pidsqid,
                                                         g1,
                                                         g2,
                                                         kmer_score,
                                                         mean_dg,
                                                         dg1,
                                                         dg2))
        
    def getDistHisto(self):
        for score in self.Dist_dict.keys():
            print "\t".join([str(score),str(self.Dist_dict[score][0]),str(self.Dist_dict[score][1]),str(self.Dist_dict[score][2])])
        
    def printDict( self ):
        for uid in self.lgtGenomes:
            print "\t".join([uid,self.lgtGenomes[uid][0],self.lgtGenomes[uid][1]])

    def __str__( self ):
        """print function"""
        return "GID1: %s GID2: %s" % (i for i in self.genomeTmers.keys())

class KmerData(object):
    def __init__(self, paths_file):
        self.LIS            = LGTInfoStore()
        self.lgt_genomes    = {}  
        self.PFD            = PathsFileData(paths_file)
    
    def kmerDataWrapper(self, kmer_dirs, outfile):
        
        # capture kmer directories
        kmer_directories = glob.glob('%s/*/*' % kmer_dirs)
        
        # loop through kmer directories
        for kmer_dir in kmer_directories:
            
            # capture kmer files
            kmer_files = glob.glob('%s/*.kmer*' % kmer_dir)
            
            pidsqid = ''
            genome1 = ''
            gemome2 = ''
            
            # loop through kmer files
            for kmer_file in kmer_files:
                
                # check if genome or lgt file
                lenny = self.checkKmerFile(kmer_file)
                
                # genome
                if lenny: 
                    if len(genome1) > 0:
                        genome2 = self.PFD.img_to_gid[kmer_file.rstrip().split("/")[-1].split('.')[0].split("_")[-1]]
                        self.grabGenomeKmerInfo(kmer_file,
                                                genome2)
                    else:
                        genome1 = self.PFD.img_to_gid[kmer_file.rstrip().split("/")[-1].split('.')[0].split("_")[-1]]
                        self.grabGenomeKmerInfo(kmer_file,
                                                genome1)

                # lgt file
                else:
                    pidsqid = kmer_file.rstrip().split("/")[-2]
                    self.grabLGTKmerInfo(kmer_file,
                                         pidsqid)
                    
            # add to dictionary
            self.LIS.addLGT(pidsqid, genome1, genome2)
        
        f = open(outfile,'w')
        
        # get kmer scores!
        self.LIS.getClosestGID(False)
        self.LIS.printInfoHeader(f)
        self.LIS.printInfo(f)

    def grabLGTKmerInfo(self, kmer_file, pidsqid):
        # tmp array to store kmer data
        tmp = []
        
        with open(kmer_file) as fh:
            for l in fh:
                fields = l.rstrip().split("\t")
                if l[0:2] != "ID":
                    tmp.append([float(i) for i in fields[1:]])
            
            # average kmer columns
            tmer = np.mean(tmp, axis=0)
            self.LIS.addLGTTmer(pidsqid, tmer)
    
    def grabGenomeKmerInfo(self, kmer_file, genome):
        # tmp array to store kmer data
        tmp = []
        
        with open(kmer_file) as fh:
            for l in fh:
                fields = l.rstrip().split("\t")
                if l[0:2] != "ID":
                    tmp.append([float(i) for i in fields[1:]])
                
            # average kmer columns
            tmer = np.mean(tmp, axis=0)
            self.LIS.addGenomeTmer(genome, tmer)
                    
    def checkKmerFile(self, kmer_file):
        if "genome" in kmer_file:
            return True
        else: 
            return False

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
    KD = KmerData(args.paths_file)
    KD.kmerDataWrapper(args.kmer_dir,
                       args.outfile)
            
  

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kmer_dir', help="Directory containing kmer data.")
    parser.add_argument('paths_file', help="Path to file containing gid -> genome path.")
    parser.add_argument('-o','--outfile', default='kmer_scores_collated.csv', help="Output filename.")
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
