#!/usr/bin/env python
###############################################################################
#
# __parse_blast_output__.py - description!
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

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BlastData(object):
    def __init__(self, transfer_group_file):
        self.blast_data             = {}
        self.blast_data_per_pidsqid = {}
        self.vector_list            = ['cloning vector',
                                       'cosmid vector',
                                       'expression vector',
                                       'phix174',
                                       'shuttle vector',
                                       'vector']
        self.transfer_groups        = {}
        self.pidsqids               = {}
        self.GD                     = TFP.GroupData(transfer_group_file)
    
    def wrapper(self, blastfile, evalue_cutoff, print_by):
        # parse blast file
        self.parseBLASTFile(blastfile,
                            evalue_cutoff)
        
        # print out blast data
        if print_by.lower() == 'transfer_group':
            self.printOutBlastDataTG()
        elif print_by.lower() == 'pidsqid':
            self.printOutBlastDataPidsqid()
            #print self.blast_data_per_pidsqid
        
    
    def parseBLASTFile(self, blastfile, evalue_cutoff):
        with open(blastfile) as fh:
            for l in fh:
                if l[0:6] == "Query=":
                    pidsqid = l.rstrip().split(" ")[-1]
                    
                    # add pidsqid to TG
                    transfer_group = self.GD.group_membership[pidsqid]
                    self.addToTG(transfer_group)
                    
                    # add pidsqid to dict
                    
                    
                    # initialise dictionaries
                    if self.transfer_groups[transfer_group] == 1:
                        self.initiliseDict(transfer_group, 'total')
                        self.initiliseDict(transfer_group, 'vector')
                    self.initiliseDictPerPidsqid(transfer_group, pidsqid, 'total')
                    self.initiliseDictPerPidsqid(transfer_group, pidsqid, 'vector')
                    
                    for l in fh:
                        if '|' in l:
                            # i.e. gb|APEK01000004.1|
                            whitespace      = l.rstrip().split()
                            straightline    = whitespace[0].rstrip().split('|')
                            blast_hit       = straightline[1]
                            score           = float(whitespace[-2])
                            evalue          = float(whitespace[-1])
                            description     = " ".join(whitespace[1:-2])
                            
                            if evalue <= evalue_cutoff:
                                
                                # add to total 
                                self.buildDict(transfer_group, 'total')
                                self.buildDictPerPidsqid(transfer_group, pidsqid, 'total')
                                
                                # search through 'vector terms'
                                for search_term in self.vector_list:
                                    if search_term in description.lower():
                                        self.buildDict(transfer_group, 'vector')
                                        self.buildDictPerPidsqid(transfer_group, pidsqid, 'vector')
                                        break
                                    
                        elif l[0:5] == '*****':
                            # no hits! 
                            break
                        elif l[0:6] == "Query_":
                            break
                            
    def addToTG(self, TG):
        try:
            self.transfer_groups[TG] += 1 
        except KeyError:
            self.transfer_groups[TG] = 1
    
    def initiliseDictPerPidsqid(self, TG, pidsqid, attribute):
        try:
            self.blast_data_per_pidsqid[TG][pidsqid][attribute] = 0
        except KeyError:
            try:
                self.blast_data_per_pidsqid[TG][pidsqid] = {attribute:0}
            except KeyError:
                self.blast_data_per_pidsqid[TG] = {pidsqid:{attribute:0}}
    
    def buildDictPerPidsqid(self, TG, pidsqid, attribute):
        try:
            self.blast_data_per_pidsqid[TG][pidsqid][attribute] += 1 
        except KeyError:
            try:
                self.blast_data_per_pidsqid[TG][pidsqid][attribute] = 1
            except KeyError:
                try:
                    self.blast_data_per_pidsqid[TG][pidsqid] = {attribute:1}
                except KeyError:
                    self.blast_data_per_pidsqid[TG] = {pidsqid:{attribute:1}}
    
    def initiliseDict(self, TG, attribute):
        try:
            self.blast_data[TG][attribute] = 0
        except KeyError:
            self.blast_data[TG] = {attribute:0}
            
    def buildDict(self, TG, attribute):
        try:
            self.blast_data[TG][attribute] += 1 
        except KeyError:
            try:
                self.blast_data[TG][attribute] = 1
            except KeyError:
                self.blast_data[TG] = {attribute:1}
    
    def printOutBlastDataTG(self):
        for TG in self.blast_data.keys():
            total       = self.blast_data[TG]['total']
            vector      = self.blast_data[TG]['vector'] 
            if vector == 0:
                percentage = 0
            else:
                percentage  = float(vector)/total
            print "%s\t%d\t%d\t%f" % (TG,
                                      total,
                                      vector,
                                      percentage)               
                 
    def printOutBlastDataPidsqid(self):
        for TG in self.blast_data_per_pidsqid.keys():
            for pidsqid in self.blast_data_per_pidsqid[TG]:
                total       = self.blast_data_per_pidsqid[TG][pidsqid]['total']
                vector      = self.blast_data_per_pidsqid[TG][pidsqid]['vector'] 
                if vector == 0 or total == 0:
                    percentage = 0
                else:
                    percentage  = float(vector)/total
                print "%s\t%s\t%d\t%d\t%f" % (TG,
                                              pidsqid,
                                              total,
                                              vector,
                                              percentage) 

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
    BD = BlastData(args.transfer_group_file)
    BD.wrapper(args.blast_file,
               args.evalue_cutoff,
               args.print_by)
                
            
            

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file', help="")
    parser.add_argument('transfer_group_file', help="")
    parser.add_argument('-ec','--evalue_cutoff', type=float, default = 0, help="")
    parser.add_argument('-pb','--print_by', default = 'transfer_group', help="")
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
