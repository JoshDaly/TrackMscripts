#!/usr/bin/env python
###############################################################################
#
# __KO_definitions__.py - description!
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
import shlex
import subprocess
import time
#import os
#import errno
#import glob
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq
#import networkx as nx

# local imports
import query_kegg_db as QKDB

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KeggQueryRunner(object):
    def __init__(self):
        self.kegg_orthology     = {}
        self.cmds               = []
        self.result_list        = []
        self.kegg_query         = QKDB.KeggQuery()
    
    def wrapper(self, kegg_list, type, num_threads, outfile):
        self.parseKeggList(kegg_list)
        
        self.keggQuery(num_threads, outfile, type)
    
    def keggQuery(self, num_threads, outfile, type):
        count = 0
        
        pool = Pool(num_threads)
        
        async_objects = []
        
        # build commands
        for ko in self.kegg_orthology.keys():
            if type == 'test':
                if count > 100:
                    break
                count += 1 
            async_obj = pool.apply_async(work, args=(self.kegg_query,ko))
            async_objects.append(async_obj)
            
        # check commands
        while True:
            allready = True
            self.checkAsyncObj(async_objects, allready)
            if allready:
                break
            time.sleep(1)
        
        # open outfile
        of = open(outfile,'w')
        
        # get data
        for obj in async_objects:
            string = "\t".join([obj.get()[0],
                                obj.get()[1]])
            of.write("%s\n" % string)
    
    def checkAsyncObj(self, async_obj, check):
        for obj in async_obj:
                if not obj.ready():
                    check = False
                    break
    
    def parseKeggList(self, kegg_list):
        with open(kegg_list) as fh:
            for l in fh:
                ko = l.rstrip()
                self.kegg_orthology[ko] = 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def work(kegg_query, kegg_id):
    return kegg_query.queryDB(kegg_id)

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
    KQR = KeggQueryRunner()
    KQR.wrapper(args.KO_list,
                args.type,
                args.num_threads,
                args.outfile)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('KO_list', help="")
    parser.add_argument('-t','--type', help="")
    parser.add_argument('-num_threads','--num_threads', type=int, default=1, help="")
    parser.add_argument('-o','--outfile', help="Path to outfile")
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
