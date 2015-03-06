#!/usr/bin/env python
###############################################################################
#
# __contig_size_distribution_plot__.py - description!
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
import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigData(object):
    def __init__(self, hitdata_file):
        self.hitdata                = TFP.HitData(hitdata_file)
        self.contig_len             =  self.hitdata.contig_size
        self.contig_hits            = {}
        self.contig_hits_no_ecoli   = {} 
        self.window                 = []

    def wrapper(self):
        # get longest contig 
        longest_contig = self.getLongestContig()
        
        # count number of hits per 50bp window
        self.countNumberOfHitsPer50bpWindow(longest_contig)

    def getLongestContig(self):
        longest_contig = 0
        for pidsqid in self.contig_len.keys():
            if self.contig_len[pidsqid][0] > longest_contig:
                longest_contig = self.contig_len[pidsqid][0]
            elif self.contig_len[pidsqid][1] > longest_contig:
                longest_contig = self.contig_len[pidsqid][1]
                
        return longest_contig
            
    def countNumberOfHitsPer50bpWindow(self, longest_contig):
        # break largest contig into 50 bp window
        for i in range(500,longest_contig+1,50):
            self.window.append(i)
            self.contig_hits[i] = 0
            self.contig_hits_no_ecoli[i] = 0
            if i >= 10000:
                break
            
        # loop through hits
        for window in self.contig_hits.keys():
            for pidsqid in self.contig_len.keys():
                if self.hitdata.intra_or_inter[pidsqid] == 'inter':
                    contig_len_1 = self.contig_len[pidsqid][0]
                    contig_len_2 = self.contig_len[pidsqid][1]
                    
                    if contig_len_1 >= window and contig_len_1 <= (window + 49):
                            self.contig_hits[window] += 1 
                            
                    if contig_len_2 >= window and contig_len_2 <= (window + 49):
                        self.contig_hits[window] += 1
                    
                    if self.hitdata.genus[pidsqid][0] != 'Escherichia' and self.hitdata.genus[pidsqid][1] != 'Escherichia':
                    
                        if contig_len_1 >= window and contig_len_1 <= (window + 49):
                            self.contig_hits_no_ecoli[window] += 1 
                            
                        if contig_len_2 >= window and contig_len_2 <= (window + 49):
                            self.contig_hits_no_ecoli[window] += 1 
                    
            
class Plot(object):
    def __init__(self, contig_data):
        self.contig_data = contig_data
    
    def wrapper(self):
        # turn data into x,y
        x,y,ye = self.convertDataToCoords()
        
        # frequency plot
        self.frequencyPlot(x,y,ye)
        
    def convertDataToCoords(self):
        y   = []
        ye  = []
        x   = self.contig_data.window
        for i,window in enumerate(self.contig_data.window):
            y.append(self.contig_data.contig_hits[window])
            ye.append(self.contig_data.contig_hits_no_ecoli[window])
        return x,y,ye

    def frequencyPlot(self,x,y,ye):
        plt.scatter(x, y, marker='|')
        plt.scatter(x, ye, marker='|')
        plt.plot(x, y, linestyle='-')
        plt.plot(x, ye, linestyle='-')
            
        #plt.plot(x,y,c='#FFFFFF')
        #plt.fill_between(x, y, 1e-6, facecolor = '#C0C0C0')
        
        plt.show()
        
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
    CD = ContigData(args.hitdata_file)
    CD.wrapper()
    P = Plot(CD)
    P.wrapper()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata_file', help="")
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
