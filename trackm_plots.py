#!/usr/bin/env python
###############################################################################
#
# __trackm_plots__.py - description!
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
import numpy as np
np.seterr(all='raise')
import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Plot(object):
    def __init__(self):
        pass
    
    def wrapper(self, type, transfer_group_data, outfile, outfmt):
        if type == 'transfer_group_bar':
            # get data
            data = self.getTransferGroupBarData(transfer_group_data)
            
            # plot data
            self.barPlot(data, outfile, outfmt)
            
    def getTransferGroupBarData(self, transfer_group_data):
        data = TFP.TransferGroupSummaryData(transfer_group_data)
        return data
    
    def barPlot(self, data, outfile, outfmt):
        # covert data to coords
        x, y = self.makeCoordsBar(data.transfer_group_data)
        
        # sort arrays by y
        x, y = self.sortArray(x,y)
        
        # initialise graph
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # necessary variables
        width = 0.35
        
        # the bars
        plt.bar(x,
                y,
                width,
                color='blue',
                linewidth = 0
                )
        
        # label, title and axis ticks
        ax.set_xlim(-width-3,len(x)+width)
        ax.set_ylim(-1,max(y)+10)
        ax.set_ylabel('Number of Pidsqids')
        #ax.set_title('Pidsqid Distribution Amongst Transfer Groups')
        plt.tight_layout()
        
        # output directory/file
        print 'Writing to file %s' % (outfile)
        plt.savefig(outfile,format=outfmt)
        plt.close()

    def sortArray(self, x, y):
        x = np.array(x)
        y = np.array(y)
        sorted_x = []
        sorted_y = []
        
        # get indices for sorted array
        indices = np.argsort(y)
        count = 0 
        for index in indices:
            sorted_x.append(count)
            sorted_y.insert(0,y[index])
            count+= 1 
        
        return sorted_x, sorted_y
    
    def makeCoordsBar(self, data):
        x = []
        y = []
        
        for tg in data.keys():
            pidsqids = int(data[tg][7])
            x.append(int(tg))
            y.append(pidsqids)
        
        return x,y

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
    P = Plot()
    P.wrapper(args.type,
              args.transfer_group_data,
              args.outfile,
              args.outfmt)
                


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-tgb','--transfer_group_data', help="")
    parser.add_argument('-t','--type', help="transfer_group_bar")
    parser.add_argument('-fmt','--outfmt', type=str, default = 'png',help="")
    parser.add_argument('-o','--outfile', help="")
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
