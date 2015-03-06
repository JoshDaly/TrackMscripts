#!/usr/bin/env python
###############################################################################
#
# __contigsize_vs_linkage_scatterplot__.py - description!
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
import glob
#import os
#import errno
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq
#import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ScatterPlot():
    def __init__(self):
        self.bamm_cov_links_files   = {}
        self.bamm_data              = {}
        self.finished_genomes = ['A00001638',
                                 'A00003224',
                                 'A00003790',
                                 'A00004488',
                                 'A00004498']
    
    def wrapper(self, draft_cov_link_dir, finished_cov_link_dir, type):
        # grab bamm stats files
        self.getCovLinkFiles(draft_cov_link_dir)
        self.getCovLinkFiles(finished_cov_link_dir)
        
        # grab bamm cov links data
        self.getCovLinkData()
        
        # Plot data on scatter plot
        if type.lower() == 'len_vs_linkage':
            self.scatterPlotContigLenVsLinkage()
        
        elif type.lower() == 'cov_vs_linkage': 
            self.scatterPlotCoverageVsLinkage()
        
        
    def getCovLinkFiles(self, cov_link_dir):
        if cov_link_dir[-1] == '/':
            # remove last element
            cov_link_dir = cov_link_dir[0:-1]
            
        cov_link_files = glob.glob("%s/*/*cov_links_stats.csv" % cov_link_dir)
        for file in cov_link_files:
            gid = file.split("/")[-2]
            self.bamm_cov_links_files[gid] = file
    
    def getCovLinkData(self):
        for gid in self.bamm_cov_links_files.keys():
            cov_link_file = self.bamm_cov_links_files[gid]
            CLD = TFP.CovLinksData(cov_link_file)
            for contig in CLD.cov_link_data[gid]:
                try:
                    link_data = CLD.cov_link_data[gid][contig][2]
                    contig_len  = CLD.cov_link_data[gid][contig][0]
                    self.bamm_data[gid][contig] = [contig_len, link_data]
                except (IndexError, KeyError):
                    try:
                        self.bamm_data[gid] = {contig:[contig_len, link_data]}
                    except KeyError:
                        pass
    
    def scatterPlotCoverageVsLinkage(self):
    
    def scatterPlotContigLenVsLinkage(self):
        # initialise figure
        #fig = plt.figure(figsize=(20,10),dpi=300)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # initialise coords
        x   = [] # contigsize
        y   = [] # linkage
        xf  = []
        xy  = []
        
        # initialise colour array
        colours     = []
        coloursf    = []
        
        # build plot data
        x,y,colours, xf, yf, coloursf = self.buildPlotData()
        
        plt.scatter(x,y, s=10, alpha=0.5, c=colours, lw=0)
        plt.scatter(xf,yf, s=50, alpha=1, c=coloursf, lw=0)
        
        # line of best fit (regression)
        xl, yl, Rsqr    = self.getLineOfRegression(x, y)
        xfl, yfl, fRsqr = self.getLineOfRegression(xf, yf)
        
        # plot line of best fit
        plt.plot(xl, yl, '-k')
        plt.plot(xfl, yfl, '--k')
    
        plt.text(.9*max(x)+.1*min(x),.9*max(y)+.1*min(y),'$R^2 = %0.2f$'% Rsqr, fontsize=20)
        plt.text(.9*max(xf)+.1*min(xf),.9*max(yf)+.1*min(yf),'$R^2 = %0.2f$'% fRsqr, fontsize=20)
        
        plt.xlabel('Contig Length (bp)')
        plt.ylabel('Total Linkage')
        
        # set axis to logarithmic
        ax.set_xscale('log')
        #ax.set_yscale('log')
        
        plt.show()
        
    def getLineOfRegression(self, x, y):
        # line of best fit (regression)
        par  = np.polyfit(x, y, 1, full=True)
        
        slope       = par[0][0]
        intercept   = par[0][1]
        xl = [min(x), max(x)]
        yl = [slope*xx + intercept  for xx in xl]
        
        # coefficient of determination, plot text
        variance = np.var(y)
        residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(x,y)])
        Rsqr = np.round(1-residuals/variance, decimals=2)
        #plt.text(.9*max(xd)+.1*min(xd),.9*max(yd)+.1*min(yd),'$R^2 = %0.2f$'% Rsqr, fontsize=30)

        return xl, yl, Rsqr
    
    def buildPlotData(self):
        x           = []
        y           = []
        colours     = []
        xf          = []
        yf          = []
        coloursf    = []
        
        for gid in self.bamm_data.keys():
            for contig in self.bamm_data[gid]:
                contig_len = int(self.bamm_data[gid][contig][0])
                coverage   = float(self.bamm_data[gid][contig][1])

                if gid in self.finished_genomes:
                    xf.append(contig_len)
                    yf.append(coverage)
                    coloursf.append('#377eb8')
                else:
                    x.append(contig_len)
                    y.append(coverage)
                    colours.append('#4daf4a')
        
        return x,y,colours, xf, yf, coloursf

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
    SP = ScatterPlot()
    SP.wrapper(args.draft_cov_link_dir,
               args.finished_cov_link_dir)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('draft_cov_link_dir', help="")
    parser.add_argument('finished_cov_link_dir', help="")
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
