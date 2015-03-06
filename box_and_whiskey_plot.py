#!/usr/bin/env python
###############################################################################
#
# __box_and_whisker_plot__.py - description!
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
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import scipy.stats
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

class BoxWhiskerPlot(object):
    def __init__(self):
        pass
        
    def wrapper(self, coverage_file, xeno_contig_cutoff, type, outdir, outfile, outfmt):
        # grab data from coverage file
        CD = self.getCoverageData(coverage_file)
        
        # turn data into one long array
        boxplot_data = self.prepareData(CD, xeno_contig_cutoff, type)
        
        # create boxplot
        if type == 'coverage':
            self.makeBoxPlot(boxplot_data, type, outdir, outfile, outfmt)
        elif type == 'len_vs_coverage':
            self.makeBoxPlot(boxplot_data, type, outdir, outfile, outfmt)
            
    def getCoverageData(self, coverage_file):
        coverage_data = TFP.CoverageData(coverage_file)
        return coverage_data
    
    def prepareData(self,cov_data, xeno_contig_cutoff, type):
        data            = []
        xeno_contig     = [] 
        good_contigs    = []
        avg_coverage    = []
        low_coverage    = []
        high_coverage   = []
        average_coverage= []
        
        # What is a xenocontig? 
        # Based on Linkage
        # If below a specified cutoff, then it is a xenocontig!
        if type == 'coverage':
            for gid in cov_data.coverage_data.keys():
                for contig in cov_data.coverage_data[gid]:
                    coverage = cov_data.coverage_data[gid][contig][0]
                    if cov_data.coverage_data[gid][contig][1] > xeno_contig_cutoff:
                        good_contigs.append(coverage)
                    else:
                        xeno_contig.append(coverage)
                    
            data = [xeno_contig, good_contigs]        
        elif type == 'len_vs_coverage':
            for gid in cov_data.coverage_data.keys():
                for contig in cov_data.coverage_data[gid]:
                    coverage    = float(cov_data.coverage_data[gid][contig][0])
                    contig_len  = cov_data.coverage_data[gid][contig][2]
                    average_coverage.append(coverage)
                    if coverage <= 1.1 and coverage >= 0.9:
                        avg_coverage.append(contig_len) 
                    elif coverage >1.1:
                        high_coverage.append(contig_len)
                    elif coverage <0.9:
                        low_coverage.append(contig_len)
            average_coverage, coverage_std = self.averageArray(average_coverage)
            print "Average coverage: %f(+-%f)" % (average_coverage,coverage_std)
            print "<0.9: %d, 0.9-1.1: %d, >1.1: %d" % (len(low_coverage), len(avg_coverage), len(high_coverage))
            data = [low_coverage,avg_coverage,high_coverage]
        return data

    def averageArray(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)
    
    def makeBoxPlot(self, data, type, outdir, outfile, outfmt):
        #fig, ax = plt.subplots(figsize=(10,6))
        # initialise figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # fig titles
        #ax.set_ylabel('Coverage Relative to Largest Contig')
        
        if type == 'coverage':
            # box plot stats
            data1 = data[0]
            data2 = data[1]
            for i in data1:
                print i
            for i in data2:
                print i 
            z, p = scipy.stats.ttest_ind(data1,data2)
            #z, p = scipy.stats.mannwhitneyu(data1, data2)
            p_value = p
            #p_value = p * 2
            print '###### Pvalue ######'
            print str(p_value)
            print '####################'
            
            labels = ['Xeno Contigs','Nativus Contigs']
        
        elif type == 'len_vs_coverage':
            data1 = data[0]
            data2 = data[1]
            data3 = data[2]
            
            # low coverage vs average
            z, p = scipy.stats.mannwhitneyu(data1, data2)
            p_value = p * 2
            print "Summary Stats:"
            print 'low coverage vs average'
            print str(p_value)
            print '####################'
            # high coverage vs average
            z, p = scipy.stats.mannwhitneyu(data3, data2)
            p_value = p * 2
            print 'high coverage vs average'
            print str(p_value)
            print '####################'
            # low coverage vs high coverage 
            z, p = scipy.stats.mannwhitneyu(data1, data3)
            p_value = p * 2
            print 'low coverage vs high coverage'
            print str(p_value)
            print '####################'
            
            labels = ['Coverage <1', 'Coverage 1.0', 'Coverage >1']
            
        params = {
           'axes.labelsize': 8,
           'text.fontsize': 8,
           'legend.fontsize': 10,
           'xtick.labelsize': 10,
           'ytick.labelsize': 10,
           'text.usetex': False,
           'figure.figsize': [2.5, 4.5]
        }
        plt.rcParams.update(params)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.tick_params(axis='x', direction='out')
        ax.tick_params(axis='y', length=0)
        
        ax.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
        ax.set_axisbelow(True)
        
        bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)#, labels=labels)
        bp_colours = ['#FF0000', '#37B700']
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['fliers'], color='red', marker='+')
        ax.set_yscale('log')
        plt.ylabel('Contig Size (bp)')
        
        # write to file
        outfile = os.path.join(outdir, outfile)
        outfile = "%s.%s" % (outfile,outfmt)
        print outfile
        plt.savefig('%s' % outfile,format=outfmt)


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
    BWP = BoxWhiskerPlot()
    BWP.wrapper(args.coverage_file,
                args.xeno_contig_cutoff,
                args.type,
                args.outdir,
                args.outfile,
                args.outfmt)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('coverage_file', help="")
    parser.add_argument('-x','--xeno_contig_cutoff', type=float, default = 0.001, help="")
    parser.add_argument('-t','--type', default = 'coverage', help="")
    parser.add_argument('-o','--outdir', help="")
    parser.add_argument('-f','--outfile', default = 'boxplot.png', help="")
    parser.add_argument('-fmt','--outfmt', default = 'png', help="")
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
