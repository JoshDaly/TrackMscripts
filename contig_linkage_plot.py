#!/usr/bin/env python
###############################################################################
#
# __contig_linkage_plot__.py - description!
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
import matplotlib.pyplot as plt
import os
import errno
import numpy as np
np.seterr(all='raise')
import math as math
import glob
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

class Plot(object):
    def __init__(self):
        self.bamm_files     = {}
        self.blast_files    = {}
        #self.BD             = TFP.BamMData(bamm_links_file, bamm_cov_file)
        #if blast_file:
        #    self.BFP        = TFP.BlastFileParser(blast_file)
        self.plot_data      = {}
        self.contig_len     = {}
        self.coords         = {}
        self.contig_linkage = {}
    
    def wrapper(self, bamm_dir, blast_dir, outdir, type, force, cov_file, links_file):
        # get bamm_links, bamm_cov and blast files
        if bamm_dir:
            self.getBammFiles(bamm_dir)
        if blast_dir:
            self.getBLASTFiles(blast_dir)
        if cov_file and links_file:
            self.getIndividualBammFiles(cov_file, links_file)    
        
        # loop through files
        for gid in self.bamm_files.keys():
            
            if type.lower()=='scatter':
            
                if self.doesXYFileExist(gid, outdir, force):
                    bamm_links_file = self.bamm_files[gid]['links']
                    bamm_cov_file   = self.bamm_files[gid]['coverages']
                    if blast_dir:
                        blast_file      = self.blast_files[gid]
                    
                    try:
                        BD  =  TFP.BamMData(bamm_links_file, bamm_cov_file)
                        if blast_dir:
                            BFP =  TFP.BlastFileParser(blast_file)
                        else:
                            BFP = False
                        
                        print gid
                    
                        # get plot data
                        self.buildPlotData(BD)
                        
                        # try to run xyscatter
                        self.checkXYScatter(outdir, gid, BD, BFP, force)
                        
                    except (ValueError,IndexError):
                        pass
                
            elif type.lower()=='bar':
                bamm_links_file = self.bamm_files[gid]['links']
                bamm_cov_file   = self.bamm_files[gid]['coverages']
                blast_file      = self.blast_files[gid]
                        
                BD  =  TFP.BamMData(bamm_links_file, bamm_cov_file)
                BFP =  TFP.BlastFileParser(blast_file)
                    
                # get plot data
                self.buildPlotData(BD)
                
                # build bar chart for each contig
                for contig in self.plot_data.keys():
                    self.barChart(contig, outdir)
        
    def getIndividualBammFiles(self,cov_file, links_file):
        gid = cov_file.split("/")[-1].split("_")[0]
        print gid
        self.bamm_files[gid] = {'coverages':cov_file}
        self.bamm_files[gid]['links'] = links_file
    
    def checkXYScatter(self, outdir, gid, BD, BFP, force):
        try:
            self.xyScatter(outdir, gid, BD, BFP, force)
        except ZeroDivisionError:
            print 'No coverage data available for %s' % gid
            pass
            
    def doesXYFileExist(self, gid, outdir, force):
        # check if file already exists
            outfile = os.path.join(outdir,"%s/%s.scatter.ylog.png" % (gid,gid))
            if force:
                print 'Writing to file %s' % (outfile)
                return True
            else:
                if os.path.isfile(outfile):
                    print "File %s already exists, set force=True to overwrite" % (outfile)
                    return False
                else:
                    print 'Writing to file %s' % (outfile)
                    return True
    
    def getBammFiles(self, bamm_dir):
        if bamm_dir[-1] == '/':
            # remove last element
            bamm_dir = bamm_dir[0:-1]
            
        bamm_files = glob.glob("%s/*/*" % bamm_dir)
        for file in bamm_files:
            gid = file.split("/")[-2]
            if 'coverages' in file:
                try:
                    self.bamm_files[gid]['coverages'] = file
                except KeyError:
                    self.bamm_files[gid] = {'coverages':file}
            if 'links.tsv' in file:
                try:
                    self.bamm_files[gid]['links'] = file
                except KeyError:
                    self.bamm_files[gid] = {'links':file}
    
    def getBLASTFiles(self, blast_dir):
        if blast_dir[-1] == '/':
            # remove last element
            blast_dir = blast_dir[0:-1]
        
        blast_files = glob.glob("%s/*/*_vs_nr.blast20151802" % blast_dir)
        for file in blast_files:
            gid = file.split("/")[-2]
            self.blast_files[gid] = file 
            
    def buildPlotData(self, bamm_data):
        # loop through contigs
        for contig1 in bamm_data.bamm_links.keys():
            for contig2 in bamm_data.bamm_links[contig1]:
                for link in bamm_data.bamm_links[contig1][contig2]:
                    try:
                        pos = link[2]
                        try:
                            self.plot_data[contig1][pos] += 1 
                        except KeyError:
                            try:
                                self.plot_data[contig1][pos] = 1
                            except KeyError:
                                self.plot_data[contig1]= {pos:1}
                    except KeyError:
                        pass
                
    def xyScatter(self, outdir, gid, bamm_data, blast_data, force):
        
        # get the longest contig
        longest_contig_len  = 0
        longest_contig      = ''
        most_coverage       = 0
        chromosome_links    = 0
        total_links         = 0
        
        # initialise figure
        fig = plt.figure(figsize=(20,10),dpi=300)
        ax = fig.add_subplot(111)
        
        # initialise coords
        x = [] # contigsize
        y = [] # coverage
        
        # initialise colour array
        colours = []
        
        # contig size vs coverage + linkage
        for contig in bamm_data.bamm_cov_data.keys():
            # colour contigs by blasthit
            if blast_data:
                colour = self.colourByAnnotation(blast_data,contig)
                colours.append(colour)
            
            if int(bamm_data.contig_len[contig]) > longest_contig_len:
                longest_contig_len  = int(bamm_data.contig_len[contig])
                longest_contig      = contig
            x.append(int(bamm_data.contig_len[contig]))
            y.append(self.returnAverage(bamm_data.bamm_cov_data[contig]))
            self.coords[contig] = [int(bamm_data.contig_len[contig]), self.returnAverage(bamm_data.bamm_cov_data[contig])]
            
        for i,v in enumerate(y):
            chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
            y[i] = float(v)/float(chromosome_coverage)
            if y[i] > most_coverage:
                most_coverage = y[i]
        
        # Grab longest contigs link count
        #for contig in bamm_data.bamm_links_data[longest_contig]:
        #    links = bamm_data.bamm_links_data[longest_contig][contig]
        #    chromosome_links += links
        
        # loop through contigs, and link with lines weighted by linkage
        for contig1 in bamm_data.bamm_links_data.keys():
            for contig2 in bamm_data.bamm_links_data[contig1]:
                chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
                total_links += bamm_data.bamm_links_data[contig1][contig2]
                linkage = math.sqrt(bamm_data.bamm_links_data[contig1][contig2])/25
                coords1 = self.coords[contig1]
                coords2 = self.coords[contig2]
                x2 = []
                y2 = []
                x2.append(coords1[0])
                x2.append(coords2[0])
                y2.append(coords1[1]/float(chromosome_coverage))
                y2.append(coords2[1]/float(chromosome_coverage))
                plt.plot(x2,y2,'k-',lw=linkage, alpha=0.5)
        
        # print out data
        try:
            self.writeStatsDataToFileWrapper(bamm_data, gid, force, outdir, longest_contig, total_links)
        except ZeroDivisionError:
            stats_file = os.path.join(outdir, '%s/%s_cov_links_stats.csv' % (gid,gid))
            print 'No coverage data available for %s' % stats_file
            pass
        
        if blast_data:
            plt.scatter(x,y, s=200, alpha=0.5, c=colours)
        else:
            plt.scatter(x,y, s=200, alpha=0.5)
            
        # set axis limits
        #plt.ylim(-0.2, most_coverage+1)
        #plt.xlim(-(longest_contig_len/10), longest_contig_len+(longest_contig_len/10))
        
        # set axis labels
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Coverage Relative to Largest Contig (%)')
        
        # set axis to logarithmic
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        # output directory/file
        outfile = os.path.join(outdir,"%s/%s.scatter.ylog.png" % (gid,gid))
        plt.savefig(outfile,format='png')
        plt.close()
                
    def writeStatsDataToFileWrapper(self, bamm_data, gid, force, outdir, longest_contig, total_links):
        stats_file = os.path.join(outdir, '%s/%s_cov_links_stats.csv' % (gid,gid))
    
        if force:
            # print header 
            self.writeStatsDataToFile(bamm_data, stats_file, longest_contig, total_links)
            
        else:
            if os.path.isfile(stats_file):
                print "File %s already exists, set force=True to overwrite" % (stats_file)
                    
            else:
                self.writeStatsDataToFile(bamm_data, stats_file, longest_contig, total_links)
                
    def writeStatsDataToFile(self, bamm_data, stats_file, longest_contig, total_links):
        chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
        try:
            chromosome_linkage  = self.returnAverage(bamm_data.bamm_links_total[longest_contig])
            # print header 
            f = open(stats_file, 'w')
            f.write("\t".join(['contig',
                               'contig_len',
                               'contig_cov',
                               'contig_links',
                               'contig_links_div_chromo',
                               'contig_links_div_total\n']))
            
            for contig in bamm_data.bamm_cov_data.keys():
                contig_cov              = str(self.returnAverage(bamm_data.bamm_cov_data[contig])/float(chromosome_coverage))
                contig_len              = str(bamm_data.contig_len[contig])
                contig_links            = '0'
                contig_links_div_total  = '0'
                contig_links_div_chromo = '0'
                try:
                    contig_links            = str(bamm_data.bamm_links_total[contig])
                    contig_links_div_total  = str(bamm_data.bamm_links_total[contig]/float(total_links))
                    contig_links_div_chromo = str(bamm_data.bamm_links_total[contig]/float(chromosome_linkage))
                except KeyError:
                    pass
                f.write("\t".join([contig,
                                   contig_len,
                                   contig_cov,
                                   contig_links,
                                   contig_links_div_chromo,
                                   "%s\n" % contig_links_div_total ]))
            
        except KeyError:
            # print header 
            f = open(stats_file, 'w')
            f.write("\t".join(['contig',
                               'contig_len',
                               'contig_cov\n']))
            for contig in bamm_data.bamm_cov_data.keys():
                contig_cov              = str(self.returnAverage(bamm_data.bamm_cov_data[contig])/float(chromosome_coverage))
                contig_len              = str(bamm_data.contig_len[contig])
                f.write("\t".join([contig,
                                   contig_len,
                                   "%s\n"  % contig_cov]))
            
    
    def colourByAnnotation(self, blast_data, contig):
        if blast_data.topHits[contig] == 'MuyBueno':
            return 'green'
        elif blast_data.topHits[contig] == 'xeno':
            return 'red'
        elif blast_data.topHits[contig] == 'NH':
            return 'grey'
        
    def returnAverage(self,array):
        a = np.array(array)
        return np.mean(a)
    
    def barChart(self, contig, outdir):
        contig_len = int(self.BD.contig_len[contig])
        
        # covert data to coords
        x, y = self.makeCoords(contig, contig_len)
        
        # initialise graph
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # necessary variables
        index = np.arange(contig_len)
        width = 0.35
        
        # the bars
        plt.bar(index,
                y,
                width)
        
        # label, title and axis ticks
        ax.set_xlim(-width,len(index)+width)
        ax.set_ylim(0,max(y))
        ax.set_ylabel('Links')
        ax.set_title('Links distribution across contig')
        plt.tight_layout()
        
        # output directory/file
        outfile = os.path.join(outdir,"%s.linkage_profile.png" % (contig))
        print 'Writing to file %s' % (outfile)
        plt.savefig(outfile,format='png')
        plt.close()
        
    def makeCoords(self, contig, contig_len):
        #coords
        x = []
        y = []
        
        for i in range(1,contig_len+1):
            x.append(str(i))
            try:
                y.append(self.plot_data[contig][str(i)])
            except KeyError:
                y.append(0)
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
    P.wrapper(args.bamm_dir,
              args.blast_dir,
              args.outdir,
              args.type,
              args.force,
              args.bamm_cov,
              args.bamm_links) 
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-bd','--bamm_dir', default=False, help="Bamm output directory containing files in fmt: A00001235/bamm_cov, bamm_links")
    parser.add_argument('outdir', help="")
    parser.add_argument('-b','--blast_dir', default=False, help="Directory containing blast output in fmt: A00001235/blast.file")
    parser.add_argument('-t','--type', default='scatter',help="Set display type: bar or scatter")
    parser.add_argument('-f','--force', default=False,help="Set force=True to overwrite file")
    parser.add_argument('-bc','--bamm_cov', default=False,help="Path to bamm coverage file")
    parser.add_argument('-bl','--bamm_links', default=False,help="Path to bamm links file")
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
