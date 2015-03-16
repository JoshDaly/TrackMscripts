#!/usr/bin/env python
###############################################################################
#
# __parse_blastoutput_contigs__.py - description!
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
import os
import errno
#import numpy as np
#np.seterr(all='raise')
#import matplotlib as mpl
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
#from Bio import SeqIO
#from Bio.Seq import Seq
#import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP
from cb2cols import Cb2Cols as CB2

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Contigs(object):
    def __init__(self, taxonomy_file, group_file, hitdata):
        self.TD                 = TFP.TaxonomyData(taxonomy_file)
        self.GD                 = TFP.GroupData(group_file)
        self.HD                 = TFP.HitData(hitdata)
        self.cb2                = CB2()
        self.blast_data         = {}
        self.isVector           = {}
        self.contig_len         = {}
        self.count              = 0
        self.transfer_groups    = {} 
        self.contam_contigs     = {}
    
    def wrapper(self, blastdir, type, outfile, contaminated_contigs_file):
        # get contaminated contigs
        self.getContaminatedContigs(contaminated_contigs_file)
        
        # read in blast dir, grab blast files
        self.readBLASTDir(blastdir)
        
        if type.lower() == 'pie_chart':
            # create pie chart
            self.plotPieChart()
        
        elif type.lower() == 'tg_status':
            self.returnTGStatus(outfile)
        
    def getContaminatedContigs(self, contaminated_contigs_file):
        with open(contaminated_contigs_file) as fh:
            for l in fh:
                tabs    = l.rstrip().split("\t")
                gid     = tabs[0]
                contig  = tabs[1]
                try:
                    self.contam_contigs[gid][contig] = 1 
                except KeyError: 
                    self.contam_contigs[gid] = {contig:1} 
    
    def readBLASTDir(self, blastdir):
        blast_search = os.path.join(blastdir,"*/*vs_nr.blast20151802")
        blast_files = glob.glob(blast_search)
        
        count = 0
        
        for blast_file in blast_files:
            gid = blast_file.split("/")[-2]
            
            count += 1 
            
            # check gid
            if gid in self.contam_contigs:
            
                # remove E.coli from analysis
                if self.TD.taxon_genus[gid] != 'Escherichia':
                    
                    BFP = TFP.BlastFileParser(blast_file)
                    
                    for query_id in BFP.isvector.keys():
                        contig_id = self.getContigID(query_id)
                        if contig_id in self.contam_contigs[gid]:
                            self.isVector[contig_id]     = BFP.isvector[query_id]
                    
                    for query_id in BFP.topHits2.keys():
                        contig_id = self.getContigID(query_id)
                        if contig_id in self.contam_contigs[gid]:
                            self.blast_data[contig_id]   = BFP.topHits2[query_id]
                        
                    for query_id in BFP.contig_len.keys():
                        contig_id = self.getContigID(query_id)
                        if contig_id in self.contam_contigs[gid]:
                            self.contig_len[contig_id]   = BFP.contig_len[query_id]
    
            if count > 2000:
                break
    
    def getContigID(self, query_id):
        underscore = query_id.split("_")
        contigid   = "_".join(underscore[1:])
        return contigid
                    
    def plotPieChart(self):
        # get pie chart data
        labels,sizes = self.getPieChartData()
        
        col_set = "qualSet3"
        ColBrewColours = self.cb2.maps[col_set].values()[0:len(labels)+1]
        print ColBrewColours
        
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=ColBrewColours)
        
        plt.axis('equal')
        
        plt.show()
        
    def getPieChartData(self):
        common_contams = ['pseudomonas', 'propionibacterium', 'Streptococcus', 'Micrococcus', 'Stenotrophomonas', 'Xanthomonas', 'Pseudoxanthomonas', 'Burkholderia', 'Deinococcus', 'Corynebacterium']
        
        piechartdata            = {}
        piechartdataCorrected   = {}
        total                   = 0
        labels                  = []
        sizes                   = []
        for contig in self.blast_data.keys():
            if self.contig_len[contig] >= 500:
                total += 1 
                description = self.blast_data[contig][0]
                evalue      = self.blast_data[contig][1]
                if  self.isVector[contig] > 0:
                    description = 'Vector'
                    self.addData(description, piechartdata)
                elif 'escherichia' in description.lower() or 'e. coli' in description.lower():
                    description = 'Escherichia'
                    self.addData(description, piechartdata)
                elif 'uncultured' in description.lower():
                    description = 'Uncultured'
                    self.addData(description, piechartdata)
                elif 'streptomyces' in description.lower():
                    description = 'Streptomyces'
                    self.addData(description, piechartdata)
                elif 'staphylococcus' in description.lower():
                    description = 'Staphylococcus'
                    self.addData(description, piechartdata)
                elif 'mycobacterium' in description.lower():
                    description = 'Mycobacterium'
                    self.addData(description, piechartdata)
                elif 'bacteroides' in description.lower():
                    description = 'Bacteroides'
                    self.addData(description, piechartdata)
                for contam in common_contams:
                    if contam.lower() in description.lower():
                        description = 'Common contaminants'
                        self.addData(description, piechartdata)
                else:
                    self.addData(description, piechartdata)
        
        for description in  piechartdata.keys():
            count = piechartdata[description]
            percentage = count/float(total)
            if percentage > 0.01: 
                self.addCorrectedData(description, count, piechartdataCorrected)
            else:
                description = 'other'
                self.addCorrectedData(description, count, piechartdataCorrected)
                
        for description in  piechartdataCorrected.keys():
            count = piechartdataCorrected[description]
            labels.append(description)
            sizes.append(count)
        
        return labels, sizes
    
    def addCorrectedData(self, description, count, dict):
        try:
            dict[description] += count
        except KeyError:
            dict[description] = count
    
    def addData(self, description, dict):
        try:
            dict[description] += 1 
        except KeyError:
            dict[description] = 1 

    def returnTGStatus(self, outfile):
        for TG in self.GD.group_data.keys():
            total = 0
            contam= 0 
            for pidsqid in self.GD.group_data[TG]:
                total += 1
                if self.checkPidsqid(pidsqid):
                    # pidsqid OK
                    pass
                else:
                    contam += 1 
            percentage = contam / float(total)
            self.transfer_groups[TG] = percentage
        
        # intialise outfile and print out TG status data
        outfile = open(outfile, 'w')
        
        for TG in self.transfer_groups.keys():
            outfile.write("%s\t%f\n" % (TG, self.transfer_groups[TG]))
                    
    def checkPidsqid(self, pidsqid):
        contig1 = self.HD.contig_name[pidsqid][0]
        contig2 = self.HD.contig_name[pidsqid][1]
        if contig1 in self.blast_data or contig2 in self.blast_data:
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
    #p = Popen(cmd.split(' '), stdout=PIPE)
    #return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    C = Contigs(args.taxonomy_file,
                args.group_file,
                args.hitdata)
    C.wrapper(args.blastdir,
              args.type,
              args.outfile,
              args.contaminated_contigs_file)
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('contaminated_contigs_file', help="")
    parser.add_argument('blastdir', help="")
    parser.add_argument('taxonomy_file', help="")
    parser.add_argument('group_file', help="")
    parser.add_argument('hitdata', help="")
    parser.add_argument('-t','--type', default = 'Pie_chart',help="Set the type of analysis to be completed. (default)Pie_chart: pie chart of specified contigs. tg_status: returns percentage of TGs contaminated.")
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
