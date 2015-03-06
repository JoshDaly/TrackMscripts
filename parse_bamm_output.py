#!/usr/bin/env python
###############################################################################
#
# __parse_bamm_output__.py - description!
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
import numpy as np
np.seterr(all='raise')
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

class BamM(object):
    def __init__(self, hitdata_file):
        self.HD                     = TFP.HitData(hitdata_file)
        self.bamm_files             = {}
        self.contig_data            = {}
        self.transfer_group_data    = {}
        self.coverage_data          = {}
        self.active_contigs         = {}
        self.average_coverage       = {}
        self.std_coverage           = {}
        self.average_links          = {}
        self.std_links              = {}
        
        
    def wrapper(self, bamm_dir):       
        # grab files 
        self.getBammFiles(bamm_dir)
        
        # loop through files and get data
        self.getBammData()
        
        # get coverage and linkage average
        self.getAverageCoverageLinks()
        
        self.printOutData()
        
    def getBammFiles(self, bamm_dir):
        files = glob.glob("%s/*/*" % bamm_dir)
        for file in files:
            gid = file.split("/")[-2]
            if 'coverages' in file:
                try:
                    self.bamm_files[gid]['coverages'] = file
                except KeyError:
                    self.bamm_files[gid] = {'coverages':file}
            if 'links' in file:
                try:
                    self.bamm_files[gid]['links'] = file
                except KeyError:
                    self.bamm_files[gid] = {'links':file}
                
    def getBammData(self):
        for gid in self.bamm_files.keys():
            
            # get bamm files
            coverages_file = self.bamm_files[gid]['coverages'] 
            links_file = self.bamm_files[gid]['links']
            
            # parse BamM files
            BammData = TFP.BamMData(links_file, coverages_file)
            
            # parse BamM data
            self.getCoverageData(BammData, gid)
            self.getLinkData(BammData, gid)
    def getAverageCoverageLinks(self):
        for gid in self.contig_data.keys():
            contig_coverages = []
            contig_linkages = []
            for contig in self.contig_data[gid]:
                genus           = self.HD.gid_to_genus[gid]
                phylum          = self.HD.genus_to_phylum[genus]
                contig_length   = self.contig_data[gid][contig]['contig_length']
                coverages        = self.contig_data[gid][contig]['coverages']
                try:
                    links           = self.contig_data[gid][contig]['links']
                except KeyError:
                    links           = 0
                
                # build arrays
                contig_coverages.append(coverages)
                contig_linkages.append(links)

                TGs = self.getTransferGroupsByContig(contig)
                
                if TGs != '':
                    try:
                        self.active_contigs[gid][contig] = 1
                    except KeyError:
                        self.active_contigs[gid] = {contig:1}

            average_coverage, std_coverage = self.averageArray(contig_coverages)
            average_linkage, std_linkage  = self.averageArray(contig_linkages)
            
            self.average_coverage[gid]  = average_coverage
            self.std_coverage[gid]      = std_coverage
            self.average_links[gid]     = average_linkage
            self.std_links[gid]         = std_linkage
            
                      
    def averageArray(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)            
    
    def getLinkData(self, BammData, gid):
        for contig in BammData.bamm_links_total.keys():
            links = BammData.bamm_links_total[contig]
            self.addContigInfo(gid,
                               contig,
                               'links',
                               links)
    
    def getCoverageData(self,BammData, gid): 
        for contig in BammData.bamm_cov_data.keys():
            contig_length       = BammData.contig_len[contig]
            cov_array           = BammData.bamm_cov_data[contig]
            average_coverage    = self.averageArray(cov_array)[0]
            
            self.addContigInfo(gid,
                               contig,
                               'contig_length',
                               contig_length)
            self.addContigInfo(gid,
                               contig,
                               'coverages',
                               average_coverage)
    
    def addContigInfo(self, gid, contig, attribute, value):
        try:
            self.contig_data[gid][contig][attribute] = value
        except KeyError:
            try:
                self.contig_data[gid][contig] = {attribute:value}
            except KeyError:
                    self.contig_data[gid] = {contig:{attribute:value}}
    
    def printOutData(self):
        # print header
        print "\t".join(["gid",
                         "phylum",
                         "genus",
                         "contig_id",
                         "contig_length",
                         "contig_coverage",
                         "avg_genome_coverage",
                         "std_genome_coverage",
                         "contig_links",
                         "avg_genome_links",
                         "std_genome_links",
                         "transfer_groups"])
        
        for gid in self.active_contigs.keys():
            for contig in self.active_contigs[gid]:
                try:
                    links           = str(self.contig_data[gid][contig]['links'])
                except KeyError:
                    links           = '0'
                genus           = self.HD.gid_to_genus[gid]
                phylum          = self.HD.genus_to_phylum[genus]
                contig_length   = self.contig_data[gid][contig]['contig_length']
                coverages       = str(self.contig_data[gid][contig]['coverages'])
                avg_coverages   = str(self.average_coverage[gid])
                std_coverages   = str(self.std_coverage[gid])
                avg_links       = str(self.average_links[gid])
                std_links       = str(self.std_links[gid])
                
                
                TGs = self.getTransferGroupsByContig(contig)
                print "\t".join([gid,
                                 phylum,
                                 genus,
                                 contig,
                                 str(contig_length),
                                 coverages,
                                 avg_coverages,
                                 std_coverages,
                                 links,
                                 avg_links,
                                 std_links,
                                 TGs])
    
    def getTransferGroupsByContig(self, contig):
        TGs = ''
        try:
            for TG in self.HD.transfer_group_to_contigs[contig]:
                if TG != 'intra_phyla':
                    TGs += '%s;' % TG
        except KeyError:
            pass
        return TGs
            

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
    B = BamM(args.hitdata)
    B.wrapper(args.bamm_dir)
                


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('bamm_dir', help="")
    parser.add_argument('hitdata', help="")
    #parser.add_argument('bamm_coverage_file', help="")
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
