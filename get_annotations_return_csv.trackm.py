#!/usr/bin/env python
###############################################################################
#
# __get_annotations_return_csv__.py - description!
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


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TaxonomyFileParser(object):
    def __init__(self, l):
        self.readTaxonFile(l)
        
    def readTaxonFile(self,l):
        tabs = l.rstrip().split("\t")
        self.gid            = tabs[0]
        self.taxonomy       = tabs[1]
        self.organism       = tabs[2]
        
class TaxonomyData(object):
    def __init__(self, taxon_data):
        self.taxon_genus    = {}
        self.taxon_phylum   = {}
        self.taxon_string   = {}
        self.buildTaxonData(taxon_data)
        
    def buildTaxonData(self, taxon_data):
        with open(taxon_data) as fh:
            for l in fh:
                TFP = TaxonomyFileParser(l)
                TFP.readTaxonFile(l)
                
                # phylum-level taxonomy
                phylum  = TFP.taxonomy.split(";")[1][1:]
                self.taxon_phylum[TFP.gid]  = phylum
                
                # try to grab genus-level taxonomy
                if "g__" in TFP.taxonomy:
                    taxons = TFP.taxonomy.split()
                    for taxon in taxons:
                        if taxon[0:3] == "g__":
                            if ";" in taxon:
                                genus = taxon[:-1]
                            else:
                                genus = taxon
                            self.taxon_genus[TFP.gid]   = genus[3:]
                            self.taxon_string[genus] = TFP.taxonomy
                else:
                    organism = TFP.organism.split()[0]
                    if '"' in organism:
                        organism = organism[1:]
                    self.taxon_genus[TFP.gid] = organism
                    self.taxon_string[organism] = TFP.taxonomy

class GroupFileParser(object):
    def __init__(self,l):
        self.readGroupFile(l)
        
    def readGroupFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num = int(tabs[0].split(':')[0].split()[1])
        self.pid_sqids = tabs[2].split(',')

class GroupData(object):
    def __init__(self,group_data):
        self.group_data         = {}
        self.group_membership   = {}
        self.buildGroupData(group_data)
    
    def buildGroupData(self, group_data):
        with open(group_data, 'r') as fh:
            for l in fh:
                GFP = GroupFileParser(l)
                GFP.readGroupFile(l)
                for pid_sqid in GFP.pid_sqids:
                    try:
                        self.group_data[GFP.group_num] += [pid_sqid]
                    except KeyError:
                        self.group_data[GFP.group_num] = [pid_sqid]
                    
                    # add pid_sqid as members
                    self.group_membership[pid_sqid] = GFP.group_num

class HitDataParser(object):
    def __init__(self,l):
        self.readHitData(l)
        
    def readHitData(self,l):
        tabs = l.rstrip().split("\t")
        self.hid=                        tabs[0]
        self.pid=                        tabs[1]
        self.ani_comp=                   tabs[2]
        self.ident=                      tabs[3]
        self.gid_1=                      tabs[4]
        self.habitat_1=                  tabs[5]
        self.phylum_1=                   tabs[6]
        self.genus_1=                    tabs[7]
        self.status_1=                   tabs[8]
        self.sequencingMethod_1=         tabs[9]
        self.sequencingCentre_1=         tabs[10]
        self.horizontalTransferred_1=    tabs[11]
        self.genomeSize_1=               tabs[12]
        self.scaffoldCount_1=            tabs[13]
        self.len_1=                      tabs[14]
        self.strand_1=                   tabs[15]
        self.cid_1=                      tabs[16]
        self.contig_1=                   tabs[17]
        self.contigLength_1=             tabs[18]
        self.sqid_1=                     tabs[19]
        self.start_1=                    tabs[20]
        self.gid_2=                      tabs[21]
        self.habitat_2=                  tabs[22]
        self.phylum_2=                   tabs[23]
        self.genus_2=                    tabs[24]
        self.status_2=                   tabs[25]
        self.sequencingMethod_2=         tabs[26]
        self.sequencingCentre_2=         tabs[27]
        self.horizontalTransferred_2=    tabs[28]
        self.genomeSize_2=               tabs[29]
        self.scaffoldCount_2=            tabs[30]
        self.len_2=                      tabs[31]
        self.strand_2=                   tabs[32]
        self.cid_2=                      tabs[33]
        self.contig_2=                   tabs[34]
        self.contigLength_2=             tabs[35]
        self.sqid_2=                     tabs[36]
        self.start_2=                    tabs[37]
        self.dirty=                      tabs[38]

class HitData(object):
    def __init__(self, hitdata_file, taxon_file, target_genus1, target_genus2):
        self.TD                  = TaxonomyData(taxon_file)
        self.target_genus1       = target_genus1.lower()
        self.target_genus2       = target_genus2.lower()
        self.target_pidsqids     = {}
        self.buildHitData(hitdata_file)
        
    def buildHitData(self, hitdata_file):
        with open(hitdata_file) as fh:
            for l in fh:
                if l[0:3] != 'hid':
                    HDP = HitDataParser(l)
                    HDP.readHitData(l)
                    
                    # convert from gid to genus, and check if target genus
                    status = self.checkGIDs(HDP.gid_1,
                                            HDP.gid_2)
                    
                    if status:
                        # add to dictionary of target pidsqids
                        pidsqid = "%s_%s" % (HDP.pid,
                                             HDP.sqid_1)
                        self.target_pidsqids[pidsqid] = 1
                
    def checkGIDs(self, gid1, gid2):
        # get updated taxonomy
        genus1 = self.TD.taxon_genus[gid1].lower() 
        genus2 = self.TD.taxon_genus[gid2].lower()

        if genus1 == self.target_genus1: 
            if genus2 == self.target_genus2:
                return True
            
        elif genus1 == self.target_genus2:
            if genus2 == self.target_genus1:
                return True
        else:
            return False

class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnoFile(l)
        
    def readAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqidgene        = tabs[0]
        self.pidsqid            = "_".join(tabs[0].split("_")[0:2])
        self.identity           = tabs[1]
        self.evalue             = tabs[2]
        self.length             = tabs[3]
        self.cog_id             = tabs[5]
        self.cog_description    = tabs[6]
        self.img_annotation     = tabs[8]
        
class AnnotationData(object):
    def __init__(self, group_file, taxon_file, genus1, genus2, hitdata_file):
        self.phage_anno             = {}
        self.GD                     = GroupData(group_file) 
        self.HD                     = HitData(hitdata_file, taxon_file, genus1, genus2)
        self.genus1                 = genus1
        self.genus2                 = genus2
    
    def buildAnnoData(self,anno_file):
        # read in anno file
        with open(anno_file) as fh:
            for l in fh:
                AFP = AnnotationFileParser(l)
                AFP.readAnnoFile(l)
                
                # check if pidsqid in transfer group
                if self.checkPidSqidGroupMembership(AFP.pidsqid):
                    # check if target pidsqid
                    if self.checkPidSqidIsTarget(AFP.pidsqid): 
                        transfer_group = self.GD.group_membership[AFP.pidsqid]
                        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (AFP.pidsqidgene,
                                                                  transfer_group,
                                                                  AFP.cog_id,
                                                                  AFP.cog_description,
                                                                  AFP.identity,
                                                                  AFP.evalue,
                                                                  AFP.length,
                                                                  AFP.img_annotation
                                                                  )
    
    def checkPidSqidIsTarget(self, pidsqid):
        try:
            lenny = self.HD.target_pidsqids[pidsqid]
            return True
        except KeyError:
            return False
                
    def checkPidSqidGroupMembership(self, pidsqid):
        try:
            lenny = self.GD.group_membership[pidsqid]
            return True
        except KeyError:           
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
    # call class
    AD = AnnotationData(args.group_file,
                        args.taxon_file,
                        args.genus1,
                        args.genus2,
                        args.hitdata_file)
    
    # build annotation data
    AD.buildAnnoData(args.anno_file)
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('anno_file', help="File containing tab delimited annotation information")
    parser.add_argument('group_file', help="File containing transfer groups")
    parser.add_argument('taxon_file', help="File containing Phil's improved taxonomy")
    parser.add_argument('hitdata_file', help="File containing hitdata")
    parser.add_argument('-g1','--genus1', help="Please provide genus 1")
    parser.add_argument('-g2','--genus2', help="Please provide genus 2")
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
