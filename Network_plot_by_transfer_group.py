#!/usr/bin/env python
###############################################################################
#
# __Network_plot_by_transfer_group__.py - description!
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
import networkx as nx
import matplotlib.pyplot as plt
#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure
from cb2cols import Cb2Cols as CB2

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnoFile(l)
        
    def readAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqidgene    = tabs[0]
        self.pidsqid        = "_".join(tabs[0].split("_")[0:2])
        self.evalue         = tabs[2]
        self.img_annotation = tabs[8]
        
class AnnotationData(object):
    def __init__(self,anno_file):
        self.phage_anno             = {}
        self.plasmid_anno           = {}
        self.transposon_anno        = {}
        self.antibiotic_anno        = {}
        self.master_list_pidsqid    = {}
        self.buildAnnoData(anno_file)
        
    def generalChecker(self, search_terms, pidsqid, annotation, dictionary):
        for term in search_terms:
            if term.lower() in annotation.lower():
                try:
                    dictionary[pidsqid] += 1 
                except:
                    dictionary[pidsqid] = 1 
                    
    def plasmidChecker(self, pidsqid, annotation):
        plasmid_search_terms    = ["plasmid"]
        self.generalChecker(plasmid_search_terms,
                            pidsqid,
                            annotation,
                            self.plasmid_anno)
        
    def phageChecker(self, pidsqid, annotation):
        phage_search_terms      = ["phage",
                                   "head",
                                   "terminase",
                                   "capsid",
                                   "tail",
                                   "integrase"]
        self.generalChecker(phage_search_terms,
                            pidsqid,
                            annotation,
                            self.phage_anno)
                
        
    def transposonChecker(self, pidsqid, annotation):
        transposon_search_terms = ["transposon"]
        self.generalChecker(transposon_search_terms,
                            pidsqid,
                            annotation,
                            self.transposon_anno)
        
    def antibioticChecker(self, pidsqid, annotation):
        antibiotic_search_terms = ["tetracycline"]
        self.generalChecker(antibiotic_search_terms,
                            pidsqid,
                            annotation,
                            self.antibiotic_anno)
        
    def initialiseDicts(self, pidsqid, dict):
        try:
            lenny = dict[pidsqid]
        except KeyError:
            dict[pidsqid] = 0
        
    def buildAnnoData(self,anno_file):
        # read in anno file
        with open(anno_file) as fh:
            for l in fh:
                AFP = AnnotationFileParser(l)
                AFP.readAnnoFile(l)
                
                # add pidsqid to master list
                self.master_list_pidsqid[AFP.pidsqid] = 1
                
                # initialise dicts
                self.initialiseDicts(AFP.pidsqid, self.plasmid_anno)
                self.initialiseDicts(AFP.pidsqid, self.phage_anno)
                self.initialiseDicts(AFP.pidsqid, self.transposon_anno)
                self.initialiseDicts(AFP.pidsqid, self.antibiotic_anno)
                
                # check annotations
                self.plasmidChecker(AFP.pidsqid, AFP.img_annotation)
                self.phageChecker(AFP.pidsqid, AFP.img_annotation)
                self.transposonChecker(AFP.pidsqid, AFP.img_annotation)
                self.antibioticChecker(AFP.pidsqid, AFP.img_annotation)

class GroupFileParser(object):
    def __init__(self,l):
        self.readGroupFile(l)
        
    def readGroupFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num = int(tabs[0].split(':')[0].split()[1])
        self.pid_sqids = tabs[2].split(',')
        
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
        self.bodysite_1=                 tabs[5]
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
        self.bodysite_2=                 tabs[22]
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

class TransferGroupPropertiesParser(object):
    def __init__(self, l):
        self.readTransGroupPropFile(l)
        
    def readTransGroupPropFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num      = int(tabs[0].split("_")[1])
        self.gid_count      = int(tabs[1])
        self.habitat_count  = int(tabs[2])
        self.phylum_count   = int(tabs[3])

class VisualiseTransferGroups(object):
    def __init__(self, anno_file):
        self.group_data             = {} # group_num = [pid_sqid1, pid_sqid2,...]
        self.group_membership       = {}
        self.hit_data               = {} 
        self.pidsqid_lookup         = {} # pid_sqid = [genus1, genus2]
        self.genus_list             = {}
        self.phylum_data            = {}  
        self.body_site_gid          = {} 
        self.body_site              = {}
        self.phylum_counts          = {}
        self.transfer_group_props   = {}  
        self.dodgey_transfer_groups = [1,14,40] 
        self.habitat_data           = {}
        self.anno_data_plasmid      = {}
        self.anno_data_phage        = {}
        self.anno_data_transposon   = {}
        self.anno_data_antibiotic   = {}
        self.AD                     = AnnotationData(anno_file)
        
    def transferGroupProperties(self, transfer_group_file):
        with open(transfer_group_file, 'r') as fh:
            for l in fh:
                TGPP = TransferGroupPropertiesParser(l)
                TGPP.readTransGroupPropFile(l)
                self.transfer_group_props[TGPP.group_num] = [TGPP.habitat_count, TGPP.phylum_count]
        
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
                        
    def addGroupToHitData(self, genus1, genus2, group_num):
        try:
            self.hit_data[genus1][genus2][group_num] = 1
        except KeyError: 
            try: 
                self.hit_data[genus1][genus2] = {group_num:1}
            except KeyError:
                self.hit_data[genus1] = {genus2:{group_num:1}}
                    
    def addToGenusList(self, genus1, gid1):
        try:
            self.genus_list[genus1][gid1] = 1 
        except KeyError:
            self.genus_list[genus1] = {gid1:1}
             
    def checkHabitatHomogeneity(self, genus, habitat):
        try:
            if habitat != self.habitat_data[genus]:
                self.habitat_data[genus] = 'mixed'
        except KeyError:
            self.habitat_data[genus] = habitat
            
    def addAnnotation(self, master_anno_dict, slave_anno_dict, pidsqid1, pidsqid2, genus1, genus2):
        try:
            # check pidsqid is in master annotation list
            lenny = self.AD.master_list_pidsqid[pidsqid1]
            
            # convert to genus -> genus
            if master_anno_dict[pidsqid1] > 0: 
                try: 
                    slave_anno_dict[genus1][genus2] = 1 
                except KeyError:
                    slave_anno_dict[genus1] = {genus2:1}
        except KeyError: 
            pass
            
        try:
            # check pidsqid is in master annotation list 
            carl = self.AD.master_list_pidsqid[pidsqid2]
            
            # convert to genus -> genus
            if master_anno_dict[pidsqid2] > 0: 
                try: 
                    slave_anno_dict[genus2][genus1] = 1 
                except KeyError:
                    slave_anno_dict[genus2] = {genus1:1}
        except KeyError:
            pass
        
    def buildPidSqidLookUp(self, pidsqid, genus1, genus2, dict):
        try:
            dict[(genus1,genus2)] += [pidsqid] 
        except KeyError:
            dict[(genus1,genus2)] = [pidsqid]
            
    def buildHitDataAll(self, hit_data, phylum_thresh, habitat_thresh, dodgey_groups):

        with open(hit_data, 'r') as fh:
            for l in fh:
                if l[0:3] != 'hid':
                    HDP = HitDataParser(l)
                    HDP.readHitData(l)
                    
                    # grab pid_sqid 
                    pid_sqid1 = "%s_%s" % (HDP.pid, HDP.sqid_1)
                    pid_sqid2 = "%s_%s" % (HDP.pid, HDP.sqid_2)
                    
                    # build pidsqid lookup
                    self.buildPidSqidLookUp(pid_sqid1, HDP.genus_1, HDP.genus_2, self.pidsqid_lookup)
                    self.buildPidSqidLookUp(pid_sqid2, HDP.genus_2, HDP.genus_1, self.pidsqid_lookup)
                    
                    # check group membership
                    try:
                        group_member  = self.group_membership[pid_sqid1]
                        if phylum_thresh:
                            if self.transfer_group_props[group_member][1] <= phylum_thresh:
                                
                                # build hit data by group
                                self.addGroupToHitData(HDP.genus_1,HDP.genus_2,group_member)
                                
                                # add gids to genus list
                                self.addToGenusList(HDP.genus_1, HDP.gid_1)
                                self.addToGenusList(HDP.genus_2, HDP.gid_2)

                                # build phylum dictionary
                                self.phylum_data[HDP.genus_1] = HDP.phylum_1
                                self.phylum_data[HDP.genus_2] = HDP.phylum_2
                                
                                # build habitat dictionary
                                self.checkHabitatHomogeneity(HDP.genus_1, HDP.bodysite_1)
                                self.checkHabitatHomogeneity(HDP.genus_2, HDP.bodysite_2)
                                
                        elif habitat_thresh:
                            if self.transfer_group_props[group_member][0] <= habitat_thresh:
                                
                                # build hit data by group
                                self.addGroupToHitData(HDP.genus_1,HDP.genus_2,group_member)
                                
                                # add gids to genus list
                                self.addToGenusList(HDP.genus_1, HDP.gid_1)
                                self.addToGenusList(HDP.genus_2, HDP.gid_2)

                                # build phylum dictionary
                                self.phylum_data[HDP.genus_1] = HDP.phylum_1
                                self.phylum_data[HDP.genus_2] = HDP.phylum_2
                        
                                # build habitat dictionary
                                self.checkHabitatHomogeneity(HDP.genus_1, HDP.bodysite_1)
                                self.checkHabitatHomogeneity(HDP.genus_2, HDP.bodysite_2)
                       
                        elif dodgey_groups:
                            if group_member not in self.dodgey_transfer_groups:
                                # build hit data by group
                                self.addGroupToHitData(HDP.genus_1,HDP.genus_2,group_member)
                                
                                # add gids to genus list
                                self.addToGenusList(HDP.genus_1, HDP.gid_1)
                                self.addToGenusList(HDP.genus_2, HDP.gid_2)

                                # build phylum dictionary
                                self.phylum_data[HDP.genus_1] = HDP.phylum_1
                                self.phylum_data[HDP.genus_2] = HDP.phylum_2
                                
                                # build habitat dictionary
                                self.checkHabitatHomogeneity(HDP.genus_1, HDP.bodysite_1)
                                self.checkHabitatHomogeneity(HDP.genus_2, HDP.bodysite_2)
                        else: 
                            # build hit data by group
                            self.addGroupToHitData(HDP.genus_1,HDP.genus_2,group_member)
                            
                            # add gids to genus list
                            self.addToGenusList(HDP.genus_1, HDP.gid_1)
                            self.addToGenusList(HDP.genus_2, HDP.gid_2)

                            # build phylum dictionary
                            self.phylum_data[HDP.genus_1] = HDP.phylum_1
                            self.phylum_data[HDP.genus_2] = HDP.phylum_2
                        
                            # build habitat dictionary
                            self.checkHabitatHomogeneity(HDP.genus_1, HDP.bodysite_1)
                            self.checkHabitatHomogeneity(HDP.genus_2, HDP.bodysite_2)
                            
                    except KeyError:
                        pass
            
    def buildHitDataByGroup(self, hit_data, group_number):
        # reset data 
        self.hit_data       = {}
        self.genus_list     = {}
        self.phylum_data    = {}
        self.body_site      = {}
        self.body_site_gid  = {}
        self.phylum_counts  = {}
        
        with open(hit_data, 'r') as fh:
            for l in fh: 
                if l[0:3] != 'hid': 
                    HDP = HitDataParser(l)
                    HDP.readHitData(l)
                    pid_sqid1 = "%s_%s" % (HDP.pid, HDP.sqid_1)
                    pid_sqid2 = "%s_%s" % (HDP.pid, HDP.sqid_2)
                    
                    # only add to self.hit_data if in user-specific group
                    if pid_sqid1 in self.group_data[group_number] or pid_sqid2 in self.group_data[group_number]:
                        try:
                            self.hit_data[HDP.genus_1][HDP.genus_2] += 1
                        except KeyError:
                            try:
                                self.hit_data[HDP.genus_1][HDP.genus_2] = 1
                            except KeyError:
                                self.hit_data[HDP.genus_1] = {HDP.genus_2:1}
                        
                        # add genus to genus node list 
                        try:
                            self.genus_list[HDP.genus_1][HDP.gid_1] = 1
                        except KeyError:
                            self.genus_list[HDP.genus_1] = {HDP.gid_1:1}
                        try:
                            self.genus_list[HDP.genus_2][HDP.gid_2] = 1
                        except KeyError:
                            self.genus_list[HDP.genus_2] = {HDP.gid_2:1}
                        
                        # build phylum dictionary
                        self.phylum_data[HDP.genus_1] = HDP.phylum_1
                        self.phylum_data[HDP.genus_2] = HDP.phylum_2
    
                        # count bodysites
                        # check for uniqueness
                        try:
                            self.body_site_gid[HDP.gid_1] += 1 
                        except KeyError:
                            self.body_site_gid[HDP.gid_1] = 1 
                            # add count to bodysite
                            try:
                                self.body_site[HDP.bodysite_1] += 1 
                            except KeyError:
                                self.body_site[HDP.bodysite_1] = 1 
                            
                            # add count to phylum
                            try:
                                self.phylum_counts[HDP.phylum_1] += 1 
                            except KeyError:
                                self.phylum_counts[HDP.phylum_1] = 1
                        
                        try:
                            self.body_site_gid[HDP.gid_2] += 1 
                        except KeyError:
                            self.body_site_gid[HDP.gid_2] = 1 
                            # add count to bodysite
                            try:
                                self.body_site[HDP.bodysite_2] += 1 
                            except KeyError:
                                self.body_site[HDP.bodysite_2] = 1
                                
                            # add count to phylum
                            try:
                                self.phylum_counts[HDP.phylum_2] += 1 
                            except KeyError:
                                self.phylum_counts[HDP.phylum_2] = 1 

        ##############################
        print "Group_%d\t%d\t%d\t%d" % (group_number,len(self.body_site_gid),len(self.body_site),len(self.phylum_counts))
        #for body_site in self.body_site.keys():
        #    print "%s : %d" % (body_site, self.body_site[body_site])                    
 
    def buildNetworkPlotAll(self, outfile):
        """Create network graph"""
        # objects 
        edgewidth                       = []
        genus_node_size_values          = []
        phylum_colour_values            = []
        phylumCols                      = {}
        phylumDict                      = {}
        genus_node_size                 = {}
        habitatCols                     = {}
        habitatDict                     = {}
        habitat_colour_values           = []
        edgeColsPlasmid                 = {}
        edgeColsPhage                   = {}
        edgeColsTransposon              = {}
        edgeColsAntibiotic              = {}
        edge_colour_values_plasmid      = {}
        edge_colour_values_phage        = {}
        edge_colour_values_transposon   = {}
        edge_colour_values_antibiotic   = {}
        
        
        # Get ColorBrewer colous
        cb2                     = CB2()
        col_set                 = "qualSet1"
        ColBrewColours          = cb2.maps[col_set].values()[0:10]
        
        # initialise graph
        G=nx.Graph()
        
        # create list of genus nodes
        genus_nodes = self.genus_list.keys()
        
        # add nodes 
        G.add_nodes_from(genus_nodes)
        
        # build node properties
        for genus in genus_nodes:
            
            ##############################
            # Size nodes by no. of genomes
            ##############################
            genus_node_size[genus] = len(self.genus_list[genus])*100
            
            ##############################
            # Colour nodes by phylum
            ##############################
            try:
                phylumCols[genus] = phylumDict[self.phylum_data[genus]]
            except KeyError:
                if self.phylum_data[genus] == "Bacteroidetes":
                    phylumDict[self.phylum_data[genus]] =  ColBrewColours[0] #"#ea00ff"
                    phylumCols[genus] = ColBrewColours[0] #"#ea00ff"

                elif self.phylum_data[genus] == "Lentisphaerae":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[1] #"#ffb700"
                    phylumCols[genus] = ColBrewColours[1] #"#ffb700"
                    
                elif self.phylum_data[genus] == "Firmicutes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[2] #"#0047ff"
                    phylumCols[genus] = ColBrewColours[2] #"#0047ff"

                elif self.phylum_data[genus] == "Spirochaetes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[3] #"#14ff00"
                    phylumCols[genus] = ColBrewColours[3] #"#14ff00"

                elif self.phylum_data[genus] == "Synergistetes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[4] #"#6600CC"
                    phylumCols[genus] = ColBrewColours[4] #"#6600CC"

                elif self.phylum_data[genus] == "Actinobacteria":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[5] #"#ffff00"
                    phylumCols[genus] = ColBrewColours[5] #"#ffff00"

                elif self.phylum_data[genus] == "Tenericutes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[6] #"#006600"
                    phylumCols[genus] = ColBrewColours[6] #"#0080ff"

                elif self.phylum_data[genus] == "Fusobacteria":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[7] #"#00e0ff"
                    phylumCols[genus] = ColBrewColours[7] #"#00e0ff"

                elif self.phylum_data[genus] == "Proteobacteria":
                    phylumDict[self.phylum_data[genus]] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                    phylumCols[genus] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
            
            ##############################
            # Colour nodes by habitat
            ##############################       
            """
            Eye,Airways,internal_organs,Gastrointestinal tract,Blood,skin,Urogenital tract,Ear,Oral
            """
            try:
                habitatCols[genus] = habitatDict[self.phylum_data[genus]]
            except KeyError:
                if self.habitat_data[genus] == "Gastrointestinal tract":
                    habitatDict[self.habitat_data[genus]] =  ColBrewColours[0] #"#ea00ff"
                    habitatCols[genus] = ColBrewColours[0] #"#ea00ff"

                elif self.habitat_data[genus] == "Eye":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[1] #"#ffb700"
                    habitatCols[genus] = ColBrewColours[1] #"#ffb700"
                    
                elif self.habitat_data[genus] == "Oral":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[2] #"#0047ff"
                    habitatCols[genus] = ColBrewColours[2] #"#0047ff"

                elif self.habitat_data[genus] == "Airways":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[3] #"#14ff00"
                    habitatCols[genus] = ColBrewColours[3] #"#14ff00"

                elif self.habitat_data[genus] == "internal_organs":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[4] #"#6600CC"
                    habitatCols[genus] = ColBrewColours[4] #"#6600CC"

                elif self.habitat_data[genus] == "Blood":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[5] #"#ffff00"
                    habitatCols[genus] = ColBrewColours[5] #"#ffff00"

                elif self.habitat_data[genus] == "skin":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[6] #"#006600"
                    habitatCols[genus] = ColBrewColours[6] #"#0080ff"

                elif self.habitat_data[genus] == "Urogenital tract":
                    habitatDict[self.habitat_data[genus]] = ColBrewColours[7] #"#00e0ff"
                    habitatCols[genus] = ColBrewColours[7] #"#00e0ff"

                elif self.habitat_data[genus] == "Ear":
                    habitatDict[self.habitat_data[genus]] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                    habitatCols[genus] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                
                elif self.habitat_data[genus] == "mixed":
                    habitatDict[self.habitat_data[genus]] = '#C9C9C9' #ColBrewColours[8] #"#ff1e00"
                    habitatCols[genus] = '#C9C9C9' #ColBrewColours[8] #"#ff1e00"
            
        # build node data arrays for plotting
        genus_node_size_values = [genus_node_size.get(node) for node in G.nodes()]
        phylum_colour_values   = [phylumCols.get(node) for node in G.nodes()]
        habitat_colour_values   = [habitatCols.get(node) for node in G.nodes()]
        
        # add edges
        for genus1 in self.hit_data.keys():
            for genus2 in self.hit_data[genus1]:
                
                ##############################
                # Colour edges by annotation
                ##############################
                                
                # add edges between two nodes
                G.add_edge(genus1,
                           genus2,
                           capacity = len(self.hit_data[genus1][genus2]))
                
                try:
                    for pidsqid in self.pidsqid_lookup[(genus1,genus2)]:
                        if self.AD.antibiotic_anno[pidsqid] >0: 
                            # contains phage annotation
                            edgeColsAntibiotic[(genus1,genus2)] = "#FF6969"
                            edgeColsAntibiotic[(genus2,genus1)] = "#FF6969"
                        else:
                            edgeColsAntibiotic[(genus1,genus2)] = "#E1E1E1"
                            edgeColsAntibiotic[(genus2,genus1)] = "#E1E1E1"
                except KeyError:
                    edgeColsAntibiotic[(genus1,genus2)] = "#E1E1E1"
                    edgeColsAntibiotic[(genus2,genus1)] = "#E1E1E1"
                
                
        # build edge data arrays for plotting
        
        edge_colour_values_antibiotic = [edgeColsAntibiotic.get(edge) for edge in G.edges()]
        
        # Edgewidth is # of transfer groups
        for (u,v,d) in G.edges(data=True):
            edgewidth.append(int(str(d).split(':')[-1].split()[0].split('}')[0]))
        
        # Build network plot
        pos= nx.spring_layout(G,k=0.30,iterations=500)
        #pos= nx.circular_layout(G)
        fig = plt.figure(figsize=(30,15),dpi=300)
        
        #plt.subplot(1,1,1,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               node_size= genus_node_size_values,
                               node_color = phylum_colour_values, 
                               alpha=0.7,
                               with_labels=True)
        #node_color = habitat_colour_values,
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = edge_colour_values_antibiotic,
                               width=edgewidth)
        #edge_color = edge_colour_values_phage,
        #edge_color = "#E1E1E1",
        nx.draw_networkx_labels(G,
                                pos,
                                font_size=8)
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')

        plt.savefig("%s" % (outfile),format='png')
    
    def buildNetworkPlotByTransferGroup(self, group_number):
        """Create network graph"""
        G=nx.Graph()
    
        # objects
        genus_nodes             = self.genus_list.keys()
        genus_node_size         = {}
        phylumDict              = {}
        phylumCols              = {}
        edgewidth               = []
        genus_node_size_values  = []
        phylum_colour_values    = []
        # Get ColorBrewer colous
        cb2                     = CB2()
        col_set                 = "qualSet1"
        ColBrewColours          = cb2.maps[col_set].values()[0:10]
        
        # add genus (nodes)
        G.add_nodes_from(genus_nodes)  
        
        # size nodes by number of genomes per genus
        for genus in genus_nodes:
            
            ##############################
            # Size nodes by no. of genomes
            ##############################
            genus_node_size[genus] = len(self.genus_list[genus])*100
            
            ##############################
            # Colour nodes by phylum
            ##############################
            try:
                phylumCols[genus] = phylumDict[self.phylum_data[genus]]
            except KeyError:
                if self.phylum_data[genus] == "Bacteroidetes":
                        phylumDict[self.phylum_data[genus]] =  ColBrewColours[0] #"#ea00ff"
                        phylumCols[genus] = ColBrewColours[0] #"#ea00ff"

                elif self.phylum_data[genus] == "Lentisphaerae":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[1] #"#ffb700"
                    phylumCols[genus] = ColBrewColours[1] #"#ffb700"
                    
                elif self.phylum_data[genus] == "Firmicutes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[2] #"#0047ff"
                    phylumCols[genus] = ColBrewColours[2] #"#0047ff"

                elif self.phylum_data[genus] == "Spirochaetes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[3] #"#14ff00"
                    phylumCols[genus] = ColBrewColours[3] #"#14ff00"

                elif self.phylum_data[genus] == "Synergistetes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[4] #"#6600CC"
                    phylumCols[genus] = ColBrewColours[4] #"#6600CC"

                elif self.phylum_data[genus] == "Actinobacteria":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[5] #"#ffff00"
                    phylumCols[genus] = ColBrewColours[5] #"#ffff00"

                elif self.phylum_data[genus] == "Tenericutes":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[6] #"#006600"
                    phylumCols[genus] = ColBrewColours[6] #"#0080ff"

                elif self.phylum_data[genus] == "Fusobacteria":
                    phylumDict[self.phylum_data[genus]] = ColBrewColours[7] #"#00e0ff"
                    phylumCols[genus] = ColBrewColours[7] #"#00e0ff"

                elif self.phylum_data[genus] == "Proteobacteria":
                    phylumDict[self.phylum_data[genus]] = '#00CCCC' #ColBrewColours[8] #"#ff1e00"
                    phylumCols[genus] = '#00CCCC' #ColBrewColours[8] #"#ff1e00" 
            
        
        # build data arrays for plotting
        genus_node_size_values = [genus_node_size.get(node) for node in G.nodes()]    
        phylum_colour_values   = [phylumCols.get(node) for node in G.nodes()]
        
        # add edges
        for genus1 in self.hit_data.keys():
            for genus2 in self.hit_data[genus1]:
                G.add_edge(genus1,
                           genus2,
                           capacity = self.hit_data[genus1][genus2])
                
        for (u,v,d) in G.edges(data=True):
            edgewidth.append(int(str(d).split(':')[-1].split()[0].split('}')[0]))
        
        # Build network plot
        pos= nx.spring_layout(G,iterations=500)
        fig = plt.figure(figsize=(21,10),dpi=300)
        
        
        #plt.subplot(1,1,1,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,
                               pos,
                               linewidths=0,
                               alpha=0.7,
                               node_size=genus_node_size_values,
                               node_color = phylum_colour_values,
                               with_labels=True)
        nx.draw_networkx_edges(G,
                               pos,
                               edge_color = "#E1E1E1",
                               width=edgewidth)
        nx.draw_networkx_labels(G,
                                pos,
                                font_size=10)
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')
        #plt.show()
        plt.savefig("Group_%d.png" % (group_number),format='png')
        
    def makeNetworkForEachGroup(self,hit_data, print_only):
        """Output individual files for group"""
        for group in self.group_data.keys():
            
            self.buildHitDataByGroup(hit_data, group)
            
            if print_only:
                pass
            #elif type.lower() == 'group':
            #    self.buildNetworkPlotByTransferGroup(group)
            #elif type.lower() == 'all': 
            #    pass
        

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
    VTG = VisualiseTransferGroups(args.anno_file)
    #AD  = AnnotationData() 
    
    # build annotation data
    #AD.buildAnnoData(args.anno_file)
    
    # grab group data
    VTG.buildGroupData(args.group_data)
    
    # build group properties
    VTG.transferGroupProperties(args.trans_group_props)
    
    # build hit data ALL
    VTG.buildHitDataAll(args.hit_data, args.phylum_thresh, args.habitat_thresh, args.dodgey_groups)
    
    # build network plot
    #VTG.buildNetworkPlotAll(args.out_file)
    
    # build network plots for each group
    VTG.makeNetworkForEachGroup(args.hit_data,args.print_only)
    
    # build group-hit data
    #VTG.buildHitDataByGroup(args.hit_data, args.group_number)
    
    # build network plot
    #VTG.buildNetworkPlot()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hit_data', help="File containing hit data")
    parser.add_argument('group_data', help="File containing transfer groups")
    parser.add_argument('trans_group_props', help="File containing phylum and habitat counts")
    parser.add_argument('anno_file', help="File containing tab delimited annotation information")
    parser.add_argument('-p','--print_only', default=False,help="Print outs information, but does not create network plots")
    parser.add_argument('-pt','--phylum_thresh', type=int,default=False,help="Set group phylum member limit")
    parser.add_argument('-ht','--habitat_thresh', type=int,default=False,help="Set group habitat member limit")
    parser.add_argument('-dg','--dodgey_groups', default=False,help="Remove dodgey groups from analysis")
    parser.add_argument('-o','--out_file', default='Network_plot.png',help="Set output file, and directory.")
    #parser.add_argument('group_number', type=int,help="User-specific group number")
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
