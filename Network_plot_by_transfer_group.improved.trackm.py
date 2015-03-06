#!/usr/bin/env python
###############################################################################
#
# __network_plot_by_transfer_groups__.py - description!
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

#global imports
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
np.seterr(all='raise')

# local imports
from cb2cols import Cb2Cols as CB2

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
                    
    def generalCheckerPhage(self, search_terms, pidsqid, annotation, dictionary):
        for term in search_terms:
            if term.lower() in annotation.lower():
                try:
                    dictionary[pidsqid][term] = 1  
                except:
                    dictionary[pidsqid] = {term:1} 
                    
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
                                   "integrase",
                                   "portal",
                                   "spike",
                                   "virion formation",
                                   "coat"]
        
        self.generalCheckerPhage(phage_search_terms,
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
        self.pidsqid=                    tabs[0]
        self.transfer_group=             tabs[1]
        self.hid=                        tabs[2]
        self.pid=                        tabs[3]
        self.ani_comp=                   tabs[4]
        self.ident=                      tabs[5]
        self.gid_1=                      tabs[6]
        self.habitat_1=                  tabs[7]
        self.phylum_1=                   tabs[8]
        self.genus_1=                    tabs[9]
        self.status_1=                   tabs[10]                       
        self.sequencingMethod_1=         tabs[11]
        self.sequencingCentre_1=         tabs[12]
        self.horizontalTransferred_1=    tabs[13]
        self.genomeSize_1=               tabs[14]
        self.scaffoldCount_1=            tabs[15]
        self.len_1=                      tabs[16]
        self.strand_1=                   tabs[17]
        self.cid_1=                      tabs[18]
        self.contig_1=                   tabs[19]
        self.contigLength_1=             int(tabs[20])
        self.sqid_1=                     tabs[21]
        self.start_1=                    tabs[22]
        self.gid_2=                      tabs[23]
        self.habitat_2=                  tabs[24]
        self.phylum_2=                   tabs[25]
        self.genus_2=                    tabs[26]
        self.status_2=                   tabs[27]
        self.sequencingMethod_2=         tabs[28]
        self.sequencingCentre_2=         tabs[29]
        self.horizontalTransferred_2=    tabs[30]
        self.genomeSize_2=               tabs[31]
        self.scaffoldCount_2=            tabs[32]
        self.len_2=                      tabs[33]
        self.strand_2=                   tabs[34]
        self.cid_2=                      tabs[35]
        self.contig_2=                   tabs[36]
        self.contigLength_2=             int(tabs[37])
        self.sqid_2=                     tabs[38]
        self.start_2=                    tabs[39]

class TransferGroupPropertiesParser(object):
    def __init__(self, l):
        self.readTransGroupPropFile(l)
        
    def readTransGroupPropFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num      = int(tabs[0].split("_")[1])
        self.gid_count      = int(tabs[1])
        self.habitat_count  = int(tabs[2])
        self.phylum_count   = int(tabs[3])

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
        
class KmerFileParser(object):
    def __init__(self,l):
        self.readKmerFile(l)
        
    def readKmerFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqid    = tabs[0] 
        self.gid1       = tabs[1]
        self.gid2       = tabs[2]
        self.kmerscore  = float(tabs[3])
        self.dgg_mean   = float(tabs[4])
        self.dg1        = float(tabs[5])
        self.dg2        = float(tabs[6])
    
class KmerData(object):
    def __init__(self, kmer_file):
        self.kmer_score         = {}
        self.mean_dgg           = {}
        self.kmer_data          = {} 
        self.buildKmerData(kmer_file)
    
    def buildKmerData(self, kmer_file):
        with open(kmer_file) as fh:
            for l in fh:
                if l[0] != 'p':
                    KFP = KmerFileParser(l)
                    try:
                        self.kmer_data[KFP.gid1][KFP.gid2] += [KFP.kmerscore]
                    except KeyError: 
                        try:
                            self.kmer_data[KFP.gid1][KFP.gid2] = [KFP.kmerscore]
                        except KeyError:
                            self.kmer_data[KFP.gid1] = {KFP.gid2: [KFP.kmerscore]}
                    self.kmer_score[KFP.pidsqid]    = KFP.kmerscore
                    self.mean_dgg[KFP.pidsqid]      = KFP.dgg_mean
            
class HitData(object):
    def __init__(self,group_file, anno_file, taxon_file, kmer_file):
        self.GD                             = GroupData(group_file)
        self.AD                             = AnnotationData(anno_file)
        self.TD                             = TaxonomyData(taxon_file)
        self.KD                             = KmerData(kmer_file)
        self.pidsqid_lookup                 = {}      
        self.hit_data                       = {}
        self.genus_list                     = {}
        self.habitat_homog_data             = {}
        self.habitat_data                   = {}
        self.transfer_group_props           = {}
        self.phylum_data                    = {}
        self.inter_dodgey_transfer_groups   = [2,4,7,14,16,48,1,3,20,26,65] # lenient # [2,4,7,14,16,32,48,1,3,9,10,17,20,26,55,65] # strict # blast scoree 160 cutoff [3, 7, 38] # old taxonomy[1,14,40]  
        self.all_dodgey_transfer_groups     = [168,8,5]
        self.genus_by_habitat               = {}
        self.kmer_directions                = {}
        self.genus_flux                     = {}
        self.contig_length                  = {}
        self.tgs_containing_ecoli           = {}
        self.genus_data                     = {}
        
    def transferGroupProperties(self, transfer_group_file):
        with open(transfer_group_file, 'r') as fh:
            for l in fh:
                TGPP = TransferGroupPropertiesParser(l)
                TGPP.readTransGroupPropFile(l)
                self.transfer_group_props[TGPP.group_num] = [TGPP.habitat_count, TGPP.phylum_count]
    
    def checkGroupMembership(self,pidsqid):
        try:
            group_member = self.GD.group_membership[pidsqid]
            return group_member
        except KeyError:
            return False
        
    def buildPidSqidLookUp(self, pidsqid, gid1, gid2, dict):
        genus1 = self.TD.taxon_genus[gid1]
        genus2 = self.TD.taxon_genus[gid2]
        try:
            dict[(genus1,genus2)] += [pidsqid] 
        except KeyError:
            dict[(genus1,genus2)] = [pidsqid]
            
    def buildHitData(self, gid1, gid2, habitat1, habitat2, transfer_group, phylum_thresh, habitat_thresh, inter_dodgey_thresh, all_dodgey_thresh, status1, status2, genome_status_thresh, ecoli_removed):
        if phylum_thresh:
            if self.transfer_group_props[transfer_group][1] <= phylum_thresh:
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
        
        elif habitat_thresh:
            if self.transfer_group_props[transfer_group][0] <= habitat_thresh:
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
        
        elif inter_dodgey_thresh:
            if transfer_group not in self.inter_dodgey_transfer_groups:
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
        
        elif all_dodgey_thresh:
            if transfer_group not in self.all_dodgey_transfer_groups:
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
        
        elif genome_status_thresh:
            if self.checkGenomeStatus(status1, status2, genome_status_thresh):
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
            
        elif ecoli_removed:
            if self.removeEcoli(gid1, gid2, ecoli_removed, transfer_group):
                self.buildHitDataWrapper(gid1,
                                         gid2,
                                         habitat1,
                                         habitat2,
                                         transfer_group)
        
        else:
            self.buildHitDataWrapper(gid1,
                                     gid2,
                                     habitat1,
                                     habitat2,
                                     transfer_group)
               
    def removeEcoli(self, gid1, gid2, type, transfer_group):
        if type.lower() == 'pidsqid':
            if self.genus_data[gid1] == 'Escherichia' or self.genus_data[gid2] == 'Escherichia':
                return False
            else:
                return True
        elif type.lower() == 'transfer_group':
            try:
                lenny = self.tgs_containing_ecoli[transfer_group]
                print 'Ecoli containing transfer group %s' % (lenny)
                return False
            except KeyError:
                # TG does not contain E.coli
                return True
    
    def checkGenomeStatus(self, status1, status2, user_set_status):
        if user_set_status.lower() == 'all':
            return True
        elif user_set_status.lower() in status1.lower() and user_set_status.lower() in status2.lower():
            return True
    
    def buildKmerDirectionDict(self, gid1, gid2, genus1, genus2, group_member):
        kmer_scores = None
        
        try:
            kmer_scores = self.KD.kmer_data[gid1][gid2]
            
            try:
                self.kmer_directions[genus1][genus2][group_member] = kmer_scores
            except KeyError:
                try:
                    self.kmer_directions[genus1][genus2] = {group_member:kmer_scores}
                except KeyError:
                    self.kmer_directions[genus1] = {genus2:{group_member:kmer_scores}}
        except KeyError:
            kmer_scores = self.KD.kmer_data[gid2][gid1]
            corr_kmer_scores = []
            # reverse kmer score i.e. 1 - 0.7 -> 0.3. # 
            # loop through kmer scores and minus 1
            for kmer_score in kmer_scores:
                corr_kmer_scores.append(1-kmer_score)
            
            try:
                self.kmer_directions[genus1][genus2][group_member] = corr_kmer_scores
            except KeyError:
                try:
                    self.kmer_directions[genus1][genus2] = {group_member:corr_kmer_scores}
                except KeyError:
                    self.kmer_directions[genus1] = {genus2:{group_member:corr_kmer_scores}}
            
    def calculateFlux(self, genus1,genus2,group_member):
        
        # loop through hitdata
        for genus1 in self.hit_data.keys():
            for genus2 in self.hit_data[genus1]:
                # tmp array
                kmer_scores = []
                for group in self.hit_data[genus1][genus2]:
                    kmer_scores.append(self.kmer_directions[genus1][genus2][group])
                # calculate average kmer score
                total_sum = 0
                total_len = 0
                for array in kmer_scores:
                    total_sum += sum(array)
                    total_len += len(array)
                kmer_average = total_sum/float(total_len)
                # add to genus_flux
                try:
                    self.genus_flux[genus1][genus2] = kmer_average
                except KeyError:
                    self.genus_flux[genus1]= {genus2:kmer_average}
    
    def buildHitDataWrapper(self, gid1, gid2, habitat1, habitat2, transfer_group):
        
        # call new taxonomy data from Taxononmy object
        genus1  = self.TD.taxon_genus[gid1] 
        genus2  = self.TD.taxon_genus[gid2]
        phylum1 = self.TD.taxon_phylum[gid1]
        phylum2 = self.TD.taxon_phylum[gid2]
        
        # build kmer direction dictionary
        self.buildKmerDirectionDict(gid1, gid2, genus1, genus2, transfer_group)
        
        # build hit data by group
        self.addGroupToHitData(genus1,genus2,transfer_group)
        
        # calculate flux between genus
        self.calculateFlux(genus1,genus2,transfer_group)
        
        # add gids to genus list
        self.addToGenusList(genus1, gid1)
        self.addToGenusList(genus2, gid2)

        # build phylum dictionary
        self.phylum_data[genus1] = phylum1
        self.phylum_data[genus2] = phylum2
        
        # check habitat homogeneity for each genus
        self.checkHabitatHomogeneity(genus1, habitat1)
        self.checkHabitatHomogeneity(genus2, habitat2)
    
    def checkContigLength(self, contig_len_1, contig_len_2, contig_len_thresh):
        if contig_len_1 > contig_len_thresh and contig_len_2 > contig_len_thresh:
            return True
        else: 
            return False
            
    def addGroupToHitData(self, genus1, genus2, transfer_group):
        kmer_score = self.kmer_directions[genus1][genus2][transfer_group]
        try:
            self.hit_data[genus1][genus2][transfer_group] = kmer_score
        except KeyError: 
            try: 
                self.hit_data[genus1][genus2] = {transfer_group:kmer_score}
            except KeyError:
                self.hit_data[genus1] = {genus2:{transfer_group:kmer_score}}
                    
    def addToGenusList(self, genus1, gid1):
        try:
            self.genus_list[genus1][gid1] = 1 
        except KeyError:
            self.genus_list[genus1] = {gid1:1}
             
    def checkHabitatHomogeneity(self, genus, habitat):
        try:
            if habitat != self.habitat_data[genus]:
                self.habitat_homog_data[genus] = 'mixed'
        except KeyError:
            self.habitat_homog_data[genus] = habitat

    def groupByHabitat(self, group_by_habitat, gid1, gid2):
        """only works with single habitats currently"""
        if self.habitat_data[gid1].lower() == group_by_habitat.lower() and self.habitat_data[gid2].lower() == group_by_habitat.lower():
            
            # add genus to habitat dict
            genus1  = self.TD.taxon_genus[gid1] 
            genus2  = self.TD.taxon_genus[gid2] 
            self.genus_by_habitat[genus1] = 1
            self.genus_by_habitat[genus2] = 1
            
            return True
        
    def getTGsContainingEcoli(self, genus1, genus2, transfer_group):
        if genus1 == 'Escherichia' or genus2 == 'Escherichia': 
            self.tgs_containing_ecoli[transfer_group] = 1
    
    def mainWrapper(self, hit_file, transfer_group_props, phylum_thresh, habitat_thresh, inter_dodgey_thresh, all_dodgey_thresh, group_by_habitat, contig_len_thresh, genome_status_thresh, remove_ecoli):
        with open(hit_file) as fh:
            for l in fh: 
                if l[0:3] != 'pid':
                    # parse hit data file
                    HDP = HitDataParser(l)
                    HDP.readHitData(l)
                    
                    # grab pid_sqid 
                    pid_sqid1 = "%s_%s" % (HDP.pid, HDP.sqid_1)
                    pid_sqid2 = "%s_%s" % (HDP.pid, HDP.sqid_2)
                    
                    # genus data
                    self.genus_data[HDP.gid_1] = HDP.genus_1
                    self.genus_data[HDP.gid_2] = HDP.genus_2
    
                    # build pidsqid lookup 
                    self.buildPidSqidLookUp(pid_sqid1, HDP.gid_1, HDP.gid_2, self.pidsqid_lookup)
                    #self.buildPidSqidLookUp(pid_sqid2, HDP.gid_1, HDP.gid_2, self.pidsqid_lookup)
                    
                    # build transfer group properties
                    self.transferGroupProperties(transfer_group_props)
                    
                    # build habitat dictionary
                    self.habitat_data[HDP.gid_1] = HDP.habitat_1
                    self.habitat_data[HDP.gid_2] = HDP.habitat_2
                    
                    # check group membership, and build data
                    transfer_group = self.checkGroupMembership(pid_sqid1)
                    
                    # get list of transfer groups containing Ecoli
                    self.getTGsContainingEcoli(HDP.genus_1, HDP.genus_2, transfer_group)
                    
                    if transfer_group:
                        if group_by_habitat:
                            if contig_len_thresh:
                                if self.checkContigLength(HDP.contigLength_1, HDP.contigLength_2, contig_len_thresh):
                                    if self.groupByHabitat(group_by_habitat, HDP.gid_1, HDP.gid_2):
                                        self.buildHitData(HDP.gid_1,
                                                          HDP.gid_2,
                                                          HDP.habitat_1,
                                                          HDP.habitat_2,
                                                          transfer_group,
                                                          phylum_thresh,
                                                          habitat_thresh,
                                                          inter_dodgey_thresh,
                                                          all_dodgey_thresh,
                                                          HDP.status_1, 
                                                          HDP.status_2,
                                                          genome_status_thresh,
                                                          remove_ecoli
                                                          )
                            else:
                                if self.groupByHabitat(group_by_habitat, HDP.gid_1, HDP.gid_2):
                                    self.buildHitData(HDP.gid_1,
                                                      HDP.gid_2,
                                                      HDP.habitat_1,
                                                      HDP.habitat_2,
                                                      transfer_group,
                                                      phylum_thresh,
                                                      habitat_thresh,
                                                      inter_dodgey_thresh,
                                                      all_dodgey_thresh,
                                                      HDP.status_1, 
                                                      HDP.status_2,
                                                      genome_status_thresh,
                                                      remove_ecoli
                                                      )
                        elif contig_len_thresh:
                            if self.checkContigLength(HDP.contigLength_1, HDP.contigLength_2, contig_len_thresh):
                                self.buildHitData(HDP.gid_1,
                                                  HDP.gid_2,
                                                  HDP.habitat_1,
                                                  HDP.habitat_2,
                                                  transfer_group,
                                                  phylum_thresh,
                                                  habitat_thresh,
                                                  inter_dodgey_thresh,
                                                  all_dodgey_thresh,
                                                  HDP.status_1, 
                                                  HDP.status_2,
                                                  genome_status_thresh,
                                                  remove_ecoli
                                                  )
                        
                        else: 
                            self.buildHitData(HDP.gid_1,
                                              HDP.gid_2,
                                              HDP.habitat_1,
                                              HDP.habitat_2,
                                              transfer_group,
                                              phylum_thresh,
                                              habitat_thresh,
                                              inter_dodgey_thresh,
                                              all_dodgey_thresh,
                                              HDP.status_1, 
                                              HDP.status_2,
                                              genome_status_thresh,
                                              remove_ecoli
                                              )
                    
class Plot(object):
    def __init__(self,hitdata):
        self.HD = hitdata
        # plotting objects
        # node properties
        self.genus_node_size                 = {}
        self.genus_node_size_values          = []
        self.phylumCols                      = {}
        self.phylumDict                      = {}
        self.phylum_colour_values            = []
        self.habitatCols                     = {}
        self.habitatDict                     = {}
        self.habitat_colour_values           = []
        # edge properties
        self.edgewidth                       = []
        self.edgeColsPlasmid                 = {}
        self.edge_colour_values_plasmid      = []
        self.edgeColsPhage                   = {}
        self.edge_colour_values_phage        = []
        self.edgeColsTransposon              = {}
        self.edge_colour_values_transposon   = []
        self.edgeColsAntibiotic              = {}
        self.edge_colour_values_antibiotic   = []
        self.edgeColsTransferGroup           = {}
        self.edge_colour_transfer_group      = []
        # directed edge properites            
        self.bidirected_edges                = []
        self.directed_edges                  = []
        self.bidirected_edges_cols           = []
        self.directed_edges_cols             = []
        # network plot
        #self.G                              = nx.Graph()
        # ColorBrewer colours
        cb2                         = CB2()
        col_set                     = "qualSet1"
        col_set_gradient            = "seqReds"
        self.ColBrewColours         = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient = cb2.maps[col_set_gradient].values()[0:10]
        
    def colourNodesByPhylum(self, genus):
        try:
            self.phylumCols[genus] = self.phylumDict[self.HD.phylum_data[genus]]
        except KeyError:
            if self.HD.phylum_data[genus] == "p__Bacteroidetes":
                self.phylumDict[self.HD.phylum_data[genus]] =  self.ColBrewColours[0] #"#ea00ff"
                self.phylumCols[genus] = self.ColBrewColours[0] #"#ea00ff"

            elif self.HD.phylum_data[genus] == "p__Lentisphaerae":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[1] #"#ffb700"
                self.phylumCols[genus] = self.ColBrewColours[1] #"#ffb700"
                
            elif self.HD.phylum_data[genus] == "p__Firmicutes":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[2] #"#0047ff"
                self.phylumCols[genus] = self.ColBrewColours[2] #"#0047ff"
                
            elif self.HD.phylum_data[genus] == "p__Spirochaetes":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[3] #"#14ff00"
                self.phylumCols[genus] = self.ColBrewColours[3] #"#14ff00"   
                
            elif self.HD.phylum_data[genus] == "p__Synergistetes":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[4] #"#6600CC"
                self.phylumCols[genus] = self.ColBrewColours[4] #"#6600CC"
                
            elif self.HD.phylum_data[genus] == "p__Actinobacteria":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[5] #"#ffff00"
                self.phylumCols[genus] = self.ColBrewColours[5] #"#ffff00"
                
            elif self.HD.phylum_data[genus] == "p__Tenericutes":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[6] #"#006600"
                self.phylumCols[genus] = self.ColBrewColours[6] #"#0080ff"
                
            elif self.HD.phylum_data[genus] == "p__Fusobacteria":
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[7] #"#00e0ff"
                self.phylumCols[genus] = self.ColBrewColours[7] #"#00e0ff"
                
            elif self.HD.phylum_data[genus] == "p__Epsilonmicrobia": 
                self.phylumDict[self.HD.phylum_data[genus]] = self.ColBrewColours[8]
                self.phylumCols[genus] = self.ColBrewColours[8] 
            
            elif self.HD.phylum_data[genus] == "p__Proteobacteria":
                self.phylumDict[self.HD.phylum_data[genus]] = "#01B5B5" # self.ColBrewColours[8] #"#ff1e00"
                self.phylumCols[genus] = "#01B5B5" #self.ColBrewColours[8] #"#ff1e00"
              
    def colourNodesByHabitat(self,genus):
        try:
            self.habitatCols[genus] = self.habitatDict[self.HD.habitat_homog_data[genus]]
        except KeyError:
            if self.HD.habitat_homog_data[genus] == "Gastrointestinal tract":
                self.habitatDict[self.HD.habitat_homog_data[genus]] =  self.ColBrewColours[0] #"#ea00ff"
                self.habitatCols[genus] = self.ColBrewColours[0] #"#ea00ff"

            elif self.HD.habitat_homog_data[genus] == "Eye":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[1] #"#ffb700"
                self.habitatCols[genus] = self.ColBrewColours[1] #"#ffb700"
                
            elif self.HD.habitat_homog_data[genus] == "Oral":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[2] #"#0047ff"
                self.habitatCols[genus] = self.ColBrewColours[2] #"#0047ff"

            elif self.HD.habitat_homog_data[genus] == "Airways":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[3] #"#14ff00"
                self.habitatCols[genus] = self.ColBrewColours[3] #"#14ff00"

            elif self.HD.habitat_homog_data[genus] == "internal_organs":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[4] #"#6600CC"
                self.habitatCols[genus] = self.ColBrewColours[4] #"#6600CC"

            elif self.HD.habitat_homog_data[genus] == "Blood":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[5] #"#ffff00"
                self.habitatCols[genus] = self.ColBrewColours[5] #"#ffff00"

            elif self.HD.habitat_homog_data[genus] == "skin":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[6] #"#006600"
                self.habitatCols[genus] = self.ColBrewColours[6] #"#0080ff"

            elif self.HD.habitat_homog_data[genus] == "Urogenital tract":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[7] #"#00e0ff"
                self.habitatCols[genus] = self.ColBrewColours[7] #"#00e0ff"
                
            elif self.HD.habitat_homog_data[genus] == "Plant":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = self.ColBrewColours[8] #"#00e0ff"
                self.habitatCols[genus] = self.ColBrewColours[8] #"#00e0ff"

            elif self.HD.habitat_homog_data[genus] == "Ear":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = '#00CCCC' #self.ColBrewColours[8] #"#ff1e00"
                self.habitatCols[genus] = '#00CCCC' #self.ColBrewColours[8] #"#ff1e00"
            
            elif self.HD.habitat_homog_data[genus] == "mixed":
                self.habitatDict[self.HD.habitat_homog_data[genus]] = '#C9C9C9' #self.ColBrewColours[8] #"#ff1e00"
                self.habitatCols[genus] = '#C9C9C9' #self.ColBrewColours[8] #"#ff1e00
   
    def colourEdgeByWrapper(self, G, colour_edge_by, genus1, genus2, transfer_group):
        if colour_edge_by.lower() == "plasmid": 
            self.colourEdgeBy(self.edgeColsPlasmid,
                              self.HD.AD.plasmid_anno,
                              genus1,
                              genus2)
            self.edge_colour_values_plasmid = [self.edgeColsPlasmid.get(edge) for edge in G.edges()]
            
        elif colour_edge_by.lower() == "transposon":
            self.colourEdgeBy(self.edgeColsTransposon,
                              self.HD.AD.transposon_anno,
                              genus1,
                              genus2)
            self.edge_colour_values_transposon = [self.edgeColsTransposon.get(edge) for edge in G.edges()]
            
        elif colour_edge_by.lower() == "phage":
            self.colourEdgeByPhage(self.edgeColsPhage,
                                   self.HD.AD.phage_anno,
                                   genus1,
                                   genus2)
            self.edge_colour_values_phage = [self.edgeColsPhage.get(edge) for edge in G.edges()]
            
        elif colour_edge_by.lower() == "antibiotic":
            self.colourEdgeBy(self.edgeColsAntibiotic,
                              self.HD.AD.antibiotic_anno,
                              genus1,
                              genus2)
            self.edge_colour_values_antibiotic = [self.edgeColsAntibiotic.get(edge) for edge in G.edges()]
            
        elif colour_edge_by.lower() == "transfer_groups":
            self.colourEdgeByTransferGroup(self.edgeColsTransferGroup,
                                           genus1,
                                           genus2,
                                           transfer_group)
            self.edge_colour_transfer_group = [self.edgeColsTransferGroup.get(edge) for edge in G.edges()]
              
    def colourEdgeByTransferGroup(self, edge_colour_dict, genus1, genus2, transfer_group):
        try:
            for pidsqid in self.HD.pidsqid_lookup[(genus1,genus2)]:
                if self.HD.checkGroupMembership(pidsqid) == transfer_group:
                    edge_colour_dict[(genus1,genus2)] = "#FF6969"
                    edge_colour_dict[(genus2,genus1)] = "#FF6969"
                    break
                else:
                    edge_colour_dict[(genus1,genus2)] = "#E1E1E1"
                    edge_colour_dict[(genus2,genus1)] = "#E1E1E1"
        except KeyError:
            edge_colour_dict[(genus1,genus2)] = "#E1E1E1"
            edge_colour_dict[(genus2,genus1)] = "#E1E1E1"
                         
    def colourEdgeByPhage(self, edge_colour_dict, anno_dict, genus1, genus2):
        max_phage_annotation_per_TG = 0
        try:
            for pidsqid in self.HD.pidsqid_lookup[(genus1,genus2)]:
                if anno_dict[pidsqid] == 0:
                    pass
                else:
                    if len(anno_dict[pidsqid]) > max_phage_annotation_per_TG:
                        max_phage_annotation_per_TG = len(anno_dict[pidsqid])
        except KeyError:
            edge_colour_dict[(genus1,genus2)] = "#C7C7C7"#"#E1E1E1" #"#373737"
            edge_colour_dict[(genus2,genus1)] = "#C7C7C7"#"#E1E1E1" #"#373737"
        
        if max_phage_annotation_per_TG ==1: 
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[2]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[2]
            
        elif max_phage_annotation_per_TG ==2:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[5]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[5]
            
        elif max_phage_annotation_per_TG ==3:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[3]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[3]
            
        elif max_phage_annotation_per_TG ==4:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[4]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[4]
            
        elif max_phage_annotation_per_TG ==5:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[5]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[5]
            
        elif max_phage_annotation_per_TG ==6:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[6]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[6]
            
        elif max_phage_annotation_per_TG ==7:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[7]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[7]
            
        elif max_phage_annotation_per_TG >=8:
            edge_colour_dict[(genus1,genus2)] = self.colBrewColoursGradient[8]
            edge_colour_dict[(genus2,genus1)] = self.colBrewColoursGradient[8]
            
        else:
            edge_colour_dict[(genus1,genus2)] = "#C7C7C7"#"#E1E1E1" #"#373737"
            edge_colour_dict[(genus2,genus1)] = "#C7C7C7"#"#E1E1E1" #"#373737"
                        
    def colourEdgeBy(self, edge_colour_dict, anno_dict, genus1, genus2):
        try:
            for pidsqid in self.HD.pidsqid_lookup[(genus1,genus2)]:
                if anno_dict[pidsqid] >0: 
                    edge_colour_dict[(genus1,genus2)] = "#FF6969"
                    edge_colour_dict[(genus2,genus1)] = "#FF6969"
                    break
                else:
                    edge_colour_dict[(genus1,genus2)] = "#E1E1E1"
                    edge_colour_dict[(genus2,genus1)] = "#E1E1E1"
        except KeyError:
            edge_colour_dict[(genus1,genus2)] = "#E1E1E1"
            edge_colour_dict[(genus2,genus1)] = "#E1E1E1"
        
    def buildNetworkData(self, colour_edge_by, transfer_group, G):
        """build data to be used to create network plot"""
        
        print self.HD.genus_list.keys()
        
        # create list of genus nodes
        genus_nodes = self.HD.genus_list.keys()
        
        # add nodes 
        G.add_nodes_from(genus_nodes)
        
        # build node properties
        for genus in genus_nodes:
            
            # Size nodes by no. of genomes
            self.genus_node_size[genus] = len(self.HD.genus_list[genus])*100
        
            # Colour nodes by phylum
            self.colourNodesByPhylum(genus)
            
            # colour nodes by habitat
            self.colourNodesByHabitat(genus)

            
        # build node data
        self.genus_node_size_values = [self.genus_node_size.get(node) for node in G.nodes()]
        self.phylum_colour_values   = [self.phylumCols.get(node) for node in G.nodes()]
        self.habitat_colour_values   = [self.habitatCols.get(node) for node in G.nodes()]
        
        # build edge properties
        for genus1 in self.HD.hit_data.keys():
            for genus2 in self.HD.hit_data[genus1]:
                
                # build directed/bidirected arrays
                self.buildDirectedEdges(genus1, genus2, G) 
                
                # build edge data
                self.colourEdgeByWrapper(G,
                                         colour_edge_by,
                                         genus1,
                                         genus2,
                                         transfer_group)                
        
        # Edgewidth is # of transfer groups
        for (u,v,d) in G.edges(data=True):
            self.edgewidth.append(int(str(d).split(':')[-1].split()[0].split('}')[0]))
        
    def buildDirectedEdges(self, genus1, genus2, G):
        if self.HD.genus_flux[genus1][genus2] > 0.6: 
            self.directed_edges.append((genus2,genus1))
            G.add_edge(genus2,
                       genus1,
                       capacity = len(self.HD.hit_data[genus1][genus2]))
        
        elif self.HD.genus_flux[genus1][genus2] < 0.4:
            self.directed_edges.append((genus1,genus2))
            G.add_edge(genus1,
                       genus2,
                       capacity = len(self.HD.hit_data[genus1][genus2]))
                    
        else:
            self.bidirected_edges.append((genus1,genus2))
            G.add_edge(genus1,
                       genus2,
                       capacity = len(self.HD.hit_data[genus1][genus2]))
            
    def drawNetworkEdges(self, colour_edge_by, pos, G):
        if colour_edge_by.lower() == "plasmid": 
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = self.edge_colour_values_plasmid,
                                   width=self.edgewidth)

        elif colour_edge_by.lower() == "transposon":
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = self.edge_colour_values_transposon,
                                   width=self.edgewidth)
            
        elif colour_edge_by.lower() == "phage":
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = self.edge_colour_values_phage,
                                   width=self.edgewidth)
            
        elif colour_edge_by.lower() == "antibiotic":
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = self.edge_colour_values_antibiotic,
                                   width=self.edgewidth)
            
        elif colour_edge_by.lower() == "transfer_groups":
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = self.edge_colour_transfer_group,
                                   width=self.edgewidth) 
        else:
            nx.draw_networkx_edges(G,
                                   pos,
                                   edge_color = "#E1E1E1",
                                   width=self.edgewidth)
            """
            # directed edges
            #nx.draw_networkx_edges(G,
                                   #pos,
                                   #edgelist   = self.directed_edges,
                                   #edge_color = "#E1E1E1")
            # bidirected edges
            #nx.draw_networkx_edges(G,
                                   #pos,
                                   #edgelist   = self.bidirected_edges,
                                   #edge_color = "#E1E1E1",
                                   #arrows=False
                                   #)#width=self.edgewidth)
            """

    def drawNetworkNodes(self, colour_node_by, pos, G):
        if colour_node_by.lower() == "phylum": 
            nx.draw_networkx_nodes(G,
                                   pos,
                                   linewidths=0,
                                   node_size= self.genus_node_size_values,
                                   node_color = self.phylum_colour_values, 
                                   alpha=0.7,
                                   with_labels=True)
            
        elif colour_node_by.lower() == "habitat":
            nx.draw_networkx_nodes(G,
                                   pos,
                                   linewidths=0,
                                   node_size= self.genus_node_size_values,
                                   node_color = self.habitat_colour_values, 
                                   alpha=0.7,
                                   with_labels=True)
            
        else:
            print "Please select node colour scheme: phylum or habitat"
            sys.exit()

    def networkPlot(self, outfile, colour_edge_by, colour_node_by, transfer_group, node_repulsion, outfmt, network_type, label_font_size):
        """create network plot"""
        # initialise network plot
        if network_type.lower() == 'undirected':
            G = nx.Graph()
        
        elif network_type.lower() == 'directed':
            G = nx.DiGraph()
            
        # network plot objects 
        # build network data
        self.buildNetworkData(colour_edge_by, transfer_group, G)
        
        # Build network plot
        # create figure
        fig = plt.figure(figsize=(30,15),dpi=300)
        plt.subplot(1,1,1,axisbg='white',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        
        # set position of nodes
        pos= nx.spring_layout(G,k=node_repulsion,iterations=500)
        #pos= nx.circular_layout(G)
        
        # draw network nodes
        self.drawNetworkNodes(colour_node_by,
                              pos,
                              G)
        
        if network_type.lower() == 'undirected':
            # draw network edges 
            self.drawNetworkEdges(colour_edge_by,
                                  pos,
                                  G)
        
        elif network_type.lower() == 'directed':
            # draw directed network edges
            self.drawNetworkEdges(colour_edge_by,
                                  pos,
                                  G)
            

        # draw network labels
        nx.draw_networkx_labels(G,
                                pos,
                                font_size=label_font_size,
                                font_color= 'black')
        
        # set axis properties
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')
        
        # write to file
        plt.savefig("%s" % (outfile),format="%s" % (outfmt))

    def NetworkPlotDirected(self, outfile, colour_edge_by, colour_node_by, transfer_group, node_repulsion, outfmt):
        """Create directed network plot"""
        # initialise directed network graph
        G = nx.DiGraph()

        # build network data
        self.buildNetworkData(colour_edge_by, transfer_group, G)
        
        # build node properties
        for genus in genus_nodes:
            
            # Size nodes by no. of genomes
            self.genus_node_size[genus] = len(self.HD.genus_list[genus])*100
        
            # Colour nodes by phylum
            self.colourNodesByPhylum(genus)
            
            # colour nodes by habitat
            self.colourNodesByHabitat(genus)

            
        # build node data
        self.genus_node_size_values = [self.genus_node_size.get(node) for node in self.G.nodes()]
        self.phylum_colour_values   = [self.phylumCols.get(node) for node in self.G.nodes()]
        self.habitat_colour_values   = [self.habitatCols.get(node) for node in self.G.nodes()]
        
        # build edge properties
         
    def buildNetworkDirectedData(self, colour_edge_by, transfer_group, G): 
        # create list of genus nodes
        genus_nodes = self.HD.genus_list.keys()
        
        # add nodes 
        G.add_nodes_from(genus_nodes)
         
    def networkPlotWrapper(self, outfile, colour_edge_by, colour_node_by, transfer_group, node_repulsion):
        pass

    def networkPlotGraphTool(self, outfile, colour_edge_by, colour_node_by, transfer_group, node_repulsion):
        """make network plot using Graph Tool!"""
        # initialise graph
        G = gt.Graph(directed=False)
        geneVertices   = {}
        count = 0

        # convert genes to vertices, numbered # add genes (nodes)
        for gene in self.geneList.keys():
            try:
                if len(self.blastData[gene]) > 30:
                    geneVertices[gene] = G.add_vertex()
                    self.vertexGene[count] = gene
                    count += 1
            except KeyError:
                pass

        # add nodes and edges
        for gene1 in self.blastData.keys():
            if len(self.blastData[gene1]) > 50:
                for gene2 in self.blastData[gene1]:
                    if len(self.blastData[gene2]) > 30:
                        if gene1 != gene2:
                            if gene1 in self.geneList and gene2 in self.geneList:
                                G.add_edge(geneVertices[gene1],
                                           geneVertices[gene2])
        # set pos
        pos = gt.sfdp_layout(G)
        #pos = gt.arf_layout(G, d=0.1)
        #pos = gt.get_hierarchy_control_points(G)
        #pos = gt.fruchterman_reingold_layout(G)

        # draw network plot
        gt.graph_draw(G,
                      pos=pos,
                      vertex_text=G.vertex_index,
                      output_size=(20000,20000),
                      vertex_font_size=100,
                      output=outFile,
                      nodesfirst=True,
                      fmt='png',
                      vertex_size = 1,
                      edge_marker_size = 1
                      )
     
class InterPhylaTransfers(object):
    def __init__(self, hitdata):
        self.HD = hitdata
        
    def checkPhylum(self, genus1, genus2):
        phylum1 = self.HD.phylum_data[genus1]
        phylum2 = self.HD.phylum_data[genus2]
    
        if phylum1 != phylum2:
            return True
        else: 
            return False
    
    def returnInterPhyla(self):
        for tuple in self.HD.pidsqid_lookup.keys():
            genus1 = tuple[0]
            genus2 = tuple[1]
            
            if self.checkPhylum(genus1, genus2):
                for pidsqid in self.HD.pidsqid_lookup[tuple]:
                    print "%s %s %s" % (genus1, genus2, pidsqid)
            
class PhageTransferGroups(object):
    def __init__(self, hitdata):
        self.HD                     = hitdata
        self.transfer_group_phage   = {}
        
    def buildPhageTGs(self, group_number, term):
        try: 
            self.transfer_group_phage[group_number][term] = 1
        except KeyError:
            self.transfer_group_phage[group_number] = {term: 1}
            
    def printPhageTGs(self):
        # print header
        print "transfer_group\tphage_annotations"
        
        for group in self.transfer_group_phage.keys():
            # join phage annotations with commas
            phage_annos = str(",".join(self.transfer_group_phage[group]))
            
            # print out dict
            print "%s\t%s" % (group, phage_annos)
                    
    def returnPhageTGs(self):
        for pidsqid in self.HD.AD.phage_anno.keys():
            # only pidsqids with phage annotations are present
            # check pidsqid membership
            group_number = self.HD.checkGroupMembership(pidsqid)
            
            if self.HD.AD.phage_anno[pidsqid] != 0 and group_number != False:
                for term in self.HD.AD.phage_anno[pidsqid]:
                
                    # add to dictionary
                    self.buildPhageTGs(group_number, term)
                
        # return phage TGs
        self.printPhageTGs()
                  
class TGTable(object):
    def __init__(self, hitdata):
        self.HD                 = hitdata
        self.genus_partners     = {}
        self.genus_TGs          = {}
        
    def buildGenusData(self, genus1, genus2, group_num): 
        try:
            self.genus_partners[genus1][genus2] = 1 
        except KeyError:
            self.genus_partners[genus1] = {genus2:1} 
        try: 
            self.genus_partners[genus2][genus1] = 1
        except KeyError:    
            self.genus_partners[genus2] = {genus1:1}
            
        try:
            self.genus_TGs[genus1][group_num] = 1
        except KeyError:
            self.genus_TGs[genus1] = {group_num:1}
        try:
            self.genus_TGs[genus2][group_num] = 1
        except KeyError:
            self.genus_TGs[genus2] = {group_num:1}
    
    def printGenusData(self, group_by_habitat):
        # print header 
        print "tax_string\tgenus\t#TGs\t#partner_genus"
        
        
        for genus in self.HD.TD.taxon_string.keys():
            if group_by_habitat:
                try:
                    lenny = self.HD.genus_by_habitat[genus]
                    try:
                        print "%s\t%s\t%d\t%d" % (self.HD.TD.taxon_string[genus],
                                                  genus,
                                                  len(self.genus_TGs[genus]),
                                                  len(self.genus_partners[genus]))
                    except KeyError:
                        print "%s\t%s\t%d\t%d" % (self.HD.TD.taxon_string[genus],
                                                  genus,
                                                  0,
                                                  0)
                except KeyError:
                    pass
            else: 
                try:
                    print "%s\t%s\t%d\t%d" % (self.HD.TD.taxon_string[genus],
                                              genus,
                                              len(self.genus_TGs[genus]),
                                              len(self.genus_partners[genus]))
                except KeyError:
                    print "%s\t%s\t%d\t%d" % (self.HD.TD.taxon_string[genus],
                                              genus,
                                              0,
                                              0)
        """
        # print genus data
        for genus in self.genus_partners.keys():
            print "%s\t%s\t%d\t%d" % (self.HD.TD.taxon_string[genus],
                                      genus,
                                      len(self.genus_TGs[genus]),
                                      len(self.genus_partners[genus]))
        """
    
    def returnGenusData(self, group_by_habitat):
        # build genus data
        for genus1 in self.HD.hit_data.keys():
            for genus2 in self.HD.hit_data[genus1]:
                for group_num in self.HD.hit_data[genus1][genus2]:
                    self.buildGenusData(genus1, genus2, group_num)
        
        # print genus data
        self.printGenusData(group_by_habitat)
        
                  
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
    # hit data object
    HD = HitData(args.group_file,
                 args.anno_file,
                 args.taxon_file,
                 args.kmer_file)
    
    HD.mainWrapper(args.hit_file,
                   args.transfer_group_props,
                   args.phylum_thresh,
                   args.habitat_thresh,
                   args.inter_dodgey_groups,
                   args.all_dodgey_groups,
                   args.group_by_habitat,
                   args.contig_size_thresh,
                   args.genome_status,
                   args.remove_ecoli)
    #"""
    # plot hit data
    P  = Plot(HD)
    P.networkPlot(args.out_file,
                  args.colour_edge_by,
                  args.colour_node_by,
                  args.transfer_group,
                  args.node_repulsion,
                  args.outfmt,
                  args.network_type,
                  args.label_font_size)
    #"""
    """
    # if you want to print out all of the inter-phyla pidsqids
    # Get inter_phyla data
    IPT = InterPhylaTransfers(HD)
    IPT.returnInterPhyla()
    """
    
    """
    # print out all phage transfer groups, and corresponding annotations
    PTG = PhageTransferGroups(HD)
    PTG.returnPhageTGs()
    """
    
    """
    # print out genus statistics
    TGT = TGTable(HD)
    TGT.returnGenusData(args.group_by_habitat)
    """
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hit_file', help="File containing hit data")
    parser.add_argument('group_file', help="File containing transfer groups")
    parser.add_argument('anno_file', help="File containing tab delimited annotation information")
    parser.add_argument('transfer_group_props', help="File containing phylum and habitat counts")
    parser.add_argument('taxon_file', help="File containing gid and taxonomy string")
    parser.add_argument('kmer_file', help="File containing kmer scores")
    parser.add_argument('-pt','--phylum_thresh', type=int,default=False,help="Set group phylum member limit")
    parser.add_argument('-ht','--habitat_thresh', type=int,default=False,help="Set group habitat member limit")
    parser.add_argument('-idg','--inter_dodgey_groups', default=False,help="Remove dodgey groups from inter-phylum analysis")
    parser.add_argument('-adg','--all_dodgey_groups', default=False,help="Remove dodgey groups from all transfers analysis")
    parser.add_argument('-tg','--transfer_group', type=int, default=1,help="Select which transfer group to colour edge by.")
    parser.add_argument('--group_by_habitat', default=False,help="Restrict heatmap to specific habitat: Options: 'gastrointestinal tract', 'ear', 'oral', 'eye', 'airways', 'internal_organs', 'blood', 'skin', 'urogenital tract', 'plant'.")
    parser.add_argument('-o','--out_file', default='Network_plot.png',help="Set output file, and directory.")
    parser.add_argument('--colour_edge_by', default='grey',help="Colour network edges by: plasmid; transposon; antibiotic; phage; transfer_groups. Default = grey")
    parser.add_argument('--colour_node_by', default='phylum',help="Colour network nodes by: phylum; habitat. Default = phylum.")
    parser.add_argument('-k','--node_repulsion', type=float,default=0.3,help="Set node repulsion value. (overlap) 0 - 1 (no overlap).")
    parser.add_argument('-outfmt','--outfmt', default='png',help="Set the output file format. Default: png.")
    parser.add_argument('-nt','--network_type', default='directed',help="Set the type of network plots to produce. Directed or undirected.")
    parser.add_argument('-lfs','--label_font_size', default=8, type=int, help="Node label font size.")
    parser.add_argument('-contig_size','--contig_size_thresh', default=False, type=int, help="Set a contig size threshold for transfers.")
    parser.add_argument('-gs','--genome_status',default=False, help="Set the status of genome to use: Finished, Draft or all")
    parser.add_argument('-re','--remove_ecoli',default=False, help="Remove ecoli from network: False (default), transfer_group or pidsqid")
    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
