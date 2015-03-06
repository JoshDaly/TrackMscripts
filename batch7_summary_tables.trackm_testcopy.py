#!/usr/bin/env python
###############################################################################
#
# __batch7_summary_tables__.py - description!
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

# system imports 
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
from collections import OrderedDict
#import os
#import errno
import numpy as np
np.seterr(all='raise')
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure

# local imports
import trackm_file_parser as TFP 

###############################################################################
###############################################################################
###############################################################################
###############################################################################
        
class TableData(object):
    def __init__(self, group_file, taxon_file, hitdata_file, anno_file, kmer_file):
        self.GD                     = TFP.GroupData(group_file) 
        self.HD                     = TFP.HitData(hitdata_file, taxon_file)
        self.AD                     = TFP.AnnotationData(anno_file)
        self.KD                     = TFP.KmerData(kmer_file)
        self.phylum_counts          = {} 
        self.genus_counts           = {}
        self.habitat_counts         = {}
        self.TG_counts              = {}
        self.contig_length          = {}
        self.transfer_length        = {}
        self.directionality         = {}
        self.genus_to_phylum        = {}
        self.cog_data_genus         = {}
        self.most_connected_genus   = {}
        self.pidsqid_tg_counts      = {}
        self.COG_genus              = {}
        self.COG_phylum             = {}
        self.COG_habitat            = {}
        self.COG_TGs                = {}
        self.status_FF              = 0
        self.status_FD              = 0
        self.status_DD              = 0
        self.status                 = {}
    
    def test(self):
        for TG in self.GD.group_data.keys():
            # loop through pidsqids
            for pidsqid in self.GD.group_data[TG]:
                contig1  = self.HD.contig_name[pidsqid][0]
                contig2  = self.HD.contig_name[pidsqid][1]
                status1  = self.HD.status[pidsqid][0]
                status2  = self.HD.status[pidsqid][1]
                gid1     = self.HD.pidsqid_to_gid[pidsqid][0]
                gid2     = self.HD.pidsqid_to_gid[pidsqid][1]

                try:
                    self.status[status1][gid1] = 1 
                except KeyError:
                    self.status[status1] = {gid1:1}
                try:
                    self.status[status2][gid2] = 1 
                except KeyError:
                    self.status[status2] = {gid2:1}
                
                
                #self.genomeStatus(status1, status2)
                
                #print "%s\t%s\t%s" % (pidsqid,
                #                      status1,
                #                      status2)

        #print "FF%d\tFD%d\tDD%d" % (self.status_FF,
        #                            self.status_FD,
        #                            self.status_DD)
        print "Genome Status:"
        print "D: %d\tPD: %d\t F: %d" % (len(self.status['Draft']),
                                         len(self.status['Permanent Draft']),
                                         len(self.status['Finished'])
                                         )
    
    def genomeStatus(self, status1, status2):
        if status1.lower() == 'finished' and status2.lower() == 'finished':
            self.status_FF += 1 
        elif status1.lower() == 'finished' and  'draft' in status2.lower():
            self.status_FD += 1
        elif 'draft' in status1.lower() and status2.lower() == 'finished':
            self.status_FD += 1
        elif status1.lower() == 'draft' and status2.lower() == 'draft':
            self.status_DD += 1 
    
    def tableDataByCOG(self):
        # loop through TGs
        for TG in self.GD.group_data.keys():
            # loop through pidsqids
            for pidsqid in self.GD.group_data[TG]:
                
                # get pidsqid associated data
                genus1   = self.HD.genus[pidsqid][0] 
                genus2   = self.HD.genus[pidsqid][1]
                phylum1  = self.HD.phylum[pidsqid][0]
                phylum2  = self.HD.phylum[pidsqid][1]
                habitat1 = self.HD.habitat[pidsqid][0]
                habitat2 = self.HD.habitat[pidsqid][1]
                
                # build COG data
                self.buildCOGDataWrapper(TG, pidsqid, genus1, phylum1, habitat1)
                self.buildCOGDataWrapper(TG, pidsqid, genus2, phylum2, habitat2)
        
        # print out COG data
        self.printHeaderCOGData()
        self.printCOGDAta()
    
    def printHeaderCOGData(self):
        print "\t".join(["cog",
                         "genus_no.",
                         "phylum_no.",
                         "habitat_no.",
                         "TG_no."
                         ])
     
    def printCOGDAta(self):
        # loop through COGs
        for cog in self.COG_genus.keys():
            # attributes
            num_genus   = len(self.COG_genus[cog])
            num_phylum  = len(self.COG_phylum[cog])
            num_habitat = len(self.COG_habitat[cog])
            num_TGs     = len(self.COG_TGs[cog])
            
            print "%s\t%d\t%d\t%d\t%d" % (cog,
                                          num_genus,
                                          num_phylum,
                                          num_habitat,
                                          num_TGs)
                 
    def buildCOGDataWrapper(self, TG, pidsqid, genus, phylum, habitat):
        try:
            # only passes if lenny is true!
            lenny = self.AD.COGs[pidsqid]
            for cog_id in lenny:
                # cog TGs
                try:
                    self.COG_TGs[cog_id][TG] = 1
                except KeyError:
                    self.COG_TGs[cog_id] = {TG:1}
                
                # cog genus 
                try:
                    self.COG_genus[cog_id][genus] += 1 
                except KeyError:
                    try:
                        self.COG_genus[cog_id][genus] = 1
                    except KeyError:
                        self.COG_genus[cog_id] = {genus:1}
                                                
                # cog phylum
                try:
                    self.COG_phylum[cog_id][phylum] = 1
                except KeyError:
                    self.COG_phylum[cog_id] = {phylum:1}
                    
                # cog habitat
                try:
                    self.COG_habitat[cog_id][habitat] = 1
                except KeyError:
                    self.COG_habitat[cog_id] = {habitat:1}
        except KeyError:
            pass   
    
    def tableDataByGenus(self):
        
        # loop through TGs
        for TG in self.GD.group_data.keys():
            # loop through pidsqids
            for pidsqid in self.GD.group_data[TG]:
                
                # get genus data
                genus1   = self.HD.genus[pidsqid][0] 
                genus2   = self.HD.genus[pidsqid][1]
                phylum1  = self.HD.phylum[pidsqid][0]
                phylum2  = self.HD.phylum[pidsqid][1]
                habitat1 = self.HD.habitat[pidsqid][0]
                habitat2 = self.HD.habitat[pidsqid][1]
                
                # add connection to genus
                self.addConnectionToGenus(genus1, genus2, TG)
                
                # get phylum data
                self.genus_to_phylum[genus1] = phylum1
                self.genus_to_phylum[genus2] = phylum2
                self.addPhylumGenus(genus1, phylum2)
                self.addPhylumGenus(genus2, phylum1)
                 
                # add habitat
                self.addHabitatGenus(genus1, habitat1)
                self.addHabitatGenus(genus2, habitat2)
                
                # add TG
                self.addTGGenus(genus1, TG)
                self.addTGGenus(genus2, TG)
                
                # add contig length
                self.addContigLenGenus(genus1, pidsqid, 0)
                self.addContigLenGenus(genus2, pidsqid, 1)
                
                # add transfer length 
                self.addTransferLenGenus(genus1, pidsqid, 0)
                self.addTransferLenGenus(genus2, pidsqid, 1)
                
                # get kmer scores
                kmer_score1 = self.KD.kmer_scores[pidsqid]
                kmer_score2 = 1 - self.KD.kmer_scores[pidsqid]
                self.addKmerScoreGenus(genus1, kmer_score1)
                self.addKmerScoreGenus(genus2, kmer_score2)
                
                # get COG data
                self.addCOGGenus(genus1, pidsqid)
                self.addCOGGenus(genus2, pidsqid)
                
        # print header
        self.printHeader('genus')
                
        # print out genus data
        for genus in self.genus_to_phylum.keys():
            
            # get COGs ordered by frequency
            cogs = self.getCOGsOrdered(genus)
        
            # get average/std
            transfer_length_array = self.transfer_length[genus]
            transfer_len_mean,transfer_len_std  = self.averageLength(transfer_length_array)
            
            contig_length_array   = self.contig_length[genus] 
            contig_len_mean, contig_len_std = self.averageLength(contig_length_array)
            
            kmer_score_array = self.directionality[genus]
            kmer_score_mean, kmer_score_std = self.averageLength(kmer_score_array)

            # number of habitats
            number_of_habitats          = len(self.habitat_counts[genus])
            # number of different genus interactions
            number_of_different_genus   = len(self.most_connected_genus[genus])
            most_connected_genus        = self.getMostConnectedGenus(genus, self.most_connected_genus) 
            
            # number of phylum-level interactions
            phylum_interactions = len(self.phylum_counts[genus])
            
            # number of different phylum interactions
            number_of_different_phylums = len(self.phylum_counts[genus])
            
            # number of different TGs
            TGs = len(self.TG_counts[genus])
            
            self.printOutGenusData(genus,
                                   self.genus_to_phylum[genus],
                                   phylum_interactions,
                                   number_of_habitats,
                                   number_of_different_genus,
                                   TGs,
                                   most_connected_genus,
                                   transfer_len_mean,
                                   transfer_len_std,
                                   contig_len_mean,
                                   contig_len_std,
                                   kmer_score_mean,
                                   kmer_score_std,
                                   cogs)
                    
    def addConnectionToGenus(self, genus1, genus2, TG):
        try:
            self.most_connected_genus[genus1][genus2][TG] = 1 
        except:
            try: 
                self.most_connected_genus[genus1][genus2] = {TG:1}
            except KeyError:
                self.most_connected_genus[genus1] = {genus2:{TG:1}}
                
        try:
            self.most_connected_genus[genus2][genus1][TG] = 1
        except:
            try: 
                self.most_connected_genus[genus2][genus1] = {TG:1}
            except KeyError:
                self.most_connected_genus[genus2] = {genus1:{TG:1}}
                
    def getCOGsOrdered(self, genus):
        cogs = ''
        try:
            cogs = self.findCOGs(self.cog_data_genus[genus])
            cogs = self.getTop5COGs(cogs)
            return cogs
        except KeyError:
            pass
                
    def addCOGGenus(self, genus, pidsqid):
        try:
            # only passes if lenny is true!
            lenny = self.AD.COGs[pidsqid]
            for cog_id in lenny: 
                try:
                    self.cog_data_genus[genus][cog_id] += 1 
                except KeyError:
                    try:
                        self.cog_data_genus[genus][cog_id] = 1
                    except KeyError:
                        self.cog_data_genus[genus] = {cog_id:1}
        except KeyError:
            pass
              
    def addKmerScoreGenus(self, genus, kmerscore):
        try:
            self.directionality[genus].append(kmerscore)
        except KeyError:
            self.directionality[genus] = [kmerscore]
                
    def addTransferLenGenus(self, genus, pidsqid, i):
        try:
            self.transfer_length[genus].append(self.HD.transfer_size_gid[pidsqid][i])
        except KeyError:
            self.transfer_length[genus] = [self.HD.transfer_size_gid[pidsqid][i]]
    
    def addContigLenGenus(self, genus, pidsqid, i):
        try:
            self.contig_length[genus].append(self.HD.contig_size_gid[pidsqid][i])
        except KeyError:
            self.contig_length[genus] = [self.HD.contig_size_gid[pidsqid][i]]
                
    def addTGGenus(self, genus, TG):
        try:
            self.TG_counts[genus][TG] = 1
        except: 
            self.TG_counts[genus] = {TG:1}

    def addPhylumGenus(self, genus, phylum):
        try:
            self.phylum_counts[genus][phylum] = 1
        except KeyError:
            self.phylum_counts[genus] = {phylum:1}
                
    def addHabitatGenus(self, genus, habitat):
        try:
            self.habitat_counts[genus][habitat] += 1 
        except KeyError:
            try:
                self.habitat_counts[genus][habitat] = 1 
            except KeyError:
                self.habitat_counts[genus] = {habitat: 1}
             
    def tableDataByTransferGroup(self):
        # print header 
        self.printHeader('TG')
        
        for TG in self.GD.group_data.keys():
            transfer_len_mean = 0
            contig_len_mean   = 0
            # reset COG dictionary
            Top_COG_ids             = {}
            
            for pidsqid in self.GD.group_data[TG]:
                
                # get data
                genus1   = self.HD.genus[pidsqid][0] 
                genus2   = self.HD.genus[pidsqid][1]
                phylum1  = self.HD.phylum[pidsqid][0]
                phylum2  = self.HD.phylum[pidsqid][1]
                habitat1 = self.HD.habitat[pidsqid][0]
                habitat2 = self.HD.habitat[pidsqid][1]
                
                self.buildPhylumCount(TG, phylum1, phylum2)
                self.buildGenusCount(TG, genus1, genus2)
                self.buildHabitatCount(TG, habitat1, habitat2)
                self.buildContigLength(TG, pidsqid)
                self.buildTransferLength(TG, pidsqid)
                self.buildCOGData(TG, pidsqid, Top_COG_ids)
                self.buildDirectionalityData(TG, pidsqid)
                self.buildPidsqidTG(TG, pidsqid)

            # get Cogs descending order
            cogs = self.findCOGs(Top_COG_ids)
            cogs = self.getTop5COGs(cogs)
            
            # get average lengths
            transfer_length_array = self.transfer_length[TG]
            transfer_len_mean,transfer_len_std  = self.averageLength(transfer_length_array)
            
            contig_length_array   = self.contig_length[TG] 
            contig_len_mean, contig_len_std = self.averageLength(contig_length_array)
            
            kmer_score_array = self.directionality[TG]
            kmer_score_mean, kmer_score_std = self.averageLength(kmer_score_array)
            
            # get most connected genus
            most_connected_genus = self.getMostConnectedGenus(TG, self.genus_counts)
            
            # get number of pidsqids per TG
            pidsqid_number = len(self.pidsqid_tg_counts[TG])
            
            self.printOutTGTable(cogs,
                                 TG,
                                 transfer_len_mean,
                                 transfer_len_std,
                                 contig_len_mean,
                                 contig_len_std,
                                 pidsqid_number,
                                 most_connected_genus,
                                 kmer_score_mean,
                                 kmer_score_std)
    
    def buildPidsqidTG(self, TG, pidsqid):
        try:
            self.pidsqid_tg_counts[TG][pidsqid] = 1
        except:
            self.pidsqid_tg_counts[TG] = {pidsqid:1}
    
    def printOutGenusData(self, genus, phylum, phylum_interactions, habitats, genus_partners, TGs, most_connected_genus, transfer_len_mean, transfer_len_std, contig_len_mean, contig_len_std, kmer_score_mean, kmer_score_std, cogs):
        # print data
        print "%s\t%s\t%d\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%s" % (genus,
                                                                          phylum,
                                                                          habitats,
                                                                          phylum_interactions,
                                                                          genus_partners,
                                                                          TGs,
                                                                          most_connected_genus,
                                                                          transfer_len_mean,
                                                                          transfer_len_std,
                                                                          contig_len_mean,
                                                                          contig_len_std,
                                                                          kmer_score_mean,
                                                                          kmer_score_std,
                                                                          cogs)
            
    def getMostConnectedGenus(self, key, dict):
        most_connected_genus = [0,0]
        for genus in dict[key]:
            try:
                if  len(dict[key][genus]) > most_connected_genus[1]: 
                    most_connected_genus = [genus, len(dict[key][genus])]
                elif len(dict[key][genus]) == most_connected_genus[1]: 
                    most_connected_genus[0] += ";%s" % (genus)
            except TypeError:
                if  dict[key][genus] > most_connected_genus[1]: 
                    most_connected_genus = [genus, dict[key][genus]]
                elif dict[key][genus] == most_connected_genus[1]: 
                    most_connected_genus[0] += ";%s" % (genus)
        
        return most_connected_genus[0]
            
    def averageLength(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)
                
    def buildDirectionalityData(self, TG, pidsqid):
        try:
            self.directionality[TG].append(self.KD.distance_from_neutral[pidsqid])
        except KeyError:
            self.directionality[TG] = [self.KD.distance_from_neutral[pidsqid]]
                
    def buildCOGData(self, TG, pidsqid, COG_dict):
        try:
            lenny = self.AD.COGs[pidsqid]
            for cog_id in lenny: 
                try:
                    COG_dict[cog_id] += 1 
                except KeyError:
                    COG_dict[cog_id] = 1
        except KeyError:
            pass
    
    def buildContigLength(self, TG, pidsqid):
        try:
            self.contig_length[TG].append(self.HD.contig_size[pidsqid])
        except KeyError:
            self.contig_length[TG] = [self.HD.contig_size[pidsqid]]
        
    def buildTransferLength(self, TG, pidsqid):
        try:
            self.transfer_length[TG].append(self.HD.transfer_size[pidsqid])
        except KeyError:
            self.transfer_length[TG] = [self.HD.transfer_size[pidsqid]]

    def getTop5COGs(self, cogs):
        semi_colon = cogs.rstrip().split(";")
        if len(semi_colon) >= 5:
            return ";".join(semi_colon[0:5]) 
        else:
            return ";".join(semi_colon)
            
    def findCOGs(self, cog_dict):
        d = OrderedDict(sorted(cog_dict.items(),key=lambda t: t[1]))
        cog_array = []
        cogs = ''
        for key,value in d.items():
            cog_array.append("%s:%d;" % (key,
                                        value))
        cog_array.reverse()
        for i in cog_array:
            cogs+= i
        
        return cogs
        
    def printHeader(self, type):
        # header
        if type == 'TG':
            print "\t".join(["transfer_group",
                             "transfer_size_mean",
                             "transfer_size_std",
                             "contig_size_mean",
                             "contig_size_std",
                             "habitat_no.",
                             "phylum_no.",
                             "genus_no.",
                             "pidsqids_no.",
                             "most_connected_genus",
                             "kmer_score_mean",
                             "kmer_score_std",
                             "cog_ids"])    
        elif type == 'genus':
            print "\t".join(["genus",
                             "phylum",
                             "habitats",
                             "phylum_interactions",
                             "genus_partners",
                             "transfer_groups",
                             "most_connected_genus",
                             "transfer_size_mean",
                             "transfer_size_std",
                             "contig_size_mean",
                             "contig_size_std",
                             "kmer_score_mean",
                             "kmer_score_std",
                             "cog_ids"]) 
    
    def printOutTGTable(self, cogs, TG, transfer_len_mean, transfer_len_std, contig_len_mean, contig_len_std, pidsqid_no, most_connected_genus, kmer_score_mean, kmer_score_std):
        # print data
        print "%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%s\t%f\t%f\t%s" % (TG,
                                                                  transfer_len_mean,
                                                                  transfer_len_std,
                                                                  contig_len_mean,
                                                                  contig_len_std,
                                                                  len(self.habitat_counts[TG]),
                                                                  len(self.phylum_counts[TG]),
                                                                  len(self.genus_counts[TG]),
                                                                  pidsqid_no,
                                                                  most_connected_genus,
                                                                  kmer_score_mean,
                                                                  kmer_score_std,
                                                                  cogs)
            
    def buildPhylumCount(self, TG, phylum1, phylum2):
        try:
            self.phylum_counts[TG][phylum1] = 1
        except KeyError:
            self.phylum_counts[TG] = {phylum1:1}
            
        try:
            self.phylum_counts[TG][phylum2] = 1
        except KeyError:
            self.phylum_counts[TG] = {phylum2:1}
    
    def buildGenusCount(self, TG, genus1, genus2):
        try:
            self.genus_counts[TG][genus1] += 1
        except KeyError:
            try:
                self.genus_counts[TG][genus1] = 1
            except KeyError:
                self.genus_counts[TG] = {genus1:1}
            
        try:
            self.genus_counts[TG][genus2] += 1
        except KeyError:
            try:
                self.genus_counts[TG][genus2] = 1
            except KeyError:
                self.genus_counts[TG] = {genus2:1}
    
    def buildHabitatCount(self, TG, habitat1, habitat2):
        try:
            self.habitat_counts[TG][habitat1] = 1
        except KeyError:
            self.habitat_counts[TG] = {habitat1:1}
            
        try:
            self.habitat_counts[TG][habitat2] = 1
        except KeyError:
            self.habitat_counts[TG] = {habitat2:1}
    
    
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
    TD = TableData(args.group_file,
                   args.taxon_file,
                   args.hitdata_file,
                   args.anno_file,
                   args.kmer_file)
    if args.type.lower() == "transfer_group":
        TD.tableDataByTransferGroup()
    elif args.type.lower() == "genus":
        TD.tableDataByGenus()
    elif args.type.lower() == "cog":
        TD.tableDataByCOG()
    elif args.type.lower() == 'test':
        TD.test()

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
    parser.add_argument('kmer_file', help="File containing Kmer directionality data")
    parser.add_argument('-t','--type',default='transfer_group', help="File containing Kmer directionality data")
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
