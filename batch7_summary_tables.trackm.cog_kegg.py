#!/usr/bin/env python
###############################################################################
#
# __batch7_summary_tables.trackm.cog_kegg__.py - description!
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
import os
import errno
import numpy as np
np.seterr(all='raise')
import math as math

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KEGGData(object):
    def __init__(self,hitdata_file, taxon_file, group_file, kegg_file):
        self.HD                 = TFP.HitData(hitdata_file, taxon_file)
        self.GD                 = TFP.GroupData(group_file)
        self.KD                 = TFP.KEGGData(kegg_file)
        self.keggs               = {}
        self.keggs_phylum        = {}
        self.kegg_genera         = {}
        self.kegg_habitat        = {}
        self.kegg_tgs            = {}
        self.keggs_avg_transfer  = {}
        self.keggs_avg_contigs   = {}
        self.transfer_group     = {}
        
    def wrapper(self, outdirectory, filename, genome_status):
        
        # loop through KEGGs
        for kegg in self.KD.kegg_to_pidsqid.keys():
            for pidsqid in self.KD.kegg_to_pidsqid[kegg]:
                if pidsqid in self.GD.group_membership:
                    
                    # call hitdata with pidsqid
                    phylum1         = self.HD.phylum[pidsqid][0]
                    phylum2         = self.HD.phylum[pidsqid][1]
                    genus1          = self.HD.genus[pidsqid][0]
                    genus2          = self.HD.genus[pidsqid][1]
                    habitat1        = self.HD.habitat[pidsqid][0]
                    habitat2        = self.HD.habitat[pidsqid][1]
                    transferSize    = self.HD.transfer_size[pidsqid]
                    contigSize1     = self.HD.contig_size_gid[pidsqid][0]
                    contigSize2     = self.HD.contig_size_gid[pidsqid][1]
                    transferGroup   = self.GD.group_membership[pidsqid]
                    status1         = self.HD.genome_status[pidsqid][0]
                    status2         = self.HD.genome_status[pidsqid][1]
                    
                    
                    # check if intra or inter
                    if self.HD.intra_or_inter[pidsqid] == 'inter':
                        if self.checkGenomeStatus(status1, status2, genome_status):
                            # phylum 
                            self.addData(self.keggs_phylum, kegg, phylum1, 'inter')
                            self.addData(self.keggs_phylum, kegg, phylum2, 'inter')
                            # genus
                            self.addData(self.kegg_genera, kegg, genus1, 'inter')
                            self.addData(self.kegg_genera, kegg, genus2, 'inter')
                            # habitat
                            self.addData(self.kegg_habitat, kegg, habitat1, 'inter')
                            self.addData(self.kegg_habitat, kegg, habitat2, 'inter')
                            # transfer length
                            self.addLengthData(self.keggs_avg_transfer, kegg, transferSize, 'inter')
                            # contig length
                            self.addLengthData(self.keggs_avg_contigs, kegg, contigSize1, 'inter')
                            self.addLengthData(self.keggs_avg_contigs, kegg, contigSize2, 'inter')
                            # transfer group
                            self.addData(self.transfer_group, kegg, transferGroup, 'inter')
        
        # print out data
        self.printData(outdirectory,
                       filename)         
                
    def printData(self, outdirectory, filename):
        # output files
        outfile_inter = os.path.join(outdirectory, filename)
        #outfile_intra = os.path.join(outdirectory, "cog_centric_intra.all_TGs.csv")
        
        f1 = open(outfile_inter, 'w')
        #f2 = open(outfile_intra, 'w') 
        
        # write header
        f1.write("\t".join(["kegg",
                            "description",
                            "num_genus",
                            "num_phyla",
                            "num_habitats",
                            "transfers_avg",
                            "transfer_std",
                            "contig_avg",
                            "contig_std",
                            "transfer_groups\n"
                            ]))        
        # print inter phylum transfers
        for kegg in self.kegg_genera['inter']:
            num_genus                   = len(self.kegg_genera['inter'][kegg])
            num_phyla                   = len(self.keggs_phylum['inter'][kegg])
            num_habitats                = len(self.kegg_habitat['inter'][kegg])
            transfer_len_array          = self.keggs_avg_transfer['inter'][kegg]
            transfer_avg,transfer_std   = self.averageLength(transfer_len_array)
            contig_len_array            = self.keggs_avg_contigs['inter'][kegg]
            contig_avg,contig_std       = self.averageLength(contig_len_array)
            transfer_groups             = self.getTransferGroups(kegg)
            kegg_description            = self.KD.kegg_description[kegg]
            
            f1.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (kegg,
                                                                   kegg_description,
                                                                   num_genus,
                                                                   num_phyla,
                                                                   num_habitats,
                                                                   int(math.floor(transfer_avg)),
                                                                   int(math.floor(transfer_std)),
                                                                   int(math.floor(contig_avg)),
                                                                   int(math.floor(contig_std)),
                                                                   transfer_groups
                                                                   ))
            
    def averageLength(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)
    
    def addLengthData(self, dict, id, key, type):
        try:
            dict[type][id].append(key)
        except KeyError:
            try:
                dict[type][id] = [key]
            except KeyError:
                dict[type] = {id:[key]}
                
    def addData(self, dict, id, key, type):
        try:
            dict[type][id][key] = 1 
        except KeyError:
            try:
                dict[type][id] = {key:1}
            except KeyError:
                dict[type] = {id:{key:1}}    
    
    def getTransferGroups(self, kegg):
        if len(self.transfer_group['inter'][kegg]) > 1:
            transfer_groups = ''
            for TG in self.transfer_group['inter'][kegg].keys():
                transfer_groups += '%s;' % TG
            return transfer_groups
        else:
            return self.transfer_group['inter'][kegg].keys()[0]
    
    def checkGenomeStatus(self, status1, status2, user_set_status):
        if user_set_status.lower() == 'all':
            return True
        elif user_set_status.lower() in status1.lower() and user_set_status.lower() in status2.lower():
            return True
        
class COGData(object):
    def __init__(self, hitdata_file, taxon_file, group_file, annotation_file):
        self.HD                 = TFP.HitData(hitdata_file, taxon_file)
        self.GD                 = TFP.GroupData(group_file)
        self.AD                 = TFP.AnnotationData(annotation_file)
        self.cogs               = {}
        self.cogs_phylum        = {}
        self.cog_genera         = {}
        self.cog_habitat        = {}
        self.cog_tgs            = {}
        self.cogs_avg_transfer  = {}
        self.cogs_avg_contigs   = {}
        self.transfer_group     = {}
        
    def wrapper(self, outdirectory, filename):
        
        # loop through COGs
        for cog in self.AD.cog_to_pidsqid.keys():
            for pidsqid in self.AD.cog_to_pidsqid[cog]:
                if pidsqid in self.GD.group_membership:
                    
                    # call hitdata with pidsqid
                    phylum1         = self.HD.phylum[pidsqid][0]
                    phylum2         = self.HD.phylum[pidsqid][1]
                    genus1          = self.HD.genus[pidsqid][0]
                    genus2          = self.HD.genus[pidsqid][1]
                    habitat1        = self.HD.habitat[pidsqid][0]
                    habitat2        = self.HD.habitat[pidsqid][1]
                    transferSize    = self.HD.transfer_size[pidsqid]
                    contigSize1     = self.HD.contig_size_gid[pidsqid][0]
                    contigSize2     = self.HD.contig_size_gid[pidsqid][1]
                    transferGroup   = self.GD.group_membership[pidsqid]
                    
                    # check if intra or inter
                    if self.HD.intra_or_inter[pidsqid] == 'inter':
                        # phylum 
                        self.addData(self.cogs_phylum, cog, phylum1, 'inter')
                        self.addData(self.cogs_phylum, cog, phylum2, 'inter')
                        # genus
                        self.addData(self.cog_genera, cog, genus1, 'inter')
                        self.addData(self.cog_genera, cog, genus2, 'inter')
                        # habitat
                        self.addData(self.cog_habitat, cog, habitat1, 'inter')
                        self.addData(self.cog_habitat, cog, habitat2, 'inter')
                        # transfer length
                        self.addLengthData(self.cogs_avg_transfer, cog, transferSize, 'inter')
                        # contig length
                        self.addLengthData(self.cogs_avg_contigs, cog, contigSize1, 'inter')
                        self.addLengthData(self.cogs_avg_contigs, cog, contigSize2, 'inter')
                        # transfer group
                        self.addData(self.transfer_group, cog, transferGroup, 'inter')
                    
        self.printData(outdirectory, filename)
        """                    
        else:
        
        # phylum 
        self.addData(self.cogs_phylum, cog, phylum1, 'intra')
        self.addData(self.cogs_phylum, cog, phylum2, 'intra')
        # genus
        self.addData(self.cog_genera, cog, genus1, 'intra')
        self.addData(self.cog_genera, cog, genus2, 'intra')
        # habitat
        self.addData(self.cog_habitat, cog, habitat1, 'intra')
        self.addData(self.cog_habitat, cog, habitat2, 'intra')
        # transfer length
        self.addCOGLengthData(self.cogs_avg_transfer, cog, transferSize, 'intra')
        # contig length
        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize1, 'intra')
        self.addCOGLengthData(self.cogs_avg_contigs, cog, contigSize2, 'intra')
        # transfer group
        self.addData(self.transfer_group, cog, transferGroup, 'intra')
        """     
        
    
    def printData(self, outdirectory, filename):
        # output files
        outfile_inter = os.path.join(outdirectory, filename)
        #outfile_intra = os.path.join(outdirectory, "cog_centric_intra.all_TGs.csv")
        
        f1 = open(outfile_inter, 'w')
        #f2 = open(outfile_intra, 'w') 
        
        # write header
        f1.write("\t".join(["cog",
                            "description",
                            "num_genus",
                            "num_phyla",
                            "num_habitats",
                            "transfers_avg",
                            "transfer_std",
                            "contig_avg",
                            "contig_std",
                            "transfer_groups\n"
                            ]))
        
        #f2.write("cog\tdescription\tnum_genus\tnum_phyla\tnum_habitats\ttransfers_avg\ttransfer_std\tcontig_avg\tcontig_std\ttransfer_groups\n")
       

        # print inter phylum transfers
        for cog in self.cog_genera['inter']:
            num_genus                   = len(self.cog_genera['inter'][cog])
            num_phyla                   = len(self.cogs_phylum['inter'][cog])
            num_habitats                = len(self.cog_habitat['inter'][cog])
            transfer_len_array          = self.cogs_avg_transfer['inter'][cog]
            transfer_avg,transfer_std   = self.averageLength(transfer_len_array)
            contig_len_array            = self.cogs_avg_contigs['inter'][cog]
            contig_avg,contig_std       = self.averageLength(contig_len_array)
            num_transfer_groups         = len(self.transfer_group['inter'][cog])
            description                 = self.AD.cog_description[cog]
            
            f1.write("%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\n" % (cog,
                                                                   description,
                                                                   num_genus,
                                                                   num_phyla,
                                                                   num_habitats,
                                                                   transfer_avg,
                                                                   transfer_std,
                                                                   contig_avg,
                                                                   contig_std,
                                                                   num_transfer_groups
                                                                   ))
            
        """            
        for cog in self.cog_genera['intra']:
        num_genus                   = len(self.cog_genera['intra'][cog])
        num_phyla                   = len(self.cogs_phylum['intra'][cog])
        num_habitats                = len(self.cog_habitat['intra'][cog])
        transfer_len_array          = self.cogs_avg_transfer['intra'][cog]
        transfer_avg,transfer_std   = self.averageLength(transfer_len_array)
        contig_len_array            = self.cogs_avg_contigs['intra'][cog]
        contig_avg,contig_std       = self.averageLength(contig_len_array)
        num_transfer_groups         = len(self.transfer_group['intra'][cog])
        description                 = self.AD.cog_description[cog]
        
        f2.write("%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\n" % (cog,
        description,
        num_genus,
        num_phyla,
        num_habitats,
        transfer_avg,
        transfer_std,
        contig_avg,
        contig_std,
        num_transfer_groups
        ))
        """
    
    def averageLength(self, array):
        a = np.array(array)
        return np.mean(a), np.std(a)
    
    def addLengthData(self, dict, id, key, type):
        try:
            dict[type][id].append(key)
        except KeyError:
            try:
                dict[type][id] = [key]
            except KeyError:
                dict[type] = {id:[key]}
                
    def addData(self, dict, id, key, type):
        try:
            dict[type][id][key] = 1 
        except KeyError:
            try:
                dict[type][id] = {key:1}
            except KeyError:
                dict[type] = {id:{key:1}}
    

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
    if args.type.lower() == "kegg":
        KD = KEGGData(args.hitdata_file,
                     args.taxon_file,
                     args.group_file,
                     args.annotation_file)
        KD.wrapper(args.outdirectory,
                   args.filename,
                   args.genome_status)
        
    elif args.type.lower() == "cog":
        CD = COGData(args.hitdata_file,
                     args.taxon_file,
                     args.group_file,
                     args.annotation_file)         
        CD.wrapper(args.outdirectory,
                   args.filename)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata_file', help="")
    parser.add_argument('taxon_file', help="")
    parser.add_argument('group_file', help="")
    parser.add_argument('annotation_file', help="")
    parser.add_argument('outdirectory', help="")
    parser.add_argument('filename', help="")
    parser.add_argument('type', help="COG or KEGG")
    parser.add_argument('-gs','--genome_status',default='all', help="Set the status of genome to use: Finished, Draft or all")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
