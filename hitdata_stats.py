#!/usr/bin/env python
###############################################################################
#
# __hitdata_stats__.py - description!
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
import matplotlib.pyplot as plt
#import numpy as np
#np.seterr(all='raise')

# local imports
import trackm_file_parser as TFP
from cb2cols import Cb2Cols as CB2

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HitDataStats(object):
    def __init__(self, hitdata, pairs_file, transfer_groups, contam_pidsqids, taxon_file):
        self.HD         = TFP.HitData(hitdata)
        self.TG         = TFP.GroupData(transfer_groups)
        self.PD         = TFP.PairsFileData(pairs_file)
        if taxon_file:
            self.TD         = TFP.TaxonomyData(taxon_file)
        if contam_pidsqids:
            self.CP         = TFP.ContaminatedPidsqids(contam_pidsqids) 
        self.contaminated_TGs       = {}
        self.contaminated_pidsqids  = {}
        # ColorBrewer colours
        cb2                                 = CB2()
        col_set                             = "qualSet1"
        col_set_gradient                    = "seqReds"
        self.ColBrewColours                 = cb2.maps[col_set].values()[0:10]
        self.colBrewColoursGradient         = cb2.maps[col_set_gradient].values()[0:9]
        
    def wrapper(self, type, outfile, outfmt, contam_pidsqids, data_type):
        if type == 'phylum_interactions':
            self.createInteractionMatrix(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'genus_interactions':
            self.createInteractionMatrix(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'phylum_hit_frequency':
            self.calculateInteractionFrequency(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'genus_hit_frequency':
            self.calculateInteractionFrequency(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'habitat_hit_frequency':
            self.calculateInteractionFrequency(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'genome_status_hit_frequency':
            self.calculateInteractionFrequency(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'phylum_hit_count':
            self.calculatePhylumHitCount(outfile, outfmt, contam_pidsqids, type, data_type)
        
        elif type == 'hitdata_fix':
            self.fixTransferGroupsHitData()
        
        elif type == 'PhiX':
            self.phixStats()
    
    def createInteractionMatrix(self, outfile, outfmt, contam_pidsqids, type, data_type):
        # determine contaminated transfer groups
        if contam_pidsqids:
            self.evaluateTransferGroups()
        
        compairsons     = {}
        interactions    = {}
        subjects        = {}
        matrix          = []
        
        for pidsqid in self.HD.hit_data.keys():
            gid1 = self.HD.pidsqid_to_gid[pidsqid][0]
            gid2 = self.HD.pidsqid_to_gid[pidsqid][1]
            subject1 = ''
            subject2 = ''
            
            if type == 'phylum_interactions':
                subject1 = self.HD.phylum[pidsqid][0]
                subject2 = self.HD.phylum[pidsqid][1]
            elif type == 'genus_interactions':
                #subject1 = self.HD.genus[pidsqid][0]
                #subject2 = self.HD.genus[pidsqid][1]
                subject1 = self.TD.taxon_string_genus[gid1]
                subject2 = self.TD.taxon_string_genus[gid2]
            
            self.addComparisons(gid1, subject1, compairsons)
            self.addComparisons(gid2, subject2, compairsons)
            
            # add subject
            subjects[subject1] = 1
            subjects[subject2] = 1
            
            try:
                transfer_group = int(self.TG.group_membership[pidsqid])
                if self.checkForContamination(data_type, pidsqid, transfer_group):
                    
                    # add interaction data
                    self.addSubject(subject1, subject2, transfer_group, interactions)
                    self.addSubject(subject2, subject1, transfer_group, interactions)
            except KeyError:
                pass
        
        array = subjects.keys()
        array.sort()
        
        # add zeros to matrix
        self.zeroArray(matrix, array)
        
        for i in range(len(array)):
            hits_to_append = []
            for v in range(i, len(array)):
                subject1 = array[i]
                subject2 = array[v]
                try:
                    hits = len(interactions[subject1][subject2])
                except KeyError:
                    hits = 0
                matrix[i][v] = hits
        
        # make scatter plot
        fig = plt.figure(figsize=(30,15),dpi=300)
        ax = plt.subplot(1,1,1,axisbg='white',autoscale_on=True, aspect='equal')#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        xs      = []
        ys      = []
        colours = []
        
        for i in range(len(matrix)):
            for v in range(len(matrix)):
                hits = matrix[i][v]
                xs.append(i)
                ys.append(v)
                self.colourByHits(hits, colours)
                
        #set the number of tick labels
        ticks = []
        for i in range(len(array)):
            ticks.append(i)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        
        # set labels as phyla
        labels = [item.get_text() for item in ax.get_xticklabels()]
        for i,v in enumerate(labels):
            labels[i] = array[i]
        
        ax.set_xticklabels(labels, fontsize = 5) # phylum=x-large genus=small 
        ax.set_yticklabels(labels, fontsize = 5) # phylum=x-large genus=small
        
        # rotate x labels 90 degrees
        plt.xticks(rotation=90)
        
        plt.tight_layout()
        
        plt.scatter(xs, ys, c = colours, s=20, marker='s',linewidths=0) # s=7000, s=10
        
        # plot hits over scatter boxes
        for i in range(len(matrix)):
            for v in range(len(matrix)):
                hits = matrix[i][v]
                #plt.text(i,v,"%d" % hits, color="black", ha='center', va='center', backgroundcolor='white', fontsize=1)
        
        outfile = "%s.%s" % (outfile, outfmt) 
        
        plt.savefig("%s" % (outfile),format="%s" % (outfmt), dpi = 300)
        
        """
        # print header
        header = " " # empty top left corner
        for i in phylum_array:
            header += "\t%s" % i
        print header
            
        for i in range(len(phylum_array)):
            string_to_print = "%s" % phylum_array[i]
            for v in range(len(phylum_array)):
                string_to_print += "\t%d" % matrix[i][v]
            print string_to_print
        """
        
    def addComparisons(self, gid, subject, dict):
        num_comparisons = len(self.PD.pairs_data[gid])
        try:
            dict[subject][gid] = num_comparisons
        except KeyError:
            dict[subject] = {gid:num_comparisons}
    
    def evaluateTransferGroups(self):
        for pidsqid in self.CP.contam_pidsqids.keys():
            self.contaminated_pidsqids[pidsqid] = 1
            try:
                transfer_group = int(self.TG.group_membership[pidsqid])
                self.contaminated_TGs[transfer_group] = 1
            except KeyError:
                pass
    
    def colourByHits(self, hits, colours):
        if hits == 0:
            colours.append('white')
        elif hits >0 and hits <= 5: # elif hits >0 and hits <= 10:  
            colours.append(self.colBrewColoursGradient[1])
        elif hits >5 and hits <= 10: # elif hits >10 and hits <= 20:
            colours.append(self.colBrewColoursGradient[2])
        elif hits >10 and hits <= 15: # elif hits >20 and hits <= 30:
            colours.append(self.colBrewColoursGradient[3])
        elif hits >15 and hits <= 20: #elif hits >30 and hits <= 40:
            colours.append(self.colBrewColoursGradient[4])
        elif hits >20 and hits <= 25: #elif hits >40 and hits <= 50:
            colours.append(self.colBrewColoursGradient[5])
        elif hits >25 and hits <= 30: #elif hits >50 and hits <= 60:
            colours.append(self.colBrewColoursGradient[6])
        elif hits >30 and hits <= 40: #elif hits >60 and hits <= 70:
            colours.append(self.colBrewColoursGradient[7])
        elif hits >40:  #elif hits >70:
            colours.append(self.colBrewColoursGradient[8])
    
    def checkPidsqid(self, pidsqid):
        try: 
            if pidsqid in self.CP.contam_pidsqids:
                return False
            else: 
                return True
        except AttributeError:
            return True
    
    def checkForContamination(self, data_type, pidsqid, transfer_group):
        if data_type == 'transfer_group':
            try: 
                if transfer_group in self.contaminated_TGs:
                    return False
                else:
                    return True
            except AttributeError:
                return True
        elif data_type == 'pidsqid':
            try: 
                if pidsqid in self.contaminated_pidsqids:
                    return False
                else:
                    return True
            except AttributeError:
                return True
    
    def zeroArray(self, matrix, array):
        for i in range(len(array)):
            zeros = []
            for v in range(len(array)):
                zeros.append(0)
            matrix.append(zeros)
            
    def addSubject(self, subject1, subject2, TG, dict):
        try:
            dict[subject1][subject2][TG] =1 
        except KeyError:
            try:
                dict[subject1][subject2] = {TG:1} 
            except KeyError:
                dict[subject1] = {subject2:{TG:1}}
            
    def addSubjectHit(self, subject1, gid, data_type, pidsqid, transfer_group, dict):
        intra_or_inter = self.HD.intra_or_inter[pidsqid]
        if data_type == 'pidsqid':
            try:
                dict[subject1][gid][intra_or_inter][pidsqid] =1 
            except KeyError:
                try:
                    dict[subject1][gid][intra_or_inter] = {pidsqid:1}
                except KeyError:
                    try:
                        dict[subject1][gid] = {intra_or_inter:{pidsqid:1}}
                    except KeyError:
                        dict[subject1] = {gid:{intra_or_inter:{pidsqid:1}}}
        
        elif data_type == 'transfer_group':
            try:
                dict[subject1][gid][intra_or_inter][transfer_group] =1 
            except KeyError:
                try:
                    dict[subject1][gid][intra_or_inter] = {transfer_group:1}
                except KeyError:
                    try:
                        dict[subject1][gid] = {intra_or_inter:{transfer_group:1}}
                    except KeyError:
                        dict[subject1] = {gid:{intra_or_inter:{transfer_group:1}}}
    
    def fixTransferGroupsHitData(self):
        # print header 
        print '\t'.join(['pidsqid','transfer_group','hid','pid','ani_comp','ident','gid_1','habitat_1','phylum_1','genus_1','status_1','sequencingMethod_1','sequencingCentre_1','horizontalTransferred_1','genomeSize_1','scaffoldCount_1','len_1','strand_1','cid_1','contig_1','contigLength_1','sqid_1','start_1','gid_2','habitat_2','phylum_2','genus_2','status_2','sequencingMethod_2','sequencingCentre_2','horizontalTransferred_2','genomeSize_2','scaffoldCount_2','len_2','strand_2','cid_2','contig_2','contigLength_2','sqid_2','start_2'])
        
        for pidsqid in self.HD.hit_data.keys():
            string_to_print = ''
            string_to_print += '%s' % pidsqid
            try:
                transfer_group = self.TG.group_membership[pidsqid]
            except KeyError:
                transfer_group = 'NA'
            for i,v in enumerate(self.HD.hit_data[pidsqid]):
                if i == 0:
                    string_to_print += '\t%s' % transfer_group
                else:
                    string_to_print += '\t%s' % v
            print string_to_print
                
    def calculateInteractionFrequency(self, outfile, outfmt, contam_pidsqids, type, data_type):
        if contam_pidsqids:
            self.evaluateTransferGroups()
        else:
            print "################ ERROR ####################"
            print 'Please provide contaminated Pidsqids file'
            print '-cp /PATHTOFILE'
            os.sys.exit()
        
        comparisons     = {}
        hits_clean      = {}
        hits_all        = {}
        subjects        = {}
        
        for pidsqid in self.HD.hit_data.keys():
            gid1 = self.HD.pidsqid_to_gid[pidsqid][0]
            gid2 = self.HD.pidsqid_to_gid[pidsqid][1]
            subject1 = ''
            subject2 = ''
            
            if type == 'phylum_hit_frequency':
                subject1 = self.HD.phylum[pidsqid][0]
                subject2 = self.HD.phylum[pidsqid][1]
            
            elif type == 'genus_hit_frequency':
                subject1 = self.HD.genus[pidsqid][0]
                subject2 = self.HD.genus[pidsqid][1]
            
            elif type == 'habitat_hit_frequency':
                subject1 = self.HD.habitat[pidsqid][0]
                subject2 = self.HD.habitat[pidsqid][1]
            
            elif type == 'genome_status_hit_frequency':
                subject1 = self.HD.status[pidsqid][0]
                subject2 = self.HD.status[pidsqid][1]
            
            self.addComparisons(gid1, subject1, comparisons)
            self.addComparisons(gid2, subject2, comparisons)
            
            try:
                transfer_group = int(self.TG.group_membership[pidsqid])
                
                # add phylum or genus
                subjects[subject1] = 1
                subjects[subject2] = 1
                
                # check if contaminated transfer group = cleaned data
                if self.checkForContamination(data_type, pidsqid, transfer_group):
                    
                    # add interaction data
                    self.addSubjectHit(subject1, gid1, data_type, pidsqid, transfer_group, hits_clean)
                    self.addSubjectHit(subject2, gid2, data_type, pidsqid, transfer_group, hits_clean)
            
                # all data
                # add interaction data
                self.addSubjectHit(subject1, gid1, data_type, pidsqid, transfer_group, hits_all)
                self.addSubjectHit(subject2, gid2, data_type, pidsqid, transfer_group, hits_all)
            
            except KeyError:
                pass
    
        # open outfile
        of = open(outfile,'w')
        
        # write header
        of.write("taxon_rank\tintra_all\tintra_clean\tintra_perc\tinter_all\tinter_clean\tinter_perc\n")
    
        # calculate total number of TGs and comparisons
        for subject in hits_all.keys():
            total_hits_all_intra    = 0
            total_hits_all_inter    = 0
            total_hits_clean_intra  = 0
            total_hits_clean_inter  = 0
            total_comparisons = 0
            for gid in comparisons[subject]:
                num_comparisons         = comparisons[subject][gid]
                total_comparisons      += num_comparisons
                num_hits_all_intra      = self.getHits(subject, gid, 'intra', hits_all)
                num_hits_all_inter      = self.getHits(subject, gid, 'inter', hits_all)
                num_hits_clean_intra    = self.getHits(subject, gid, 'intra', hits_clean)
                num_hits_clean_inter    = self.getHits(subject, gid, 'inter', hits_clean)
            
                # accumulate hits
                total_hits_all_intra         += num_hits_all_intra
                total_hits_all_inter         += num_hits_all_inter
                total_hits_clean_intra       += num_hits_clean_intra
                total_hits_clean_inter       += num_hits_clean_inter
            
            # calculate frequency
            hit_frequency_all_intra     = total_hits_all_intra/float(total_comparisons) * 100
            hit_frequency_all_inter     = total_hits_all_inter/float(total_comparisons) * 100
            hit_frequency_clean_intra   = total_hits_clean_intra/float(total_comparisons) * 100
            hit_frequency_clean_inter   = total_hits_clean_inter/float(total_comparisons) * 100
            try:
                intra_perc              = hit_frequency_clean_intra/float(hit_frequency_all_intra)
            except ZeroDivisionError:
                intra_perc              = 0
            try:
                inter_perc              = hit_frequency_clean_inter/float(hit_frequency_all_inter)
            except ZeroDivisionError:
                inter_perc              = 0
            
            of.write("%s\t%f\t%f\t%.2f\t%f\t%f\t%.2f\n" % (subject,
                                                           hit_frequency_all_intra,
                                                           hit_frequency_clean_intra,
                                                           intra_perc,
                                                           hit_frequency_all_inter,
                                                           hit_frequency_clean_inter,
                                                           inter_perc)) 
        
        # close outfile
        of.close()
        
    def getHits(self, subject, gid, intra_or_inter, dict):
        num_hits = 0
        try:
            num_hits = len(dict[subject][gid][intra_or_inter])
        except KeyError:
            num_hits = 0
        return num_hits
    
    def calculatePhylumHitCount(self, outfile, outfmt, contam_pidsqids, type, data_type):
        inter_phyla_TGs         = {}
        contam_inter_phyla_TGs  = {}
        
        # loop through pidsqids
        for pidsqid in self.HD.hit_data.keys():
            
            transfer_group = self.HD.hit_data[pidsqid][0]
            
            if transfer_group != 'NA':
            
                if self.HD.intra_or_inter[pidsqid] == 'inter':
                    inter_phyla_TGs[transfer_group] = 1
                    
                    if pidsqid in self.CP.contam_pidsqids:
                        contam_inter_phyla_TGs[transfer_group] = 1 
                
        # print out the number of contaminated TGs 
        print "The total number of TGs containing inter-phyla hits: %d" % len(inter_phyla_TGs)
        print "The number of contaminted inter-phyla TGs: %d" % len(contam_inter_phyla_TGs)
        
         
    def phixStats(self):
        
        # what sequencing platforms are affected?
        seq_plts    = {}
        genus       = {}
        seq_centre  = {}
        
        for pidsqid in self.HD.hit_data.keys():
            # transfer group 340 == phiX
            transfer_group = self.HD.hit_data[pidsqid][0]
            
            if transfer_group == '33':
                sequence_platform1  = self.HD.hit_data[pidsqid][10]
                sequence_platform2  = self.HD.hit_data[pidsqid][27]
                sequence_centre1    = self.HD.hit_data[pidsqid][11]
                sequence_centre2    = self.HD.hit_data[pidsqid][28]
                genus1              = self.HD.genus[pidsqid][0] 
                genus2              = self.HD.genus[pidsqid][1]
                # add to dict
                seq_plts[sequence_platform1] = 1
                seq_plts[sequence_platform2] = 1
                genus[genus1] = 1
                genus[genus2] = 1
                seq_centre[sequence_centre1] = 1
                seq_centre[sequence_centre2] = 1
                
        for centre in seq_centre.keys():
            print centre
        
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
    HDS = HitDataStats(args.hitdata,
                       args.pairs_file,
                       args.transfer_groups_file,
                       args.contaminted_pidsqids,
                       args.taxon_file)
    HDS.wrapper(args.type,
                args.outfile,
                args.outfmt,
                args.contaminted_pidsqids,
                args.data_type)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitdata', help="")
    parser.add_argument('pairs_file', help="")
    parser.add_argument('-tf','--taxon_file', default=False, help="")
    parser.add_argument('-tg','--transfer_groups_file', help="")
    parser.add_argument('-dt','--data_type', default='transfer_group', help="transfer_group[default] or pidsqid")
    parser.add_argument('-t','--type', default = 'phylum_hit_frequency', help="Type of summary table to create. phylum_interactions, genus_interactions, habitat_hit_frequency, genus_hit_frequency, phylum_hit_frequency(default), genome_status_hit_frequency.")
    parser.add_argument('-cp','--contaminted_pidsqids', default=False, help="Remove contaminated pidsqids file. Default=False")
    parser.add_argument('-o','--outfile', default='hitdata_plot', help="Output file name prefix.")
    parser.add_argument('-of','--outfmt', default='png', help="format of network plot. Default = png.")
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
