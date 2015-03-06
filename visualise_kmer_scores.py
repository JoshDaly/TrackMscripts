#!/usr/bin/env python
###############################################################################
#
# __visualise_kmer_scores__.py - Read in kmer score file and visualise on scatter plot!
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

import numpy as np
np.seterr(all='raise')

#import matplotlib as mpl
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
        self.gid_to_taxon   = {}
        self.buildTaxonData(taxon_data)
        
    def buildTaxonData(self, taxon_data):
        with open(taxon_data) as fh:
            for l in fh:
                TFP = TaxonomyFileParser(l)
                TFP.readTaxonFile(l)
                
                # phylum-level taxonomy
                phylum  = TFP.taxonomy.split(";")[1][1:]
                self.taxon_phylum[TFP.gid]  = phylum
                
                # gid to taxon
                self.gid_to_taxon[TFP.gid] = TFP.taxonomy
                
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

class KmerFileParser(object):
    def __init__(self, l):
        self.readKmerFile(l)
        
    def readKmerFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqid    = tabs[0]
        self.gid_1      = tabs[1]
        self.gid_2      = tabs[2]
        self.kmer_score = float(tabs[3])
        self.mean_dg    = float(tabs[4])
        self.dg1        = float(tabs[5])
        self.dg2        = float(tabs[6])
        
class KmerData(object):
    def __init__(self, group_file, taxon_file):
        self.kmer_score       = []
        self.mean_dg          = []
        self.kmer_data        = {}
        self.dgg_data         = {}
        self.closest_dgg_tg   = {}
        self.dg_tg            = {}
        self.dodgey_inter_tgs = [3, 7, 38]
        self.GD               = GroupData(group_file)
        self.TD               = TaxonomyData(taxon_file)

    def buildKmerDataWrapper(self, kmer_file, genus1, genus2):
        with open(kmer_file) as fh:
            for l in fh:
                if l[0] != 'p':
                    if genus1 and genus2: 
                        self.buildKmerDataByGenus(l, genus1, genus2)
                    else:
                        self.buildKmerData(l)
    
    def checkPidsqid(self, pidsqid):
        status = True
        for dodgey_group in self.dodgey_inter_tgs:
            if pidsqid in self.GD.group_data[dodgey_group]:
                status = False
                break
        return status
    
    def buildKmerDataByGenus(self, l, genus1, genus2):
        KFP = KmerFileParser(l)
        KFP.readKmerFile(l)
        
        if self.checkPidsqid(KFP.pidsqid):
            # check if gid equals target genus
            if self.TD.taxon_genus[KFP.gid_1].lower() == genus1.lower():
                if self.TD.taxon_genus[KFP.gid_2].lower() == genus2.lower():
                    # no need to reverse kmer score
                    self.kmer_data[KFP.pidsqid] = [KFP.kmer_score, KFP.mean_dg]
                    print "%s %s %f" % (self.TD.taxon_genus[KFP.gid_1],
                                        self.TD.taxon_genus[KFP.gid_2],
                                        KFP.kmer_score)
                
            elif self.TD.taxon_genus[KFP.gid_1].lower() == genus2.lower(): 
                if self.TD.taxon_genus[KFP.gid_2].lower() == genus1.lower():
                    # must reverse kmer score
                    corr_kmer_score = 1 - KFP.kmer_score
                    print "%s %s %f %f" % (self.TD.taxon_genus[KFP.gid_2],
                                           self.TD.taxon_genus[KFP.gid_1],
                                           KFP.kmer_score,
                                           corr_kmer_score)
                    
                    self.kmer_data[KFP.pidsqid] = [corr_kmer_score, KFP.mean_dg]
        
        #else: 
        #    print "%s %s not targeted" % (self.TD.taxon_genus[KFP.gid_1], self.TD.taxon_genus[KFP.gid_2])
    
    def buildKmerData(self, l):
        KFP = KmerFileParser(l)
        KFP.readKmerFile(l)
        if self.checkPidsqid(KFP.pidsqid):
            self.kmer_data[KFP.pidsqid] = [KFP.kmer_score, KFP.mean_dg]
            self.kmer_score.append(KFP.kmer_score)
            self.mean_dg.append(KFP.mean_dg)
            self.dgg_data[KFP.pidsqid] = [KFP.gid_1, KFP.dg1, KFP.gid_2, KFP.dg2]
                    
    def whichGIDIsCloser(self, gid1, gid2, dg1, dg2, transfer_group, pidsqid):
        dg_contender    = 0
        gid_contender   = ''
        if dg1 < dg2:
            dg_contender = dg1
            gid_contender = gid1
        else:
            dg_contender = dg2
            gid_contender = gid2
        
        try: 
            curr_closest_gid, curr_closet_dg, curr_pidsqid = self.closest_dgg_tg[transfer_group]
            if dg_contender < curr_closet_dg:
                self.closest_dgg_tg[transfer_group] = [gid_contender, dg_contender, pidsqid]  
        except KeyError: 
            self.closest_dgg_tg[transfer_group] = [gid_contender, dg_contender, pidsqid]
            
    def findClosestGID(self):
        """Find closest GID/taxonomy to LGT for each transfer group"""
        # loop through transfer groups
        for transfer_group in self.GD.group_data.keys():
            for pidsqid in self.GD.group_data[transfer_group]:
                gid1 = self.dgg_data[pidsqid][0] 
                gid2 = self.dgg_data[pidsqid][2]
                dg1  = self.dgg_data[pidsqid][1]
                dg2  = self.dgg_data[pidsqid][3]
                try:
                    self.dg_tg[transfer_group].append(dg1)
                    self.dg_tg[transfer_group].append(dg2)
                except KeyError: 
                    self.dg_tg[transfer_group] = [dg1, dg2]
                self.whichGIDIsCloser(gid1, gid2, dg1, dg2, transfer_group, pidsqid)
        self.printOutClosetDgs()
    
    def printOutClosetDgs(self):
        # print header
        print "\t".join(['transfer_group',
                         'dg',
                         'mean_dg',
                         'dg_std',
                         'gid',
                         'pidsqid',
                         'taxon_string'])
        
        for transfer_group in self.closest_dgg_tg.keys():
            dg              = self.closest_dgg_tg[transfer_group][1]
            gid             = self.closest_dgg_tg[transfer_group][0]
            pidsqid         = self.closest_dgg_tg[transfer_group][2]
            taxon_string    = self.TD.gid_to_taxon[gid]
            np_array        = np.array(self.dg_tg[transfer_group])
            mean            = np.mean(np_array)
            stdev           = np.std(np_array)
            print "%s\t%f\t%f\t%f\t%s\t%s\t%s" % (transfer_group,
                                              dg,
                                              mean,
                                              stdev,
                                              gid,
                                              pidsqid,
                                              taxon_string
                                              )
                           
class Plot(object):
    def __init__(self, kmer_data):
        self.KD  = kmer_data
    
    def colourDotsByWrapper(self,transfer_group, pidsqid):
        if transfer_group:
            return self.colourDotsByTransferGroups(transfer_group)
        elif pidsqid:
            return self.colourDotsByPidsqid(pidsqid) 
        else:
            return self.colourDotsByNothin()
            
    def colourDotsByTransferGroups(self, transfer_group):
        # coordinates and colour array
        # transfer group coords
        xs = [] 
        ys = []
        cs = []
        ss = []
        # other coords
        oxs = []
        oys = []
        ocs = []
        oss = []
        
        for pidsqid in self.KD.kmer_data.keys():
            if pidsqid in self.KD.GD.group_data[transfer_group]:
                xs.append(self.KD.kmer_data[pidsqid][0])
                ys.append(self.KD.kmer_data[pidsqid][1])
                ss.append(100)
                cs.append("#FF0000")
            else: 
                oxs.append(self.KD.kmer_data[pidsqid][0]) 
                oys.append(self.KD.kmer_data[pidsqid][1])
                ocs.append("#002BFF") 
                oss.append(100) 

        return xs, ys, cs, ss, oxs, oys, ocs, oss
    
    def colourDotsByPidsqid(self, target_pidsqid):
        # coordinates and colour array
        # pisqid coords
        xs = [] 
        ys = []
        cs = []
        ss = []
        # other coords
        oxs = []
        oys = []
        ocs = []
        oss = []
        
        for pidsqid in self.KD.kmer_data.keys():            
            if pidsqid == target_pidsqid:
                xs.append(self.KD.kmer_data[pidsqid][0])
                ys.append(self.KD.kmer_data[pidsqid][1])
                cs.append("#FF0000")
                ss.append(1000)
            else:
                oxs.append(self.KD.kmer_data[pidsqid][0])
                oys.append(self.KD.kmer_data[pidsqid][1])
                ocs.append("#002BFF")
                oss.append(100)
        return xs, ys, cs, ss, oxs, oys, ocs, oss
        
    def colourDotsByNothin(self):
        # coordinates and colour array
        xs = [] 
        ys = []
        cs = []
        ss = []

        
        for pidsqid in self.KD.kmer_data.keys():
            xs.append(self.KD.kmer_data[pidsqid][0])
            ys.append(self.KD.kmer_data[pidsqid][1])
            cs.append("#002BFF")
            ss.append(100)
        return xs, ys, cs, ss
        
    def scatterPlot(self, outfile, outfmt, pidsqid, transfer_group):
        """build scatter plot"""
        # initialise graph
        fig = plt.figure(figsize=(30,15),dpi=300)
        ax = fig.add_subplot(111, xlim=[0,1], ylim=[0,0.15])
        plt.xticks(np.arange(0,1,0.1))
        plt.xlabel('Kmer score', fontsize=18)
        plt.ylabel('Average dgg', fontsize=18)
        
        # pisqid coords
        xs = [] 
        ys = []
        cs = []
        ss = []
        # other coords
        oxs = []
        oys = []
        ocs = []
        oss = [] 
        
        # colour dots by
        if pidsqid or transfer_group:
            xs, ys, cs, ss, oxs, oys, ocs, oss = self.colourDotsByWrapper(transfer_group, pidsqid)
            ax.scatter(oxs,
                       oys,
                       c=ocs, 
                       alpha=0.7,
                       linewidths=0,
                       s=oss)
            
        else: 
            xs, ys, cs, ss = self.colourDotsByNothin()
            
        ax.scatter(xs,
                   ys,
                   c=cs, 
                   alpha=0.7,
                   linewidths=0,
                   s=ss)
        
        plt.savefig("%s" % (outfile),format="%s" % (outfmt))
        

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
    # build kmer data
    KD = KmerData(args.group_file,
                  args.taxon_file)
    KD.buildKmerDataWrapper(args.kmer_data,
                            args.genus1,
                            args.genus2)
    #KD.findClosestGID()
    
    
    #"""
    # plot scatter
    P = Plot(KD)
    P.scatterPlot(args.outfile,
                  args.outfmt,
                  args.pidsqid,
                  args.transfer_group)
    #"""


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kmer_data', help="File containing kmer data.")
    parser.add_argument('group_file', help="File containing transfer group information.")
    parser.add_argument('taxon_file', help="File containing Phil's improved taxonomy.")
    parser.add_argument('-o','--outfile', default = 'scatter.png',help="Output file.")
    parser.add_argument('-pidsqid','--pidsqid', default=False,help="Please specify which pidsqid you would like to highlight.")
    parser.add_argument('-tg','--transfer_group', default=False,type=int,help="Please specify which transfer group you would like to highlight.")
    parser.add_argument('-outfmt','--outfmt', default = 'png',help="Output file format.")
    parser.add_argument('-g1','--genus1', default = False,help="Select target genus 1.")
    parser.add_argument('-g2','--genus2', default = False,help="Select target genus 2.")
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
