#!/usr/bin/env python
###############################################################################
#
# __COG_bar_chart__.py - description!
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

import os
import errno

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

# group file
class GroupFileParser(object):
    def __init__(self,l):
        self.readGroupFile(l)
        
    def readGroupFile(self,l):
        tabs = l.rstrip().split("\t")
        self.group_num = int(tabs[0].split(':')[0].split()[1])
        self.pid_sqids = tabs[2].split(',')

# COG annotation file 
class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnoFile(l)
        
    def readAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqidgene    = tabs[0] 
        self.pid            = tabs[0].split("_")[0] 
        self.sqid           = tabs[0].split("_")[1]
        self.gene           = tabs[0].split("_")[2]
        self.COG            = tabs[5]

# Kegg annotation file
class KeggAnnotationFileParser(object):
    def __init__(self,l):
        self.readKeggAnnoFile(l)
    
    def readKeggAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqidgene    = tabs[0] 
        self.pid            = tabs[0].split("_")[0] 
        self.sqid           = tabs[0].split("_")[1]
        self.gene           = tabs[0].split("_")[2]
        self.geneid         = tabs[1]

# Kegg master file
class KeggDBParser(object):
    def __init__(self,l):
        self.readKeggDB(l)
        
    def readKeggDB(self,l):
        tabs = l.rstrip().split("\t")
        self.geneID = tabs[0]
        self.KO     = tabs[1] 
        
class COGCatFileParser(object):
    def __init__(self,l):
        self.readCOGCatFile(l)
        
    def readCOGCatFile(self,l):
        commas = l.rstrip().split(',')
        self.cogID      = commas[0] 
        self.catergory  = commas[1]

# Annotation data
class AnnotationData(object):
    def __init__(self):
        # TID = Transfer ID
        self.groupsToTID    = {}
        self.keggLookUp     = {}
        self.keggToTID      = {}
        self.cogToTID       = {}
        self.cogCounts      = {}
        self.keggCounts     = {}
        self.cogCategories  = {}
        
    def buildTIDGroups(self,file):
        with open(file,'r') as fh:
            for l in fh:
                GFP = GroupFileParser(l)
                GFP.readGroupFile(l)
                try:
                    self.groupsToTID[GFP.group_num] += GFP.pid_sqids 
                except KeyError:
                    self.groupsToTID[GFP.group_num] = GFP.pid_sqids
    
    def buildKeggDB(self,file):
        with open(file,'r') as fh:
            for l in fh:
                KDBP = KeggDBParser(l)
                KDBP.readKeggDB(l)
                self.keggLookUp[KDBP.KO] = KDBP.geneID
                
    def buildCOGCats(self,file):
        with open(file,'r') as fh:
            for l in fh:
                CCFP = COGCatFileParser(l)
                CCFP.readCOGCatFile(l)
                self.cogCategories[CCFP.cogID] = CCFP.catergory
                
    def getCOGData(self,file):
        with open(file,'r') as fh:
            for l in fh:
                AFP = AnnotationFileParser(l)
                AFP.readAnnoFile(l)
                
                if len(AFP.COG) > 1 :
                    self.cogToTID[AFP.pidsqidgene] = AFP.COG
                else:    
                    self.cogToTID[AFP.pidsqidgene] = 'unknown'
                    
    def getKeggData(self,file):
        with open(file,'r') as fh:
            for l in fh:
                KAFP = KeggAnnotationFileParser(l)
                KAFP.readKeggAnnoFile(l)
                self.keggToTID[KAFP.pidsqidgene] = KAFP.geneid

    def barChart(self, group_num, dbtype, outdir):
        # reset dictiaries
        self.cogCounts      = {}
        self.keggCounts     = {}
        
        # coordinates
        x = []
        y = []
        
        # 
        n_groups = 0
        
        ##############################
        # COG bar chart
        ##############################
        
        if dbtype.lower() == 'cog': 
            # reset dictiaries
            self.cogCounts      = {}
            
            # build x and y 
            for pidsqidgene in self.cogToTID.keys():
                pidsqid = '_'.join(pidsqidgene.split('_')[0:2])
                
                if pidsqid in self.groupsToTID[group_num]:
                    
                    if self.cogToTID[pidsqidgene] == 'unknown':
                        try:
                            self.cogCounts['unknown'] += 1
                        except KeyError:
                            self.cogCounts['unknown'] = 1
                    else:
                        if len(self.cogCategories[self.cogToTID[pidsqidgene]]) > 1: 
                            for cog in self.cogCategories[self.cogToTID[pidsqidgene]]: 
                                try:
                                    self.cogCounts[cog] += 1
                                except:
                                    self.cogCounts[cog] = 1
                        else:
                            try: 
                                self.cogCounts[self.cogCategories[self.cogToTID[pidsqidgene]]] += 1
                            except KeyError: 
                                self.cogCounts[self.cogCategories[self.cogToTID[pidsqidgene]]] = 1
            
            # sort cog categories        
            sorted_cog_cats = self.cogCounts.keys()
            sorted_cog_cats.sort()
            
            n_groups = len(sorted_cog_cats)
            
            for cog_cat in sorted_cog_cats:
                x.append(self.cogCounts[cog_cat])
                y.append(cog_cat)
        
        ##############################
        # KEGG bar chart 
        ##############################
        
        elif dbtype.lower() == 'kegg': 
            # reset dictionaries
            self.keggCounts     = {}
            present = False
            
            
            # build x and y 
            for pidsqidgene in self.keggToTID.keys():
                pidsqid = '_'.join(pidsqidgene.split('_')[0:2])
                
                
                if pidsqid in self.groupsToTID[group_num]:
                    #print "Group %d: %s" % (group_num, pidsqid)
                    # search through Kegg DB
                    for KO in self.keggLookUp.keys():
                        
                        # check if gene in kegg DB
                        if self.keggToTID[pidsqidgene] in self.keggLookUp[KO]:
                            present = True
                            try:
                                self.keggCounts[KO] += 1 
                            except KeyError:
                                self.keggCounts[KO] = 1
                                 
                    # gene not present in Kegg DB, set to unknown
                    if not present:
                        try:
                            self.keggCounts['unknown'] += 1 
                        except KeyError:
                            self.keggCounts['unknown'] = 1
                # reset 
                present = False 
            
            # sort cog categories        
            sorted_kegg_cats = self.keggCounts.keys()
            sorted_kegg_cats.sort()
            
            n_groups = len(sorted_kegg_cats)
            
            # populate coordinates
            for KO in sorted_kegg_cats:
                x.append(self.keggCounts[KO])
                y.append(KO)



        fig = plt.figure()
        ax = fig.add_subplot(111) 
        
        #necessary variables
        index = np.arange(n_groups)
        width = 0.35
        
        # the bars
        plt.bar(index,
                x,
                width)
        
        # axes and labels
        ax.set_xlim(-width,len(index)+width)
        xTickMarks = y
        ax.set_xticks(index+(width/2))
        xtickNames = ax.set_xticklabels(xTickMarks)
        plt.setp(xtickNames, fontsize=10)
        plt.tight_layout()
        
        if dbtype.lower() == 'cog':
            # set title
            ax.set_title('Group %d COG groups' % (group_num))
            
            # output directory/file
            outfile = os.path.join(outdir,"Group_%d.cog.png" % (group_num))
            plt.savefig(outfile,format='png')
        elif dbtype.lower() == 'kegg':
            # set title
            ax.set_title('Group %d KEGG groups' % (group_num))
            
            # output directory/file
            outfile = os.path.join(outdir,"Group_%d.kegg.png" % (group_num))
            plt.savefig(outfile,format='png')
        
        
    def buildBarCharts(self, dbtype, outdir):
        for group in self.groupsToTID.keys():
            self.barChart(int(group),dbtype, outdir)
    
    def testBuildBarCharts(self):
        count = 0 
        
        for group in self.groupsToTID.keys():
            self.barChart(int(group),'kegg')
            
            count += 1 
            if count >= 10:
                break
        
        
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
    # objects 
    AD = AnnotationData()
    
    # Kegg Look Up
    AD.buildKeggDB(args.kegg_lookup)
                
    # Transfer groups
    AD.buildTIDGroups(args.group_file)
    
    # COG categories
    AD.buildCOGCats(args.cog_cats)
    
    # COG data
    AD.getCOGData(args.cog_file)
    
    # Kegg data
    AD.getKeggData(args.kegg_file)
    
    # create bar charts 
    AD.buildBarCharts(args.dbtype, args.outdir)
    #AD.testBuildBarCharts()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kegg_lookup', help="")
    parser.add_argument('group_file', help="")
    parser.add_argument('cog_file', help="")
    parser.add_argument('cog_cats', help="")
    parser.add_argument('kegg_file', help="")
    parser.add_argument('-dbtype','--dbtype', default='cog',help="Set Database type. kegg or cog.")
    parser.add_argument('-o','--outdir', default='cog',help="Set output directory.")
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
