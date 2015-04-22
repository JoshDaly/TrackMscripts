#!/usr/bin/env python
###############################################################################
#
# __KEGG_pipeline_summary_file__.py - description!
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
import urllib
import re

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KEGGBlastFileParser(object):
    def __init__(self, l):
        self.readKeggBlastFile(l)
    
    def readKeggBlastFile(self,l):
        tabs = l.rstrip().split("\t")
        self.uid        = tabs[0]
        self.gene_id    = tabs[1]
        self.identity   = tabs[2]
        self.col4       = tabs[3]
        self.col5       = tabs[4]
        self.col6       = tabs[5]
        self.col7       = tabs[6]
        self.col8       = tabs[7]
        self.col9       = tabs[8]
        self.col10      = tabs[9]
        self.col11      = tabs[10]
        self.col12      = tabs[11]

class KEGGBlastData(object):
    def __init__(self, kegg_blast_file):
        self.kegg_blast_data    = {}
        self.wrapper(kegg_blast_file)
        
    def wrapper(self, kegg_blast_file):
        with open(kegg_blast_file) as fh:
            for l in fh:
                KBFP = KEGGBlastFileParser(l)
                self.kegg_blast_data[KBFP.uid] = [KBFP.gene_id,
                                                  KBFP.identity,
                                                  KBFP.col4,
                                                  KBFP.col5,
                                                  KBFP.col6,
                                                  KBFP.col7,
                                                  KBFP.col8,
                                                  KBFP.col9,
                                                  KBFP.col10,
                                                  KBFP.col11,
                                                  KBFP.col12
                                                  ]

    
class ProkkaTabularData(object):
    def __init__(self, prokka_tabular_file):
        self.prokka_tabular_data    = {}
        self.wrapper(prokka_tabular_file)
    
    def wrapper(self, prokka_tabular_file):
        with open(prokka_tabular_file) as fh:
            for i,l in enumerate(fh):
                if '>' in l:
                    
                    # contains pidsqid
                    pidsqid = l.split()[1]
                    
                    count = 0 
                    
                    # loop through and capture each annotated gene feature 
                    for l in fh:
                        
                        if 'locus_tag' in l:
                            count += 1 
                            # gene id 
                            gene_id         = l.rstrip().split("\t")[4]
                            pidsqid_gene_id = "%s_%d" % (pidsqid, count)
                            
                            self.prokka_tabular_data[gene_id] = pidsqid_gene_id
                            
                        elif '#' in l:
                            break

class KOLookupData(object):
    def __init__(self, ko_lookup_table):
        self.ko_table = {}
        self.grabKOData(ko_lookup_table)
        
    def grabKOData(self, ko_lookup_table):
        with open(ko_lookup_table) as fh:
            for l in fh:
                if len(l) !=0:
                    tabs        = l.rstrip().split("\t")
                    ko          = tabs[0] 
                    definition  = tabs[1]
                    
                    # add to dict
                    self.ko_table[ko] = definition
                
class KGeneLookup(object):
    def __init__(self, k_gene_lookup_table):
        self.k_gene_table = {}
        self.grabKGeneDdata(k_gene_lookup_table)
        
    def grabKGeneDdata(self, k_gene_lookup_table):
        with open(k_gene_lookup_table) as fh:
            for l in fh:
                if len(l) >1:
                    try:
                        tabs    = l.rstrip().split("\t")
                        ko      = tabs[0]
                        kgene   = tabs[1]
                    
                        # add to dict
                        # kgenes can have multiple associated KOs
                        # this sets it to the last seen KO
                        self.k_gene_table[kgene] = ko
                    except IndexError:
                        print len(l)                
    

class KOPipeline(object):
    def __init__(self, kegg_blast_file, prokka_tabular_file, ko_lookup_table, k_gene_lookup_table):
        self.KD     = KEGGBlastData(kegg_blast_file)
        self.PD     = ProkkaTabularData(prokka_tabular_file)
        self.KLD    = KOLookupData(ko_lookup_table)
        self.KGD    = KGeneLookup(k_gene_lookup_table)
        self.kegg_attributes = {}
    
    def wrapper(self, outfile, logfile, parse_log_file, something_entirely_different):
        
        # parse log file
        if parse_log_file:
            self.parseLogFile(parse_log_file, something_entirely_different, outfile)
        
        elif something_entirely_different:
            self.parseKeggData()
        
            self.printKeggData(outfile)
        
        else:
            # search for Kegg description through website
            self.grabKEGGData(logfile)
            
            # print to file
            self.printKeggData(outfile)
        
    def parseKeggData(self):
        for uid in self.KD.kegg_blast_data.keys():
            k_gene          = self.KD.kegg_blast_data[uid][0]
            try:
                ko              = self.KGD.k_gene_table[k_gene]
                ko_definition   = self.KLD.ko_table[ko]
            except KeyError:
                ko              = 'NA'
                ko_definition   = 'NA'
            self.kegg_attributes[uid] = [ko, ko_definition]
    
    def grabKEGGData(self, logfile):
        of = open(logfile, 'w')
        for uid in self.KD.kegg_blast_data.keys():
            gene_id = self.KD.kegg_blast_data[uid][0]
            kegg_search = 'http://www.genome.jp/dbget-bin/www_bget?%s' % gene_id
            of.write("%s STARTED\n" % kegg_search)
            page = urllib.urlopen(kegg_search).read()
            definition, organism  = self.getKeggAttributes(page, uid)
            of.write("%s FINISHED\n" % kegg_search)
            self.kegg_attributes[uid] = [definition, organism]
            of.write("%s\t%s\n" % (uid, self.kegg_attributes[uid]))
            of.flush()
    
    def getKeggAttributes(self, page, geneid):
        definition  = 'NA'
        organism    = 'NA'
        if re.search("No such data", page):
            return definition, organism
        else:
            tmp = page.split("\n")
            for i,l in enumerate(tmp):
                if re.findall('Definition.+', l):
                    # grab line after definition
                    definition=tmp[i+1].rstrip().split(">")[-2][:-3]
                if re.findall('Organism.+', l):
                    organism=tmp[i+1].rstrip().split(">")[-2][12:-3]
                    
            return definition, organism
    
    def printKeggData(self, outfile):
        of = open(outfile,'w')
        
        for gene_id in self.KD.kegg_blast_data.keys():
            string_to_print = '%s' % gene_id
            pidsqid_gene_id = self.PD.prokka_tabular_data[gene_id]
            k_gene          = self.KD.kegg_blast_data[gene_id][0]
            ko              = self.kegg_attributes[gene_id][0]
            ko_definition   = self.kegg_attributes[gene_id][1]
            string_to_print += "\t%s\t%s\t%s\t%s\t" % (pidsqid_gene_id,
                                                       k_gene,
                                                       ko,
                                                       ko_definition)
            string_to_print += "\t".join(self.KD.kegg_blast_data[gene_id][1:])
            of.write("%s\n" % string_to_print)
            
    def parseLogFile(self, log_file, outfile):
        kegg_data = {}
        with open(log_file) as fh:
            for l in fh:
                if 'STARTED' in l:
                    ko_id = l.rstrip().split()[0].split("?")[-1]
                    
                    # continue looping through to grab definition and organism
                    for l in fh:
                        if 'batch' in l:
                            ko_id_attributes    = l.rstrip().split("\t")
                            gene_id             = ko_id_attributes[0]
                            definition          = ko_id_attributes[1].split(",")[0][2:-1]
                            organism            = ko_id_attributes[1].split(",")[1][2:-2]
                            
                            # add data
                            kegg_data[gene_id] = [ko_id,definition,organism]
                            
                            # return to next ko_id
                            break 
        
        of = open(outfile,'w')
        
        # print out log data
        for gene_id in self.KD.kegg_blast_data.keys():
            string_to_print = '%s' % gene_id
            pidsqid_gene_id = self.PD.prokka_tabular_data[gene_id]
            ko_id           = self.KD.kegg_blast_data[gene_id][0]
            definition      = kegg_data[gene_id][1]
            organism        = kegg_data[gene_id][2]
            string_to_print += "\t%s\t%s\t%s\t%s\t" % (pidsqid_gene_id,
                                                     ko_id,
                                                     definition,
                                                     organism)
            string_to_print += "\t".join(self.KD.kegg_blast_data[gene_id][1:])
            of.write("%s\n" % string_to_print)
                            
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
    KOP = KOPipeline(args.kegg_blast_file,
                     args.prokka_tabular_file,
                     args.ko_lookup_table,
                     args.k_gene_lookup_table)
    KOP.wrapper(args.outfile,
                args.logfile,
                args.parse_log_file,
                args.something_entirely_different)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('kegg_blast_file', help="")
    parser.add_argument('prokka_tabular_file', help="")
    parser.add_argument('ko_lookup_table', help="")
    parser.add_argument('k_gene_lookup_table', help="")
    parser.add_argument('-o','--outfile', default='kegg_summary_file.csv',help="")
    parser.add_argument('-l','--logfile', default='KEGG_pipeline_summary_file.log',help="")
    parser.add_argument('-plf','--parse_log_file', default=False,help="Provide log file because Im an idiot")
    parser.add_argument('-sed','--something_entirely_different', default=False,help="")
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
