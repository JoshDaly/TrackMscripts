#!/usr/bin/env python
###############################################################################
#
# __run_bamm_multiplexed__.py - description!
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
import subprocess
import glob
import os
import errno
import shlex

# local imports
import trackm_file_parser as TFP

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamMRunner(object):
    def __init__(self, paths_file):
        self.Path               = TFP.PathsFileData(paths_file)
        self.bamm_make_cmds     = []
        self.bamm_parse_cmds    = []
        self.genome_sras        = {}
        self.working_gids       = {}
        self.genomes_to_update  = {}
        
    def wrapper(self, sra_master_list, genomes_dir, bamm_dir, process_to_run, genomes_to_update):
        
        # get list of genomes to run bamm make 
        self.readGenomesToUpdate(genomes_to_update)
        
        # read inter phyla metadata file containing sras
        self.readSRAMasterList(sra_master_list)
        
        if process_to_run.lower() == 'all':
            # make bamm make commands
            self.bammMake(genomes_dir,
                          bamm_dir,
                          genomes_to_update)
        
            # run bamm make 
            for cmd in self.bamm_make_cmds:
                runCommand(cmd)
        
            # make bamm parse commands
            self.bammParse()
    
            # run bamm parse
            for cmd in self.bamm_parse_cmds:
                runCommand(cmd)
                    
        elif process_to_run.lower() == 'make':
            # make bamm make commands
            self.bammMake(genomes_dir,
                          bamm_dir,
                          genomes_to_update)
        
            # run bamm make 
            for cmd in self.bamm_make_cmds:
                runCommand(cmd)
        
        elif process_to_run.lower() == 'parse':
            # make bamm parse commands
            self.bammParse()
    
            # run bamm parse
            for cmd in self.bamm_parse_cmds:
                runCommand(cmd)

    def readGenomesToUpdate(self, genomes_to_update):
        try:
            with open(genomes_to_update) as fh:
                for l in fh:
                    gid = l.rstrip()
                    self.genomes_to_update[gid] = 1
        except TypeError:
            print "*********************************************************************************"
            print "                             Updating all genomes                                "
            print "*********************************************************************************"
    
    def readSRAMasterList(self, sra_master_list):
        with open(sra_master_list) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                gid = tabs[0]
                for tab in tabs[1:]:
                    sra = tab
                    try:
                        self.genome_sras[gid][sra] = 1 
                    except KeyError:
                        self.genome_sras[gid] = {sra:1}

    
    def bammMake(self, genomes_dir, bamm_dir, genomes_to_update):
        # bamm make -d genomes/6 -c  -t 20 -o ./bamm/A00001904/ #
        for gid in self.genome_sras.keys():
            if self.checkGenomesToUpdate(gid, genomes_to_update):
                if genomes_dir[-1] == '/':
                    genomes_dir = genomes_dir[:-1]
                    
                if bamm_dir[-1] == '/':
                    bamm_dir = bamm_dir[:-1]
    
                img_id      = self.Path.gid_to_img[gid]
                sra_paired, sra_single = self.matchSRAwithGID(gid)
                
                # create BamM make commands
                self.makeBammCmds(genomes_dir, img_id, sra_paired, sra_single, bamm_dir, gid) 
            
    def checkGenomesToUpdate(self, gid, genomes_to_update):
        if genomes_to_update:
            try:
                lenny = self.genomes_to_update[gid]
                return True
            except KeyError: 
                return False
        else:
            return True
    def makeBammCmds(self, genomes_dir, img_id, sra_paired, sra_single, bamm_dir, gid):
        if len(sra_paired) > 0 and len(sra_single) > 0:
            self.bamm_make_cmds.append("bamm make -f -d %s/%s.fna -c %s -s %s -t 20 -o %s/%s" % (genomes_dir,
                                                                                                 img_id,
                                                                                                 sra_paired,
                                                                                                 sra_single,
                                                                                                 bamm_dir,
                                                                                                 gid))
        elif len(sra_paired) > 0 and len(sra_single) == 0:
            self.bamm_make_cmds.append("bamm make -f -d %s/%s.fna -c %s -t 20 -o %s/%s" % (genomes_dir,
                                                                                           img_id,
                                                                                           sra_paired,
                                                                                           bamm_dir,
                                                                                           gid))
        elif len(sra_paired) == 0 and len(sra_single) >0:
            self.bamm_make_cmds.append("bamm make -f -d %s/%s.fna -s %s -t 20 -o %s/%s" % (genomes_dir,
                                                                                           img_id,
                                                                                           sra_single,
                                                                                           bamm_dir,
                                                                                           gid))
        else:
            pass
        
    def matchSRAwithGID(self, gid):
        sra_paired = ''
        sra_single = ''
        sra1 = ''
        sra2 = ''
        count = 1
        sra_files = glob.glob('/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/transfer_groups/sra/*.fastq.gz')
        
        for sra_file in sra_files:
            for sra in self.genome_sras[gid]:
                if sra in sra_file:
                    if '_1.fastq.gz' in sra_file or '_2.fastq.gz' in sra_file:
                        sra_paired += '%s ' % sra_file
                    
                    else:
                        sra_single += '%s ' % sra_file
        return sra_paired, sra_single
    
    def bammParse(self):
        for gid in self.genome_sras.keys():
            bamm_dir    = '/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/transfer_groups/sra/bamm/%s' % gid
            try:
                bamm_files = glob.glob('%s/*.bam' % bamm_dir)
                #bamm_file  = glob.glob('%s/*.bam' % bamm_dir)[0]
                bamm_output, threads = self.getBamMOutputFiles(bamm_files)
                if threads >0:
                    self.bamm_parse_cmds.append('bamm parse -l %s/%s.links.tsv -c %s/%s.coverages.tsv -i %s/%s.inserts.tsv -m pmean -b %s -t %d' % (bamm_dir,
                                                                                                                                                    gid,
                                                                                                                                                    bamm_dir,
                                                                                                                                                    gid,
                                                                                                                                                    bamm_dir,
                                                                                                                                                    gid,
                                                                                                                                                    bamm_output,
                                                                                                                                                    threads
                                                                                                                                                    ))   
            except IndexError:
                pass
            
    def getBamMOutputFiles(self, bamm_files):
        threads = 0
        bamm_output = ''
        try:
            for bamm_file in bamm_files:
                threads += 1
                bamm_output += '%s ' % bamm_file
            return bamm_output, threads
        except IndexError:
            pass
        
                
            
            
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
        """Run a command and take care of stdout
    
        expects 'cmd' to be a string like "foo -b ar"
    
        returns (stdout, stderr)
        
        Must be outside class object!!!!!!!
        """
        print cmd
        args = shlex.split(cmd)
        p = subprocess.Popen(args)
        return p.communicate() 
        
def doWork( args ):
    """ Main wrapper"""
    BR = BamMRunner(args.paths_file)
    BR.wrapper(args.sra_master_list,
               args.genomes_dir,
               args.bamm_dir,
               args.process_to_run,
               args.genomes_to_update
               )
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('paths_file', help="")
    parser.add_argument('sra_master_list', help="")
    parser.add_argument('genomes_dir', help="")
    parser.add_argument('bamm_dir', help="")
    parser.add_argument('-ptr','--process_to_run', default='all',help="Parse, make or all")
    parser.add_argument('-gtu','--genomes_to_update', help="")
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
