#!/usr/bin/env python
###############################################################################
#
# __run_kmer_analysis__.py - description!
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
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import os
import errno

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

class PathsFileParser(object):
    def __init__(self,l):
        self.readPathsFile(l)
    
    def readPathsFile(self,l):
        tabs = l.rstrip().split("\t")
        self.gid             = tabs[0]
        self.path_to_file    = tabs[1]
        
class PathsFileData(object):
    def __init__(self, paths_file):
        self.gid_to_file = {}
        self.buildPathsData(paths_file)
    
    def buildPathsData(self, paths_file):
        with open(paths_file) as fh:
            for l in fh:
                PFP = PathsFileParser(l)
                if l[0] != "#":
                    PFP.readPathsFile(l)
                    self.gid_to_file[PFP.gid] = PFP.path_to_file

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
    def __init__(self, paths_file, transfer_group_file):
        self.PFD            = PathsFileData(paths_file)
        self.GD             = GroupData(transfer_group_file)
        self.kmer_cmds      = []
        self.lgt_fasta_data = {}
        
    def checkGroupMembership(self,pidsqid):
        try:
            group_member = self.GD.group_membership[pidsqid]
            return group_member
        except KeyError:
            return False
  
    def makeKmerCommand(self, pidsqid, group_member, genome_path1, genome_path2, data):
        """
        Runs these command for each lgt event
        kmer_counter.rb -w 500 -W 504 -m 500 genome1.fasta >genome1.kmer_counts.csv
        kmer_counter.rb -w 500 -W 504 -m 500 genome2.fasta >genome2.kmer_counts.csv
        kmer_counter.rb -w 500 -W 504 -m 500 lgt.fasta >lgt.kmer_counts.csv
        """
        
        if data == "inter_phyla":
            path_to_pidsqid = "/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/transfer_groups/evalue0/TG_%s/%s/%s.fna" % (group_member, pidsqid, pidsqid)
            path_to_TG      = "/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/transfer_groups/evalue0/TG_%s/%s/" % (group_member, pidsqid)
        
        elif data == "all":
            path_to_pidsqid = "/srv/projects/trackm/batch7/all_transfers/kmer_analysis/TG_%s/%s/%s.fna" % (group_member, pidsqid, pidsqid)
            path_to_TG      = "/srv/projects/trackm/batch7/all_transfers/kmer_analysis/TG_%s/%s/" % (group_member, pidsqid)
        
        # input files
        genome1_name    = genome_path1.rstrip().split("/")[-1].split('.')[0]
        genome2_name    = genome_path2.rstrip().split("/")[-1].split('.')[0]
        
        # output files
        genome1_kmer    = os.path.join(path_to_TG,"genome_%s.kmer_counts.csv" % (genome1_name))
        genome2_kmer    = os.path.join(path_to_TG,"genome_%s.kmer_counts.csv" % (genome2_name))
        lgt_kmer        = os.path.join(path_to_TG,"%s.kmer_score.csv" % (pidsqid))
        
        # genome 1 
        self.kmer_cmds.append(["kmer_counter.rb -w 500 -W 504 -m 500 %s > %s" % (genome_path1,genome1_kmer)])
        # genome 2 
        self.kmer_cmds.append(["kmer_counter.rb -w 500 -W 504 -m 500 %s > %s" % (genome_path2,genome2_kmer)])
        # lgt seq
        self.kmer_cmds.append(["kmer_counter.rb -w 500 -W 504 -m 500 %s > %s" % (path_to_pidsqid,lgt_kmer)])
    
    def placeLGTFastaInDirectory(self, lgt_fasta, data):
        # build data structure containing all fasta sequences        
        # parse fasta file using biopython
        for accession,sequence in SeqIO.to_dict(SeqIO.parse(lgt_fasta,"fasta")).items():
            self.lgt_fasta_data[accession] = sequence.seq
        
        # loop though transfer groups
        for TG in self.GD.group_data.keys():
            for pidsqid in self.GD.group_data[TG]:
                
                if data == "inter_phyla":
                    # create file to write to
                    f = open('/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/transfer_groups/evalue0/TG_%s/%s/%s.fna' % (TG, pidsqid, pidsqid),'w')
                    
                    f.write('>%s\n' % pidsqid)
                    f.write(str(self.lgt_fasta_data[pidsqid]))
                    f.close()
                    
                elif data == "all":
                    # create file to write to
                    f = open('/srv/projects/trackm/batch7/all_transfers/kmer_analysis/TG_%s/%s/%s.fna' % (TG, pidsqid, pidsqid),'w')
                    
                    f.write('>%s\n' % pidsqid)
                    f.write(str(self.lgt_fasta_data[pidsqid]))
                    f.close()
                    
                else: 
                    print " "
                    print "########################################"
                    print "ERROR: Please select --data inter_phyla or all"
                    print "########################################"
                    print " "
                    sys.exit()
                
    def mainWrapper(self, hit_file, data, num_threads):
        with open(hit_file) as fh:
            for l in fh: 
                if l[0:3] != 'hid':
                    
                    # parse hit data
                    HDP = HitDataParser(l)
                    HDP.readHitData(l)
                    
                    # grab pid_sqid 
                    pid_sqid1 = "%s_%s" % (HDP.pid, HDP.sqid_1)
                    pid_sqid2 = "%s_%s" % (HDP.pid, HDP.sqid_2)
                    
                    # grab genome paths
                    genome_path1 = self.PFD.gid_to_file[HDP.gid_1] 
                    genome_path2 = self.PFD.gid_to_file[HDP.gid_2]
                    
                    # check group membership, and build data
                    group_member = self.checkGroupMembership(pid_sqid1)
                    
                    if group_member:
                        
                        # build commands
                        self.makeKmerCommand(pid_sqid1,
                                             group_member,
                                             genome_path1,
                                             genome_path2,
                                             data)
                        
        # set the number of threads to use
        pool = Pool(num_threads)
        
        #print self.kmer_cmds 
                   
        # run commands
        print pool.map(runCommand, self.kmer_cmds) # run analysis
                        
    def buildDirectoryStructure(self, num_threads, data):
        # build array of cmds
        # of the form: 
        # /srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/TG_#/pidsqid
        cmds = [] 
        
        # loop though transfer groups
        for TG in self.GD.group_data.keys():
            
            # need to create TG_# directory first
            if data == "inter_phyla":
                cmds.append(["mkdir /srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/transfer_groups/TG_%s/" % (TG)])
                
                for pidsqid in self.GD.group_data[TG]:
                    cmds.append(["mkdir /srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/kmer_analysis/transfer_groups/TG_%s/%s" % (TG, pidsqid)])
                
            elif data == "all":
                cmds.append(["mkdir /srv/projects/trackm/batch7/all_transfers/kmer_analysis/transfer_groups/TG_%s" % (TG)])
                            
                for pidsqid in self.GD.group_data[TG]:
                    cmds.append(["mkdir /srv/projects/trackm/batch7/all_transfers/kmer_analysis/transfer_groups/TG_%s/%s" % (TG, pidsqid)])

            else:
                print " "
                print "########################################"
                print "ERROR: Please select --data inter_phyla or all"
                print "########################################"
                print " "
                sys.exit()

        # run commands
        for cmd in cmds:
            runCommand(cmd)

def runCommand(cmd):
        """Run a command and take care of stdout
    
        expects 'cmd' to be a string like "foo -b ar"
    
        returns (stdout, stderr)
        
        Must be outside class object!!!!!!!
        """
        print cmd
        p = subprocess.Popen(cmd, shell=True)
        return p.communicate()  
         
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""
    
    HD = HitData(args.paths_file,
                 args.transfer_group_file)
    if args.task.lower() == 'run':
        HD.mainWrapper(args.hit_data,
                       args.data,
                       args.num_threads)
    elif args.task.lower() == 'build':
        HD.buildDirectoryStructure(args.num_threads,
                                   args.data)
    
    elif args.task.lower() == 'place':
        HD.placeLGTFastaInDirectory(args.lgt_fasta,
                                    args.data)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hit_data', help="File containing hit data")
    parser.add_argument('paths_file', help="File containing gid and path to genome file")
    parser.add_argument('transfer_group_file', help="File containing transfer groups")
    parser.add_argument('task', help="Set task for program: run, build, or place.")
    parser.add_argument('-t','--num_threads', type=int, default = 5, help="Number of threads for multiprocessing")
    parser.add_argument('-d','--data', help="Set whether inter phyla or all data is to be used.")
    parser.add_argument('-fna','--lgt_fasta', help="File containing all lgt sequences in fasta format.")
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
