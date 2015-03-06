#!/usr/bin/env python
###############################################################################
#
# __scrape_SRA_data__.py - description!
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
import urllib
import re
import os
import errno
import shlex
import subprocess
#import numpy as np
#np.seterr(all='raise')

# local imports
import trackm_file_parser as TFP


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ShortReadArchive(object):
    def __init__(self):
        self.project_ids    = {}
        self.SRA_urls       = {}
        self.SRAs           = {}
        self.cmds           = []

    def wrapper(self, project_ids, print_to_file, outdir, num_threads):
        # get project ids
        self.getProjectIds(project_ids)
        
        # grab SRR IDs
        self.grabSRRIds()
        
        # print to file
        if print_to_file:
            self.print_to_file(outdir)
            
        # run wget to grab sra files
        self.downloadSRAFiles(num_threads, outdir)
    
    def getProjectIds(self, project_ids):
        with open(project_ids) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                try:
                    gid         = tabs[0]
                    project_id  = tabs[1]
                    if self.checkIfFileExists(gid):
                        self.project_ids[gid] = project_id
                    else:
                        print gid
                except IndexError:
                    pass
                
    def  grabSRRIds(self):
        count = 0 
        # loop through project ids
        for gid in self.project_ids.keys():
            project_id = self.project_ids[gid]
            ncbi_search = "http://www.ncbi.nlm.nih.gov/sra/?term=%s" % (project_id)
            page = urllib.urlopen(ncbi_search).read()
            # check page
            self.checkSRAPage(page, gid, project_id)
            if count >10:
                break
            #count +=1 
            
    def checkSRAPage(self, page, gid, project_id):
        # check if page contains no SRA record
        if re.search("No items found", page):
            pass
        else:
            # check if project contains single sequence run
            if re.findall("ftp.+",page):
                self.getFTPs(page, gid, project_id)
            else:
                if re.findall("SRX\d{6}",page):
                    self.searchSRXs(page, gid, project_id)
    
    def getFTPs(self, page, gid, project_id):
        FTPs = re.findall("ftp.+",page)
        
        # loop through FTPs
        for ftp in FTPs:
            
            # split ftp on whitespace
            for i in ftp.split():
                
                if re.search("ftp.+/SRR\d.+/SRR\d.+",i):
                    
                    # remove href= & add .sra
                    sra     = i.split("/")[-1][:-1]
                    # add sra
                    self.addSRA(self.SRAs, gid, project_id, sra)
                    
                    sra     = "%s.sra" % sra
                    ftp_url = i.split("=")[-1]
                    ftp_url = "%s/%s" % (ftp_url[1:-1],sra)
                    # add ftp_url
                    self.addSRA(self.SRA_urls, gid, project_id, ftp_url)
                    
                    
    def searchSRXs(self, page, gid, project_id):
        SRX_uniq    = {}
        SRXs        = re.findall("SRX\d{6}",page)
        
        # loop through SRXs
        for srx in SRXs:
            
            # get unique list of SRXs
            SRX_uniq[srx] = 1
        
        # loop through unique list of SRXs
        for srx in SRX_uniq.keys():
            srx_url = "http://www.ncbi.nlm.nih.gov/sra/%s" % srx
            
            srx_page = urllib.urlopen(srx_url).read() 
            
            # search through SRX page
            self.getFTPs(srx_page, gid, project_id)
            
        
    def addSRA(self, dict, gid, project_id, sra):
        try:
            dict[gid][project_id].append(sra)
        except KeyError:
            try:
                dict[gid][project_id] = [sra]
            except KeyError:
                dict[gid] = {project_id:[sra]}
                
    def print_to_file(self,outdir):
        output_file = os.path.join(outdir, 'inter_phyla.gids.project_ids.sras.csv')
        f = open(output_file, 'w')
        for gid in self.SRAs.keys():
            for project_id in self.SRAs[gid]:
                sras = "\t".join(self.SRAs[gid][project_id])
                string_to_print = "\t".join([gid,project_id,sras])
                f.write("%s\n" % string_to_print)
                
    def checkIfFileExists(self, gid):
        file = '/srv/projects/trackm/batch7/inter_phyla_analysis/improved_Taxonomy/transfer_groups/sra/bamm/%s' % gid
        if os.path.isfile(file):
            return False
        else:
            # file is not present 
            return True
    
    def downloadSRAFiles(self, num_threads, outdir):
        # set threads
        pool = Pool(num_threads)
        
        # loop through sra urls
        for gid in self.SRA_urls.keys():
            for project_id in self.SRA_urls[gid]:
                for sra_url in self.SRA_urls[gid][project_id]:
                    sra_file = sra_url.split("/")[-1]
                    outfile = os.path.join(outdir, sra_file)
                    cmd = "wget %s -O %s" % (sra_url, outfile) 
                    self.cmds.append(cmd)
        print pool.map(runCommand, self.cmds)
        
            
        """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""         

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
    SRA = ShortReadArchive()
    SRA.wrapper(args.project_ids,
                args.print_to_file,
                args.outdir,
                args.num_threads)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('project_ids', help="File containing gid linked to project ID eg. A00002472 PRJNA46389")
    parser.add_argument('-p','--print_to_file', default=False, help="")
    parser.add_argument('-o','--outdir', help="")
    parser.add_argument('-t','--num_threads', default=1, type=int, help="")
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
