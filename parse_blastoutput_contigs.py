#!/usr/bin/env python
###############################################################################
#
# __parse_blastoutput_contigs__.py - description!
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
import glob
import os
import errno
import matplotlib.pyplot as plt

# local imports
import trackm_file_parser as TFP
from cb2cols import Cb2Cols as CB2
import grab_gg_taxonomy as GGT

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Contigs(object):
    def __init__(self, taxonomy_file, group_file, hitdata, gg_tax_file):
        self.TD                     = TFP.TaxonomyData(taxonomy_file)
        self.GD                     = TFP.GroupData(group_file)
        self.HD                     = TFP.HitData(hitdata)
        self.GGT                    = GGT.GreenGenesTaxonomy(gg_tax_file)
        self.cb2                    = CB2()
        self.blast_data             = {}
        self.isVector               = {}
        self.contig_len             = {}
        self.count                  = 0
        self.transfer_groups        = {} 
        self.contam_contigs         = {}
        self.contam_contig_gids     = {}
        self.bamm_data              = {}
        self.bamm_cov_links_files   = {}
        self.pidsqid_count          = {}
        self.pidsqid_list           = {}
        self.ecoli                  = self.getEcoliOuttaHere()
    
    def getEcoliOuttaHere(self):
        ecoli_gids  =   []
        for gid in self.HD.gid_to_genus.keys():
            if self.HD.gid_to_genus[gid] == 'Escherichia':
                ecoli_gids.append(gid)
        return ecoli_gids
    
    def wrapper(self, blastdir, cov_links_dir, type, tg_status_outfile, contaminated_contigs_file, length_min, length_max, coverage_min, coverage_max, show_plot, pie_graph_outfile, pie_graph_outfmt, pidsqid_outfile):
        # get contaminated contigs
        self.getContaminatedContigs(contaminated_contigs_file)
        
        # grab bamm stats files
        self.getCovLinkFiles(cov_links_dir)
        self.getCovLinkData()
        
        # read in blast dir, grab blast files
        self.readBLASTDir(blastdir)
        
        if type.lower() == 'pie_chart':
            # create pie chart
            self.plotPieChart(length_min, length_max, coverage_min, coverage_max, pie_graph_outfile, pie_graph_outfmt, show_plot, pidsqid_outfile)
        
        elif type.lower() == 'tg_status':
            self.returnTGStatus(tg_status_outfile)
            
        elif type.lower() == 'list':
            self.returnListOfPidsqids()
        
    def getCovLinkFiles(self, cov_link_dir):
        if cov_link_dir[-1] == '/':
            # remove last element
            cov_link_dir = cov_link_dir[0:-1]
            
        cov_link_files = glob.glob("%s/*/*cov_links_stats.csv" % cov_link_dir)
        for file in cov_link_files:
            gid = file.split("/")[-2]
            if gid in self.contam_contigs:
                self.bamm_cov_links_files[gid] = file
            
    def getCovLinkData(self):
        for gid in self.bamm_cov_links_files.keys():
            if gid not in self.ecoli:
                try:
                    cov_link_file = self.bamm_cov_links_files[gid]
                    CLD = TFP.CovLinksData(cov_link_file)
                    for contig in CLD.cov_link_data[gid]:
                        contig_len                      = int(CLD.cov_link_data[gid][contig][0])
                        coverage                        = float(CLD.cov_link_data[gid][contig][1])
                        try:
                            link_data                       = int(CLD.cov_link_data[gid][contig][2])
                            links_relative_to_contig_length = link_data/float(contig_len)
                            links_rel_to_largest_contig     = float(CLD.cov_link_data[gid][contig][3])
                            self.bamm_data[gid][contig] = [contig_len, coverage, link_data, links_relative_to_contig_length, links_rel_to_largest_contig]
                        except (IndexError, KeyError):
                            try:
                                link_data                       = int(CLD.cov_link_data[gid][contig][2])
                                links_relative_to_contig_length = link_data/float(contig_len)
                                links_rel_to_largest_contig     = float(CLD.cov_link_data[gid][contig][3])
                                self.bamm_data[gid] = {contig:[contig_len, coverage, link_data, links_relative_to_contig_length, links_rel_to_largest_contig]}
                            except (IndexError, KeyError):
                                try:
                                    self.bamm_data[gid][contig] = [contig_len, coverage]
                                except KeyError:
                                    self.bamm_data[gid] = {contig:[contig_len, coverage]}
                except KeyError:
                    pass
    
    def getContaminatedContigs(self, contaminated_contigs_file):
        with open(contaminated_contigs_file) as fh:
            for l in fh:
                tabs    = l.rstrip().split("\t")
                gid     = tabs[0]
                contig  = tabs[1]
                self.contam_contig_gids[contig] = gid
                try:
                    self.contam_contigs[gid][contig] = 1 
                except KeyError: 
                    self.contam_contigs[gid] = {contig:1} 
                
                    
    
    def readBLASTDir(self, blastdir):
        blast_search = os.path.join(blastdir,"*/*vs_nr.blast20151802")
        blast_files = glob.glob(blast_search)
        
        count = 0
        
        for blast_file in blast_files:
            gid = blast_file.split("/")[-2]
            
            #count += 1 
            
            # check gid
            if gid in self.contam_contigs:
            
                # remove E.coli from analysis
                #if self.TD.taxon_genus[gid] != 'Escherichia':
                BFP = TFP.BlastFileParser(blast_file)
                
                for query_id in BFP.isvector.keys():
                    contig_id = self.getContigID(query_id)
                    if contig_id in self.contam_contigs[gid]:
                        self.isVector[contig_id]     = BFP.isvector[query_id]
                
                for query_id in BFP.topHits2.keys():
                    contig_id = self.getContigID(query_id)
                    if contig_id in self.contam_contigs[gid]:
                        self.blast_data[contig_id]   = BFP.topHits2[query_id]
                    
                for query_id in BFP.contig_len.keys():
                    contig_id = self.getContigID(query_id)
                    if contig_id in self.contam_contigs[gid]:
                        self.contig_len[contig_id]   = BFP.contig_len[query_id]
    
            if count > 2000:
                break
    
    def getContigID(self, query_id):
        underscore = query_id.split("_")
        contigid   = "_".join(underscore[1:])
        return contigid
                    
    def plotPieChart(self, length_min, length_max, coverage_min, coverage_max, pie_graph_outfile, pie_graph_outfmt, show_plot, pidsqid_outfile):
        # get pie chart data
        labels,sizes = self.getPieChartData(length_min, length_max, coverage_min, coverage_max)
        
        col_set = "qualSet3"
        ColBrewColours = self.cb2.maps[col_set].values()[0:len(labels)+1]
        
        # initialise outfile
        pidsqid_outfile = open(pidsqid_outfile, 'w')
        
        # print pidsqid counts data
        pos1 = -1.7
        pos2 = 0
        for category in self.pidsqid_count.keys():
            pos2 = pos2 + 0.1
            plt.text(pos1, pos2, "%s: %d" % (category, self.pidsqid_count[category]))
            try:
                pidsqid_outfile.write("%s\t%s\n" % (category, self.pidsqid_list[category]))
            except KeyError:
                pass
        # print header
        pos2 = pos2 + 0.1
        plt.text(pos1, pos2, "Pidsqids")
        
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=ColBrewColours)
        
        plt.axis('equal')
        
        if show_plot:
            plt.show()
        else:
            outfile = "%s.%s" % (pie_graph_outfile, pie_graph_outfmt)
            plt.savefig(outfile, format=pie_graph_outfmt)
            
        
        
        
    def getPieChartData(self, length_min, length_max, coverage_min, coverage_max):
        common_contams = ['pseudomonas', 'propionibacterium', 'Streptococcus', 'Micrococcus', 'Stenotrophomonas', 'Xanthomonas', 'Pseudoxanthomonas', 'Burkholderia', 'Deinococcus', 'Corynebacterium']
        
        # contig properties
        len_min = length_min
        len_max = length_max
        cov_min = coverage_min
        cov_max = coverage_max
        # function attributes
        piechartdata            = {}
        piechartdataCorrected   = {}
        total                   = 0
        labels                  = []
        sizes                   = []
        
        # create temp file for inter phyla output
        temp_output = open('temp_output_file','w')
        
        for contig in self.blast_data.keys():
            
            # check contig satisfies length minimum
            if self.contig_len[contig] >= 500:
                
                # grab gid for contig
                contig_gid  = self.contam_contig_gids[contig]
                    
                # check contig satisfies length and coverage limits
                if self.doesContigSatisfyConstraints(contig_gid, contig, len_min, len_max, cov_min, cov_max):
                    
                    total += 1 
                    
                    # get contig blast information
                    description = self.blast_data[contig][0]
                    evalue      = self.blast_data[contig][1]
                    
                    if  self.isVector[contig] > 0:
                        description = 'Vector'
                        self.addData(description, piechartdata)
                        self.addContigPidsqidCount(contig, description)
                    else:
                        # check top blast hit
                        description = self.checkBLAST(temp_output, contig, contig_gid, description)
                        self.addData(description, piechartdata)
                        self.addContigPidsqidCount(contig, description)
        
        print "Total contigs: %d" % (total)
         
        piechart_labels = ['Other','Vector','Inter Phyla', 'Genus', 'Family', 'Order', 'Class', 'Phylum']
        for description in  piechart_labels:
            count = piechartdata[description]
            labels.append(description)
            sizes.append(count)
        
        return labels, sizes       
    
    def addContigPidsqidCount(self, contig, category):
        try:
            pidsqid_count   = len(self.HD.contig_pidsqid_count_all[contig])
            pidsqid_list = ",".join(self.HD.contig_pidsqid_count_all[contig].keys())
            try:
                self.pidsqid_list[category] += pidsqid_list
            except KeyError:
                self.pidsqid_list[category] = pidsqid_list
        except KeyError:
            pidsqid_count = 0 
        
        try:
            self.pidsqid_count[category] += pidsqid_count
        except KeyError:
            self.pidsqid_count[category] = pidsqid_count
            
        
    
    def doesContigSatisfyConstraints(self, gid, contig, len_min, len_max, cov_min, cov_max):
        # get contig attributes
        contig_len = int(self.bamm_data[gid][contig][0])
        contig_cov = float(self.bamm_data[gid][contig][1])
        
        # check against constraints
        if contig_len >= len_min and contig_len <= len_max and contig_cov <= cov_max and contig_cov >= cov_min:
            return True
        else:
            #print "Contig: %s Length: %d Coverage: %f" % (contig, contig_len, contig_cov)
            return False
             
    
    def checkBLAST(self, temp_output, contig, contig_gid, description):
        # split description on whitespace
        description_split = description.split()
        
        original_description = description
        
        # get contig taxon information
        _organism   = self.TD.taxon_organism[contig_gid]
        _genus      = self.TD.taxon_genus[contig_gid]
        _family     = self.TD.taxon_family[contig_gid]
        _order      = self.TD.taxon_order[contig_gid]
        _class      = self.TD.taxon_class[contig_gid]
        _phylum     = self.TD.taxon_phylum[contig_gid]
        
        # get blastn description taxon tree
        blast_phylum, blast_class, blast_order, blast_family, blast_genus, blast_other = self.getBlastDescriptionTaxonTree(description_split) 
        
        # check contig taxa against blast taxa
        description = self.checkContigVsBlastTaxa(_phylum, blast_phylum, _class, blast_class, _order, blast_order, _family, blast_family, _genus, blast_genus, blast_other)
        
        ###########################################
        if description == "Inter Phyla":
            temp_output.write("%s\t%s\t%s\t%d\n" % (contig_gid, contig, original_description, self.contig_len[contig]))
        
        ###########################################
        
        return description
    
    def checkContigVsBlastTaxa(self, _phylum, blast_phylum, _class, blast_class, _order, blast_order, _family, blast_family, _genus, blast_genus, blast_other):
        # check at which point taxa matches up
        if _genus == blast_genus:
            return 'Genus'
        elif _family == blast_family:
            return 'Family'
        elif _order == blast_order:
            return 'Order'
        elif _class == blast_class:
            return 'Class'
        elif _phylum == blast_phylum:
            return 'Phylum'
        else:
            if len(blast_phylum) > 1:
                return 'Inter Phyla'
            else:
                return 'Other' 
            
            
    def getBlastDescriptionTaxonTree(self, description_split):
        
        _phylum = '' 
        _class  = ''
        _order  = ''
        _family = ''
        _genus  = ''
        _other  = ''
        
        for description in description_split:
            description = self.checkDescription(description)
            try:
                _phylum = self.GGT._phylum[description.lower()]
                _class  = self.GGT._class[description.lower()]
                _order  = self.GGT._order[description.lower()]
                _family = self.GGT._family[description.lower()]
                _genus  = description.lower()
            except KeyError:
                _other = " ".join(description_split)
                
        return _phylum, _class, _order, _family, _genus, _other
                    
    def checkDescription(self, description):
        if 'E.coli' in description or description == 'E.':
            description = 'Escherichia'
            return description
        elif description == 'C.coli':
            description = 'Campylobacter'
            return description
        elif 'Clostridium' in description:
            description = 'Clostridium'
            return description
        else:
            return description

    def addCorrectedData(self, description, count, dict):
        try:
            dict[description] += count
        except KeyError:
            dict[description] = count
    
    def addData(self, description, dict):
        try:
            dict[description] += 1 
        except KeyError:
            dict[description] = 1 

    def returnTGStatus(self, outfile):
        for TG in self.GD.group_data.keys():
            total = 0
            contam= 0 
            for pidsqid in self.GD.group_data[TG]:
                total += 1
                if self.checkPidsqid(pidsqid):
                    # pidsqid OK
                    pass
                else:
                    contam += 1 
            percentage = contam / float(total)
            self.transfer_groups[TG] = percentage
        
        # intialise outfile and print out TG status data
        outfile = open(outfile, 'w')
        
        for TG in self.transfer_groups.keys():
            outfile.write("%s\t%f\n" % (TG, self.transfer_groups[TG]))
                    
    def checkPidsqid(self, pidsqid):
        contig1 = self.HD.contig_name[pidsqid][0]
        contig2 = self.HD.contig_name[pidsqid][1]
        if contig1 in self.blast_data or contig2 in self.blast_data:
            return False
        else:
            return True
        
    def returnListOfPidsqids(self):
        pass
        
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
    C = Contigs(args.taxonomy_file,
                args.group_file,
                args.hitdata,
                args.gg_tax_file)
    C.wrapper(args.blastdir,
              args.cov_links_dir,
              args.type,
              args.tg_status_outfile,
              args.contaminated_contigs_file,
              args.length_min,
              args.length_max,
              args.coverage_min,
              args.coverage_max,
              args.show_plot,
              args.pie_graph_outfile,
              args.pie_graph_outfmt,
              args.pidsqid_outfile
              )
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('contaminated_contigs_file', help="")
    parser.add_argument('blastdir', help="")
    parser.add_argument('taxonomy_file', help="")
    parser.add_argument('group_file', help="")
    parser.add_argument('hitdata', help="")
    parser.add_argument('gg_tax_file', help="")
    parser.add_argument('cov_links_dir', help="")
    parser.add_argument('-t','--type', default = 'Pie_chart',help="Set the type of analysis to be completed. (default)Pie_chart: pie chart of specified contigs. tg_status: returns percentage of TGs contaminated. list: returns list of contaminated pidsqids.")
    parser.add_argument('-tgso','--tg_status_outfile', help="")
    parser.add_argument('-lmn','--length_min', type=int, default=0, help="")
    parser.add_argument('-lmx','--length_max', type=int, default=999999999, help="")
    parser.add_argument('-cmn','--coverage_min', type=float, default=0, help="")
    parser.add_argument('-cmx','--coverage_max', type=float, default=1000, help="")
    parser.add_argument('-sp','--show_plot', default=False, help="")
    parser.add_argument('-po','--pie_graph_outfile', default = 'piegraph', help="")
    parser.add_argument('-pidsqid_outfile','--pidsqid_outfile', default = 'pidsqid.csv', help="")
    parser.add_argument('-pofmt','--pie_graph_outfmt', default='png',help="")
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
