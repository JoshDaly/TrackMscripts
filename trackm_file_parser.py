#!/usr/bin/env python
###############################################################################
#
# __trackm_file_parser__.py - description!
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
# GNU General Public License fo r more details. #
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
import re as re

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class KmerFileParser(object):
    def __init__(self,l):
        self.readKmerFile(l)
        
    def readKmerFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqid        = tabs[0]
        self.gid1           = tabs[1]
        self.gid2           = tabs[2]
        self.kmer_score     = float(tabs[3])
        self.mean_dg        = tabs[4]
        self.dg1            = tabs[5]
        self.dg2            = tabs[6]
        
class KmerData(object):
    def __init__(self, kmer_file):
        self.kmer_scores            = {}
        self.dgg_values             = {}
        self.distance_from_neutral  = {}
        self.buildKmerData(kmer_file)
        
    def buildKmerData(self, kmer_file):
        with open(kmer_file) as fh:
            for l in fh:
                if l[0:3] != 'pid': 
                    KFP = KmerFileParser(l)
                    self.kmer_scores[KFP.pidsqid]           = KFP.kmer_score
                    self.dgg_values[KFP.pidsqid]            = [KFP.dg1, KFP.gid2]
                    self.distance_from_neutral[KFP.pidsqid] = abs(0.5-KFP.kmer_score)

###############################################################################                    
###############################################################################
###############################################################################
                
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
        self.taxon_genus        = {}
        self.taxon_family       = {}
        self.taxon_order        = {}
        self.taxon_class        = {}
        self.taxon_phylum       = {}
        self.taxon_string       = {}
        self.taxon_string_genus = {}
        self.taxon_organism     = {}
        self.gids               = {}
        self.buildTaxonData(taxon_data)
        
    def buildTaxonData(self, taxon_data):
        with open(taxon_data) as fh:
            for l in fh:
                TFP = TaxonomyFileParser(l)
                
                semi_colon_data = TFP.taxonomy.split(";")
                organism = TFP.organism.split()
                
                self.gids[TFP.gid] = 1
                
                _organism, _genus, _family, _order, _class, _phylum = self.getTaxonData(semi_colon_data, organism)
                
                self.addTaxonData(TFP.gid, _organism, _genus, _family, _order, _class, _phylum)         
                
                self.taxon_string[_organism] = TFP.taxonomy
                
                self.taxon_string_genus[TFP.gid] = "p__%s; c__%s; o__%s; f__%s; g__%s" % (_phylum, _class, _order, _family, _genus)  
    
    def getTaxonData(self, taxon_string, organism):
        _organism  = ''
        _genus     = ''
        _family    = ''
        _order     = ''
        _class     = ''
        _phylum    = ''
        if organism[0][0] == '"':
            _organism = organism[0][1:]
        else:
            _organism = organism[0]
        
        for taxon_rank in taxon_string:
            if 'p__' in taxon_rank:
                _phylum    = self.removeSemiColon(taxon_rank).lower()
            elif 'c__' in taxon_rank:
                _class     = self.removeSemiColon(taxon_rank).lower()
            elif 'o__' in taxon_rank:
                _order     = self.removeSemiColon(taxon_rank).lower()
            elif 'f__' in taxon_rank:
                _family    = self.removeSemiColon(taxon_rank).lower()
            elif 'g__' in taxon_rank:
                _genus     = self.removeSemiColon(taxon_rank).lower()
                
        return _organism, _genus, _family, _order, _class, _phylum
                
    def removeSemiColon(self, taxon_rank):
        if ';' in taxon_rank:
            return taxon_rank[4:-1]
        else:
            return taxon_rank[4:]
        
    def addTaxonData(self, gid, _organism, _genus, _family, _order, _class, _phylum):
        self.taxon_phylum[gid]  = _phylum
        self.taxon_class[gid]   = _class
        self.taxon_order[gid]   = _order
        self.taxon_family[gid]  = _family
        self.taxon_genus[gid]   = _genus
        self.taxon_organism[gid]= _organism

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
                for pid_sqid in GFP.pid_sqids:
                    try:
                        self.group_data[GFP.group_num].append(pid_sqid)
                    except KeyError:
                        self.group_data[GFP.group_num] = [pid_sqid]
                    
                    # add pid_sqid as members
                    self.group_membership[pid_sqid] = GFP.group_num

###############################################################################                    
###############################################################################
###############################################################################

class HitDataParser(object):
    def __init__(self,l):
        self.readHitData(l)
        
    def readHitData(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqid=                    tabs[0]
        self.transfer_group=             tabs[1]
        self.hid=                        tabs[2]
        self.pid=                        tabs[3]
        self.ani_comp=                   tabs[4]
        self.ident=                      tabs[5]
        self.gid_1=                      tabs[6]
        self.habitat_1=                  tabs[7]
        self.phylum_1=                   tabs[8]
        self.genus_1=                    tabs[9]
        self.status_1=                   tabs[10]
        self.sequencingMethod_1=         tabs[11]
        self.sequencingCentre_1=         tabs[12]
        self.horizontalTransferred_1=    tabs[13]
        self.genomeSize_1=               tabs[14]
        self.scaffoldCount_1=            tabs[15]
        self.len_1=                      tabs[16]
        self.strand_1=                   tabs[17]
        self.cid_1=                      tabs[18]
        self.contig_1=                   tabs[19]
        self.contigLength_1=             tabs[20]
        self.sqid_1=                     tabs[21]
        self.start_1=                    tabs[22]
        self.gid_2=                      tabs[23]
        self.habitat_2=                  tabs[24]
        self.phylum_2=                   tabs[25]
        self.genus_2=                    tabs[26]
        self.status_2=                   tabs[27]
        self.sequencingMethod_2=         tabs[28]
        self.sequencingCentre_2=         tabs[29]
        self.horizontalTransferred_2=    tabs[30]
        self.genomeSize_2=               tabs[31]
        self.scaffoldCount_2=            tabs[32]
        self.len_2=                      tabs[33]
        self.strand_2=                   tabs[34]
        self.cid_2=                      tabs[35]
        self.contig_2=                   tabs[36]
        self.contigLength_2=             tabs[37]
        self.sqid_2=                     tabs[38]
        self.start_2=                    tabs[39]
        
class HitData(object):
    def __init__(self, hitdata_file):
        self.hit_data                   = {}
        self.transfer_size              = {}
        self.transfer_size_gid          = {}
        self.contig_size                = {}
        self.contig_size_gid            = {}
        self.contig_name                = {}
        self.habitat                    = {}
        self.phylum                     = {}
        self.genus                      = {}
        self.pidsqid_to_gid             = {}
        self.status                     = {}
        self.gid_to_genus               = {}
        self.genus_to_phylum            = {}
        self.intra_or_inter             = {}
        self.genome_count_genus         = {}
        self.genome_status              = {}
        self.transfer_group_to_contigs  = {}
        self.tranfer_groups_per_genome  = {}
        self.gid_status                 = {}
        self.contig_to_gid              = {}
        self.contig_pidsqid_count_inter = {}
        self.contig_pidsqid_count_intra = {}
        self.contig_pidsqid_count_all   = {}
        self.buildHitData(hitdata_file)
        
        
    def buildHitData(self, hitdata_file):
        with open(hitdata_file) as fh:
            for l in fh:
                if l[0:3] != 'pid':
                    HDP = HitDataParser(l)
                    
                    # updated phylogeny
                    gid1    = HDP.gid_1
                    gid2    = HDP.gid_2
                    #pidsqid = "%s_%s" % (HDP.pid, HDP.sqid_1)
                    
                    # build data structures
                    self.phylumDistance(HDP.phylum_1, HDP.phylum_2, HDP.pidsqid)
                    self.transfer_size[HDP.pidsqid]                 = (int(HDP.len_1) + int(HDP.len_2))/float(2)
                    self.transfer_size_gid[HDP.pidsqid]             = [int(HDP.len_1), int(HDP.len_2)]
                    self.contig_size[HDP.pidsqid]                   = [int(HDP.contigLength_1), int(HDP.contigLength_2)]
                    self.contig_size_gid[HDP.pidsqid]               = [int(HDP.contigLength_1), int(HDP.contigLength_2)]
                    self.habitat[HDP.pidsqid]                       = [HDP.habitat_1, HDP.habitat_2]
                    self.phylum[HDP.pidsqid]                        = [HDP.phylum_1, HDP.phylum_2]
                    self.genus[HDP.pidsqid]                         = [HDP.genus_1, HDP.genus_2]
                    self.pidsqid_to_gid[HDP.pidsqid]                = [gid1, gid2]
                    self.contig_name[HDP.pidsqid]                   = [HDP.contig_1, HDP.contig_2]
                    self.status[HDP.pidsqid]                        = [HDP.status_1, HDP.status_2]
                    self.gid_status[HDP.gid_1]                      = HDP.status_1
                    self.gid_status[HDP.gid_2]                      = HDP.status_2
                    self.gid_to_genus[gid1]                         = HDP.genus_1
                    self.gid_to_genus[gid2]                         = HDP.genus_2
                    self.genus_to_phylum[HDP.genus_1]               = HDP.phylum_1
                    self.genus_to_phylum[HDP.genus_2]               = HDP.phylum_2
                    self.contig_to_gid[HDP.contig_1]                = HDP.gid_1
                    self.contig_to_gid[HDP.contig_2]                = HDP.gid_2
                    self.genomeCountPerGenus(HDP.genus_1, gid1)
                    self.genomeCountPerGenus(HDP.genus_2, gid2)
                    self.genome_status[HDP.pidsqid]                 = [HDP.status_1, HDP.status_2] 
                    self.contigsPerTransferGroup(HDP.transfer_group, HDP.contig_1)
                    self.contigsPerTransferGroup(HDP.transfer_group, HDP.contig_2)
                    self.addGenomeToTG(HDP.gid_1, HDP.transfer_group)
                    self.addGenomeToTG(HDP.gid_2, HDP.transfer_group)
                    self.addToContigPidsqidCount(HDP.contig_1, HDP.pidsqid)
                    self.addToContigPidsqidCount(HDP.contig_2, HDP.pidsqid)
                    # HDP.transfer_group,
    
                    self.hit_data[HDP.pidsqid]          = [HDP.transfer_group,
                                                           HDP.hid,
                                                           HDP.pid,
                                                           HDP.ani_comp,
                                                           HDP.ident,
                                                           HDP.gid_1,
                                                           HDP.habitat_1,
                                                           HDP.phylum_1,
                                                           HDP.genus_1,
                                                           HDP.status_1,
                                                           HDP.sequencingMethod_1,
                                                           HDP.sequencingCentre_1,
                                                           HDP.horizontalTransferred_1,
                                                           HDP.genomeSize_1,
                                                           HDP.scaffoldCount_1,
                                                           HDP.len_1,
                                                           HDP.strand_1,
                                                           HDP.cid_1,
                                                           HDP.contig_1,
                                                           HDP.contigLength_1,
                                                           HDP.sqid_1,
                                                           HDP.start_1,
                                                           HDP.gid_2,
                                                           HDP.habitat_2,
                                                           HDP.phylum_2,
                                                           HDP.genus_2,
                                                           HDP.status_2,
                                                           HDP.sequencingMethod_2,
                                                           HDP.sequencingCentre_2,
                                                           HDP.horizontalTransferred_2,
                                                           HDP.genomeSize_2,
                                                           HDP.scaffoldCount_2,
                                                           HDP.len_2,
                                                           HDP.strand_2,
                                                           HDP.cid_2,
                                                           HDP.contig_2,
                                                           HDP.contigLength_2,
                                                           HDP.sqid_2,
                                                           HDP.start_2
                                                           ]
                                                    
    def addToContigPidsqidCount(self, contig, pidsqid):
        
        if self.intra_or_inter[pidsqid] == 'inter':
            try:
                self.contig_pidsqid_count_inter[contig][pidsqid] = 1
            except KeyError:
                self.contig_pidsqid_count_inter[contig] = {pidsqid:1}
        elif self.intra_or_inter[pidsqid] == 'intra':
            try:
                self.contig_pidsqid_count_intra[contig][pidsqid] = 1
            except KeyError:
                self.contig_pidsqid_count_intra[contig] = {pidsqid:1}
        
        # add to total
        try:
            self.contig_pidsqid_count_all[contig][pidsqid] = 1
        except KeyError:
            self.contig_pidsqid_count_all[contig] = {pidsqid:1}
    
    def addGenomeToTG(self, gid, TG):
        try:
            self.tranfer_groups_per_genome[gid] += TG
        except KeyError:
            self.tranfer_groups_per_genome[gid] = [TG]
    
    def contigsPerTransferGroup(self, TG, contig):
        try:
            self.transfer_group_to_contigs[contig][TG] = 1
        except KeyError:
            self.transfer_group_to_contigs[contig] = {TG:1}
    
    def phylumDistance(self, phylum1, phylum2, pidsqid):
        if phylum1 == phylum2:
            self.intra_or_inter[pidsqid] = 'intra'
        else:
            self.intra_or_inter[pidsqid] = 'inter'
            
    def genomeCountPerGenus(self,genus,gid):
        try:
            self.genome_count_genus[genus][gid] = 1
        except KeyError:
            self.genome_count_genus[genus] = {gid:1}

###############################################################################                    
###############################################################################
###############################################################################

class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnoFile(l)
        
    def readAnnoFile(self,l):
        tabs = l.rstrip().split("\t")
        self.pidsqidgene        = tabs[0]
        self.pidsqid            = "_".join(tabs[0].split("_")[0:2])
        self.identity           = tabs[1]
        self.evalue             = tabs[2]
        self.length             = tabs[3]
        self.cog_id             = tabs[5]
        self.cog_description    = tabs[6]
        self.img_annotation     = tabs[8]
        
class AnnotationData(object):
    def __init__(self, anno_file):
        self.COGs               = {}
        self.cog_to_pidsqid     = {}
        self.cog_description    = {}
        self.IMG_anno_data      = {}
        self.buildAnnoData(anno_file)
    
    def buildAnnoData(self,anno_file):
        # read in anno file
        with open(anno_file) as fh:
            for l in fh:
                AFP = AnnotationFileParser(l)
                
                # check COG is present 
                if self.checkCOGStatus(AFP.cog_id):
                    try:
                        self.COGs[AFP.pidsqid].append(AFP.cog_id)
                    except KeyError:
                        self.COGs[AFP.pidsqid] = [AFP.cog_id]
                        
                    # build cog data
                    self.addCOG(AFP.cog_id, AFP.pidsqid)
                    self.cog_description[AFP.cog_id]    = AFP.cog_description
                        
                # build IMG annotation data
                try:
                    self.IMG_anno_data[AFP.pidsqid].append(AFP.img_annotation)
                except KeyError: 
                    self.IMG_anno_data[AFP.pidsqid] = [AFP.img_annotation]
                        
    def addCOG(self, cog, pidsqid):
        try:
            self.cog_to_pidsqid[cog][pidsqid] = 1
        except KeyError:
            self.cog_to_pidsqid[cog] = {pidsqid:1}
    
    def checkCOGStatus(self, COGid):
        if len(COGid) > 0:
            return True
        else:
            return False

###############################################################################                    
###############################################################################
###############################################################################

class PathsFileParser(object):
    def __init__(self,l):
        self.readPathsFile(l)
    
    def readPathsFile(self,l):
        tabs = l.rstrip().split("\t")
        self.gid             = tabs[0]
        self.path_to_file    = tabs[1]
        self.img             = tabs[1].split("/")[-2]
        
class PathsFileData(object):
    def __init__(self, paths_file):
        self.gid_to_file        = {}
        self.img_to_gid         = {}
        self.gid_to_img         = {}
        self.buildPathsData(paths_file)
    
    def buildPathsData(self, paths_file):
        with open(paths_file) as fh:
            for l in fh:
                if l[0] != "#":
                    PFP = PathsFileParser(l)
                    PFP.readPathsFile(l)
                    self.gid_to_file[PFP.gid] = PFP.path_to_file
                    self.img_to_gid[PFP.img]  = PFP.gid
                    self.gid_to_img[PFP.gid]  = PFP.img
                
        

###############################################################################                    
###############################################################################
###############################################################################

class BlastFileParser(object):
    def __init__(self, blast_file):
        self.blast_data     = {}
        self.blast_data_tg  = {}
        self.topHits        = {}
        self.topHits2       = {}
        self.isvector       = {}
        self.contig_len     = {}
        self.parseBlastFile(blast_file)
    
    def parseBlastFile(self, blast_file):
        with open(blast_file) as fh:
            for l in fh:
                if l[0:6] == "Query=":
                    query_id = l.rstrip().split("=")[-1].split()[0]
                    # initialise dict for query
                    self.isvector[query_id] = 0
                    
                    for l in fh:
                        if l[0:7] == "Length=":
                            self.contig_len[query_id] = int(l.rstrip().split("=")[-1])

                        #if re.search('[a-z]+"|"[a-z][1-9]+',l):
                        if re.search('^[a-z]+\|',l):
                            try:
                                whitespace      = l.rstrip().split()
                                straightline    = whitespace[0].rstrip().split('|')
                                subject_id      = straightline[1]
                                score = float(whitespace[-2])
                                evalue = float(whitespace[-1])
                                description = " ".join(whitespace[1:-2])
                                
                                # collect data
                                self.addToBLASTData(query_id, subject_id, score, evalue, description)
                                self.addToBLASTDataTG(query_id, subject_id, score, evalue, description)
                                
                                # get topHit
                                self.getTopHit(query_id, description, evalue)
                                self.getTopHit2(query_id, description, evalue)
                                self.getVector(query_id, description)
                            except (IndexError, ValueError):
                                pass
                            
                        elif l[0:5] == '*****':
                    
                            self.topHits[query_id]  = 'NH' # NH= No Hits
                            self.topHits2[query_id] = ['NH',0]
                            
                            # no hits! 
                            break
                        elif l[0:6] == "Query_":
                            break
        
        # set longest contig query as default taxonomy
        try:
            self.setContigQueryDefault()
        except KeyError:
            pass
        
    def setContigQueryDefault(self):
        longest_contig_len  = 0
        longest_contig      = ''
         
        for contig in self.contig_len.keys():
            if self.contig_len[contig] > longest_contig_len:
                longest_contig = contig
                longest_contig_len = self.contig_len[contig]
                
        # set default taxonomy to longest contig
        default_taxon = self.topHits[longest_contig].split(' ')[0]
        
        for contig in self.topHits.keys():
            if default_taxon in self.topHits[contig]:
                self.topHits[contig] = 'MuyBueno'
            elif self.topHits[contig] == 'NH':
                pass
            else:
                self.topHits[contig] = 'xeno'
            
    def getTopHit2(self, query_id, description, evalue):
        try:
            lenny = self.topHits2[query_id]
            pass
        except KeyError:
            self.topHits2[query_id] = [description, evalue]
    
    def getVector(self, query_id, description):
        vector_search = ['vector','homo sapien','human', 'mus musculus', 'medicago', 'xenopus', 'mouse']
        for vector in vector_search:
            if vector in description.lower():
                try:
                    self.isvector[query_id] +=1 
                except KeyError:
                    self.isvector[query_id] =1
                break 
    
    def getTopHit(self, query_id, description, evalue):
        try:
            lenny = self.topHits[query_id]
            pass
        except KeyError:
            if evalue <= 0.00001:
                self.topHits[query_id] = description
            else: 
                self.topHits[query_id] = 'NH'
    
    def addToBLASTData(self, query_id, subject_id, score, evalue, description):
        try:
            self.blast_data[query_id][subject_id].append([score, evalue, description])
        except KeyError:
            try:
                self.blast_data[query_id][subject_id] = [score, evalue, description]
            except KeyError:
                self.blast_data[query_id] = {subject_id:[score, evalue, description]}
    
    def addToBLASTDataTG(self, query_id, subject_id, score, evalue, description):
        if query_id > subject_id:
            try:
                self.blast_data_tg[query_id][subject_id] = [evalue, score, description]
            except KeyError:
                self.blast_data_tg[query_id] = {subject_id:[evalue, score, description]}
        else:
            try:
                self.blast_data_tg[subject_id][query_id] = [evalue, score, description]
            except KeyError:
                self.blast_data_tg[subject_id] = {query_id:[evalue, score, description]}

###############################################################################                    
###############################################################################
###############################################################################

class KEGGFileParser(object):
    def __init__(self,l):
        self.readKEGGFile(l)
        
    def readKEGGFile(self, l):
        tabs                    = l.rstrip().split("\t")
        self.geneid             = tabs[0]
        self.pidsqid_gene       = tabs[1]
        self.pidsqid            = "_".join(tabs[1].split("_")[:-1])
        self.kegg_id            = tabs[2]
        self.definition         = tabs[3]
        self.organism           = tabs[4]
        self.identity           = tabs[5]
        self.col4               = tabs[6]
        self.col5               = tabs[7]
        self.col6               = tabs[8]
        self.col7               = tabs[9]
        self.col8               = tabs[10]
        self.col9               = tabs[11]
        self.col10              = tabs[12]
        self.evalue             = tabs[13]
        self.col12              = tabs[14]
    
class KEGGData(object):
    def __init__(self, kegg_file):
        self.kegg_data          = {}
        self.all_kegg_data      = {}
        self.kegg_to_pidsqid    = {}
        self.kegg_definition    = {}
        self.buildKEGGData(kegg_file)
        
    def buildKEGGData(self, kegg_file):
        with open(kegg_file) as fh:
            for l in fh:
                KFP = KEGGFileParser(l)
                try:
                    self.kegg_data[KFP.pidsqid].append(KFP.kegg_id)
                except KeyError:
                    self.kegg_data[KFP.pidsqid] = [KFP.kegg_id]
                try:
                    self.kegg_to_pidsqid[KFP.kegg_id].append(KFP.pidsqid)
                except KeyError:
                    self.kegg_to_pidsqid[KFP.kegg_id] = [KFP.pidsqid]
                # add description to dict
                self.kegg_definition[KFP.kegg_id] = KFP.definition  
                
                self.all_kegg_data[KFP.pidsqid_gene] = [KFP.pidsqid_gene,
                                                        KFP.pidsqid,
                                                        KFP.kegg_id,
                                                        KFP.definition,
                                                        KFP.organism,
                                                        KFP.identity,
                                                        KFP.col4,
                                                        KFP.col5,
                                                        KFP.col6,
                                                        KFP.col7,
                                                        KFP.col8,
                                                        KFP.col9,
                                                        KFP.col10,
                                                        KFP.evalue,
                                                        KFP.col12]

###############################################################################                    
###############################################################################
###############################################################################

class KEGGMappingData(object):
    def __init__(self, kegg_mapping_file):
        self.kegg_mapping   = {}
        self.pathway_lookup = {}
        self.wrapper(kegg_mapping_file)
        
    def wrapper(self, kegg_mapping_file):
        with open(kegg_mapping_file) as fh:
            line_number = 0 
            for l in fh: 
                commas = l.rstrip().split(',')
                if line_number == 0:
                    # Kegg Pathway No.'s
                    pass
                    
                elif line_number == 1: 
                    # Kegg Pathway Descriptions
                    self.captureKEGGPathways(commas)
                    
                else:
                    # data
                    self.captureKEGGData(commas)
                line_number += 1 
                    
    def captureKEGGData(self, delimited_line):
        column_number = 0
        for col in delimited_line:
            if column_number == 0:
                pass
            else:
                if int(col) == 1:
                    pathway = self.pathway_lookup[column_number]
                    kegg_id = delimited_line[0]
                    try:
                        self.kegg_mapping[kegg_id][pathway] = 1 
                    except KeyError:
                        self.kegg_mapping[kegg_id] = {pathway:1}
            column_number += 1
                    
    def captureKEGGPathways(self, delimited_line):
        column_number = 0
        for col in delimited_line:
            if column_number >= 1:
                self.pathway_lookup[column_number] = col
            column_number += 1 

###############################################################################                    
###############################################################################
###############################################################################

class KEGGSummaryTableParser(object):
    def __init__(self, l):
        self.readKEGGSummaryFile(l)
        
    def readKEGGSummaryFile(self, l):
        tabs            = l.rstrip().split("\t")
        self.kegg            = tabs[0]
        self.num_genus       = tabs[1]
        self.num_phyla       = tabs[2]
        self.num_habitats    = tabs[3]
        self.transfers_avg   = tabs[4]
        self.transfer_std    = tabs[5]
        self.contig_avg      = tabs[6]
        self.contig_std      = tabs[7]
        self.transfer_groups = tabs[8]

class KEGGSummaryTableData(object):
    def __init__(self, kegg_summary_table_file):
        self.kegg_list = {}
        self.wrapper(kegg_summary_table_file)
        
    def wrapper(self, kegg_summary_table_file):
        line_number = 0
        with open(kegg_summary_table_file) as fh:
            for l in fh:
                if line_number == 0:
                    pass
                else:
                    KSTP = KEGGSummaryTableParser(l)
                    self.kegg_list[KSTP.kegg] = 1 
                line_number += 1 

###############################################################################                    
###############################################################################
###############################################################################

class KEGGLookupTable(object):
    def __init__(self, kegg_lookup_table): 
        self.KEGG_lookup = {}        
        self.collateKEGGLookupTable(kegg_lookup_table)
    
    def collateKEGGLookupTable(self, kegg_lookup_table):
        with open(kegg_lookup_table) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                KO      = tabs[0] 
                kegg_id = tabs[1]
                
                self.KEGG_lookup[kegg_id] = KO 
            
###############################################################################                    
###############################################################################
###############################################################################
    
class BamMLinksFileParser(object):
    def __init__(self, l):
        self.readBamMLinksFile(l)
        
    def readBamMLinksFile(self,l):
        tabs = l.rstrip().split("\t")
        #cid_1  cid_2   len_1   pos_1   rev_1   len_2   pos_2   rev_2   file
        self.cid_1 = tabs[0] 
        self.cid_2 = tabs[1]
        self.len_1 = tabs[2]
        self.pos_1 = tabs[3]
        self.rev_1 = tabs[4]
        self.len_2 = tabs[5]
        self.pos_2 = tabs[6]
        self.rev_2 = tabs[7]
        self.file  = tabs[8]
        
class BamMCoverageFileParser(object):
    def __init__(self, l):
        self.readBamMCoverageFile(l)
        
    def readBamMCoverageFile(self,l):
        tabs = l.rstrip().split("\t")
        #contig Length  649989915.SRR067990.bam 649989915.SRR072965_1.bam
        self.cid                = tabs[0]
        self.cid_len            = tabs[1]
        self.num_coverage_files = len(tabs)
        self.covs               = []
        for i in range(2, self.num_coverage_files):
            self.covs.append(float(tabs[i])) 

class BamMData(object):
    def __init__(self, bamm_links_file, bamm_cov_file):
        self.bamm_links         = {}
        self.bamm_links_data    = {}
        self.bamm_links_total   = {}
        self.bamm_cov_data      = {}
        self.contig_len         = {}
        self.wrapper(bamm_links_file,
                     bamm_cov_file)
        
    def wrapper(self, bamm_links_file, bamm_cov_file):
        with open(bamm_links_file) as fh:
            for l in fh:
                if l[0] != '#':
                    BLFP = BamMLinksFileParser(l)
                    
                    self.buildLinksDataAll(BLFP.cid_1,
                                           BLFP.cid_2,
                                           BLFP.len_1,
                                           BLFP.len_2,
                                           BLFP.pos_1,
                                           BLFP.pos_2
                                           )
                    
                    # build links data
                    self.buildLinksData(BLFP.cid_1, BLFP.cid_2)
                    self.buildLinksData(BLFP.cid_2, BLFP.cid_1)
        
        with open(bamm_cov_file) as fh:
            for l in fh:
                if l[0] != '#':
                    BCFP = BamMCoverageFileParser(l)
                    
                    self.contig_len[BCFP.cid] = BCFP.cid_len
                    
                    # build coverage data
                    for cov in BCFP.covs:
                        self.buildCovData(BCFP.cid, cov)
                    
                    
    def buildLinksDataAll(self, cid1, cid2, len1, len2, pos1, pos2):
        try:
            self.bamm_links[cid1][cid2].append([len1, len2, pos1, pos2])
        except KeyError:
            try:
                self.bamm_links[cid1][cid2] = [[len1, len2, pos1, pos2]]
            except KeyError:
                self.bamm_links[cid1] = {cid2:[[len1, len2, pos1, pos2]]}
        try:
            self.bamm_links[cid1][cid1].append(len1, len2, pos1, pos2)
        except KeyError:
            try:
                self.bamm_links[cid2][cid1] = [[len2, len1, pos2, pos1]]
            except KeyError:
                self.bamm_links[cid2] = {cid1:[[len2, len1, pos2, pos1]]}
            
    
    def buildCovData(self, cid, cov):
        try:
            self.bamm_cov_data[cid].append(cov)
        except KeyError:
            self.bamm_cov_data[cid] = [cov]
                    
    def buildLinksData(self, cid1, cid2):
        try:
            self.bamm_links_data[cid1][cid2] += 1 
        except KeyError:
            try:
                self.bamm_links_data[cid1][cid2] = 1
            except KeyError:
                self.bamm_links_data[cid1] = {cid2:1}
        try:
            self.bamm_links_total[cid1] += 1 
        except KeyError:
            self.bamm_links_total[cid1] = 1 
            
        try:
            self.bamm_links_total[cid2] += 1
        except KeyError:
            self.bamm_links_total[cid2] = 1
    
###############################################################################                    
###############################################################################
###############################################################################

class CoverageLinksFileParser(object):
    def __init__(self, l):
        self.readCoverageLinksFile(l)
        
    def readCoverageLinksFile(self,l):
        # gid     phylum  genus   contig_id       contig_length   contig_coverage avg_genome_coverage     contig_links    avg_genome_links        transfer_groups
        tabs = l.rstrip().split("\t")
        self.gid                    = tabs[0] 
        self.phylum                 = tabs[1]
        self.genus                  = tabs[2]
        self.contig_id              = tabs[3]
        self.contig_len             = tabs[4]
        self.contig_coverage        = tabs[5]
        self.avg_genome_coverage    = tabs[6]
        self.contig_links           = tabs[7]
        self.avg_genome_links       = tabs[8]
        self.transfer_groups        = []
        for i in tabs[9].split(";"):
            self.transfer_groups.append(i)
        self.call                   = tabs[10]
        
class CoverageLinksData(object):
    def __init__(self, coverage_links_file):
        self.TGs_contigs            = {}
        self.contigs_TGs            = {}
        self.contigs                = []
        self.TGs                    = []
        self.gids                   = {}
        self.contig_to_genus        = {}
        self.TG_status              = {}
        self.wrapper(coverage_links_file)
        
    def wrapper(self, coverage_links_file):
        with open(coverage_links_file) as fh:
            for l in fh:
                if l[0] != 'g':
                    CLFP = CoverageLinksFileParser(l)
                    self.contigs.append(CLFP.contig_id)
                    self.contig_to_genus[CLFP.contig_id] = CLFP.genus
                    for TG in CLFP.transfer_groups:
                        if len(TG) > 0:
                            self.addTG(TG)
                            self.addContigToTG(CLFP.contig_id, TG)
                            self.TG_status[TG] = CLFP.call

    def addContigToTG(self, contig, TG):
        try:
            self.TGs_contigs[TG][contig] = 1
        except KeyError:
            self.TGs_contigs[TG] = {contig:1}
        
        try:
            self.contigs_TGs[contig][TG] = 1
        except KeyError:
            self.contigs_TGs[contig] = {TG:1}
        
        
    def addTG(self, TG):
        if TG not in self.TGs:
            self.TGs.append(TG)

###############################################################################                    
###############################################################################
###############################################################################

class CovLinksFileParser(object):
    def __init__(self, l):
        self.readCovLinksFile(l)
    
    def readCovLinksFile(self,l):
        tabs = l.rstrip().split("\t")
        self.col_count              = len(tabs)
        self.contig_id              = tabs[0]
        self.contig_len             = tabs[1]
        self.relative_cov           = tabs[2]
        if len(tabs) == 6:
            self.contig_links       = tabs[3]
            self.links_div_chromo   = tabs[4]
            self.links_div_total    = tabs[5]
        
class CovLinksData(object):
    def __init__(self, cov_links_file):
        self.cov_link_data = {}
        self.wrapper(cov_links_file)
        
    def wrapper(self, cov_links_file):
        with open(cov_links_file) as fh:
            gid = cov_links_file.split("/")[-1].split("_")[0]
            for l in fh:
                if l[0:6] != 'contig':
                    CLFP = CovLinksFileParser(l)
                    if CLFP.col_count == 6:
                        self.addAllData(gid, CLFP.contig_id, CLFP.contig_len, CLFP.relative_cov, CLFP.contig_links, CLFP.links_div_chromo, CLFP.links_div_total)
                    else:
                        self.addCovData(gid, CLFP.contig_id, CLFP.contig_len, CLFP.relative_cov)
                        
    def addAllData(self, gid, contig, contig_len, rel_cov, contig_links, links_div_chromo, links_div_total):
        try:
            self.cov_link_data[gid][contig] = [contig_len, rel_cov, contig_links, links_div_chromo, links_div_total]
        except KeyError:
            self.cov_link_data[gid] = {contig:[contig_len, rel_cov, contig_links, links_div_chromo, links_div_total]}
        
    def addCovData(self, gid, contig, contig_len, rel_cov):
        try:
            self.cov_link_data[gid][contig] = [contig_len, rel_cov]
        except KeyError:
            self.cov_link_data[gid] = {contig: [contig_len, rel_cov]}
        
###############################################################################                    
###############################################################################
###############################################################################

class CoverageFileParser(object):
    def __init__(self,l):
        self.readCoverageFile(l)
        
    def readCoverageFile(self,l):
        tabs=l.rstrip().split("\t")
        self.gid    = tabs[0]
        self.covs   = []
        for i in tabs[1:]:
            self.covs.append(i)

class CoverageData(object):
    def __init__(self, coverage_file):
        self.coverage_data  = {}
        self.wrapper(coverage_file)
    
    def wrapper(self, coverage_file):
        with open(coverage_file) as fh:
            for l in fh:
                CFP = CoverageFileParser(l)
                for i in CFP.covs:
                    colon = i.rstrip().split(":")
                    contig_id           = colon[0] 
                    rel_cov             = float(colon[1])
                    links_divide_total  = float(colon[2])
                    contig_len          = int(colon[3])
                    try:
                        self.coverage_data[CFP.gid][contig_id] = [rel_cov, links_divide_total, contig_len]
                    except KeyError:
                        self.coverage_data[CFP.gid] = {contig_id:[rel_cov, links_divide_total, contig_len]}
                    
###############################################################################                    
###############################################################################
###############################################################################

class TransferGroupSummaryFileParser(object):
    def __init__(self, l):
        self.readTransferGroupSummaryFile(l)
        
    def readTransferGroupSummaryFile(self,l):
        tabs = l.rstrip().split("\t")
        self.transfer_group         = tabs[0]
        self.transfer_size_mean     = tabs[1]    
        self.transfer_size_std      = tabs[2]
        self.contig_size_mean       = tabs[3]    
        self.contig_size_std        = tabs[4]    
        self.habitat_no             = tabs[5]    
        self.phylum_no              = tabs[6]   
        self.genus_no               = tabs[7]   
        self.pidsqids_no            = tabs[8]   
        self.most_connected_genus   = tabs[9]    
        self.kmer_score_mean        = tabs[10]    
        self.kmer_score_std         = tabs[11] 
        try:  
            self.kegg_ids               = tabs[12]
        except IndexError:
            pass
        
class TransferGroupSummaryData(object):
    def __init__(self, transfer_group_summary_file):
        self.transfer_group_data = {}
        self.wrapper(transfer_group_summary_file)
        
    def wrapper(self, transfer_group_summary_file):
        with open(transfer_group_summary_file) as fh:
            for l in fh:
                if l[0] != 't':
                    TGSFP = TransferGroupSummaryFileParser(l)
                    self.transfer_group_data[TGSFP.transfer_group] = [TGSFP.transfer_size_mean,
                                                                      TGSFP.transfer_size_std,
                                                                      TGSFP.contig_size_mean,
                                                                      TGSFP.contig_size_std,
                                                                      TGSFP.habitat_no,
                                                                      TGSFP.phylum_no,
                                                                      TGSFP.genus_no,
                                                                      TGSFP.pidsqids_no,
                                                                      TGSFP.most_connected_genus,
                                                                      TGSFP.kmer_score_mean,
                                                                      TGSFP.kmer_score_std
                                                                      ]
                    try:
                        self.transfer_group_data[TGSFP.transfer_group].append(TGSFP.kegg_ids)
                    except AttributeError:
                        pass
                                                            
###############################################################################                    
###############################################################################
###############################################################################
###############################################################################

class PairsFileParser(object):
    def __init__(self, l): 
        self.readPairsFile(l)
        
    def readPairsFile(self, l):
        commas = l.rstrip().split(',')
        self.gid_1      = commas[0]
        self.gid_2      = commas[1]
        self.ani_comp   = float(commas[3])
            
class PairsFileData(object):
    def __init__(self, pairs_data_file):
        self.pairs_data = {}
        self.wrapper(pairs_data_file)
    
    def wrapper(self, pairs_data_file):
        with open(pairs_data_file) as fh:
            for l in fh:
                if l[0] != 'g':
                    PFP = PairsFileParser(l)
                    self.addPairsData(PFP.gid_1, PFP.gid_2, PFP.ani_comp)
                    self.addPairsData(PFP.gid_2, PFP.gid_1, PFP.ani_comp)
    
    def addPairsData(self, gid1, gid2, ani_comp):
        try:
            self.pairs_data[gid1][gid2] = ani_comp
        except KeyError:
            self.pairs_data[gid1] = {gid2:ani_comp}
            
###############################################################################                    
###############################################################################
###############################################################################
###############################################################################    
    
class NucMerParser:
    """Wrapper class for parsing nucmer output"""
    # constants to make the code more readable
    _START_1  = 0
    _END_1    = 1
    _START_2  = 2
    _END_2    = 3
    _LEN_1    = 4
    _LEN_2    = 5
    _IDENTITY = 6
    _ID_1     = 7
    _ID_2     = 8

    def __init__(self):
        self.prepped = False

    def readNuc(self, fp):
        """Read through a nucmer coords file

        this is a generator function
        """
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fp: # search for the first record
                        if l[0] == '=': # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fp:
                fields = l.split('|')
                yield ([int(i) for i in fields[0].split()] +
                       [int(i) for i in fields[1].split()] +
                       [int(i) for i in fields[2].split()] +
                       [float(i) for i in fields[3].split()] +
                       fields[4].split())
            break # done! 
        
###############################################################################                    
###############################################################################
###############################################################################
###############################################################################

class ContaminatedPidsqidsFileParser(object):
    def __init__(self, l):
        self.readContaminatedPidsqidsFile(l)
        
    def readContaminatedPidsqidsFile(self, l):
        tabs = l.rstrip().split("\t")
        self.category   = tabs[0]
        self.pidsqids   = tabs[1].split(",")

class ContaminatedPidsqids(object):
    def __init__(self, contam_pidsqids_file):
        self.contam_pidsqids        = {}
        self.all_contam_pidsqids    = {}
        self.wrapper(contam_pidsqids_file)
        
    def wrapper(self, contam_pidsqids_file): 
        with open(contam_pidsqids_file) as fh:
            for l in fh:
                CPFP = ContaminatedPidsqidsFileParser(l)
                
                for pidsqid in CPFP.pidsqids:
                    self.all_contam_pidsqids[pidsqid] = 1
                
                # conservative threshold 
                if CPFP.category == 'Inter Phyla' or CPFP.category == 'Vector' or CPFP.category == 'PhiX': 
                    for pidsqid in CPFP.pidsqids:
                        self.contam_pidsqids[pidsqid] = 1
        print 'Total number of pidsqids affected %d' % (len(self.all_contam_pidsqids))
        
###############################################################################                    
###############################################################################
###############################################################################
###############################################################################    

class ContaminatedContigsFileParser(object):
    def __init__(self, l):
        self.readContaminatedContigsFile(l)
        
    def readContaminatedContigsFile(self, l):
        tabs = l.rstrip().split("\t")
        self.contigid   = tabs[0]
        self.category   = tabs[1]
        
class ContaminatedContigData(object):
    def __init__(self, contam_contigs_file):
        self.contam_contigs = {}
        self.wrapper(contam_contigs_file)
    
    def wrapper(self, contam_contigs_file):
        with open(contam_contigs_file) as fh:
            for l in fh:
                CCFP = ContaminatedContigsFileParser(l)
                try:
                    self.contam_contigs[CCFP.category][CCFP.contigid] = 1
                except KeyError:
                    self.contam_contigs[CCFP.category] = {CCFP.contigid:1}
                    
###############################################################################                    
###############################################################################
###############################################################################
###############################################################################
        