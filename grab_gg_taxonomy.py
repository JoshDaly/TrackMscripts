#!/usr/bin/env python
###############################################################################
#
# __grab_gg_taxonomy__.py - description!
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
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GreenGenesTaxonomy(object):
    def __init__(self, gg_tax_file):
        self._phylum        = {}
        self._class         = {}
        self._order         = {}
        self._family        = {}
        self.wrapper(gg_tax_file)
    
    def wrapper(self, gg_tax_file):
        self.grabGGTaxData(gg_tax_file)
    
    def grabGGTaxData(self, gg_tax_file):
        with open(gg_tax_file) as fh:
            for l in fh:
                whitespace  = l.rstrip().split()
                gg_id       = whitespace[0]
                gg_tax      = whitespace[1:]
                try:
                    __genus, __family, __order, __class, __phylum = self.getTaxonRanks(gg_tax)
                    self.addTaxData(__genus, __family, __order, __class, __phylum)
                except IndexError:
                    pass
    
    def getTaxonRanks(self, gg_tax):
        __genus     = ''
        __family    = ''
        __order     = ''
        __class     = ''
        __phylum    = ''
        
        for taxon_rank in gg_tax:
            if 'p__' in taxon_rank:
                __phylum    = self.removeSemiColon(taxon_rank).lower()
            elif 'c__' in taxon_rank:
                __class     = self.removeSemiColon(taxon_rank).lower()
            elif 'o__' in taxon_rank:
                __order     = self.removeSemiColon(taxon_rank).lower()
            elif 'f__' in taxon_rank:
                __family    = self.removeSemiColon(taxon_rank).lower()
            elif 'g__' in taxon_rank:
                __genus     = self.removeSemiColon(taxon_rank).lower()
                
        return __genus, __family, __order, __class, __phylum

    def removeSemiColon(self, taxon_rank):
        if ';' in taxon_rank:
            return taxon_rank[3:-1]
        else:
            return taxon_rank[3:]

    def addTaxData(self, __genus, __family, __order, __class, __phylum):
        if len(__genus) > 1:
            self._phylum[__genus]   = __phylum
            self._class[__genus]    = __class
            self._order[__genus]    = __order
            self._family[__genus]   = __family
                
###############################################################################
###############################################################################
###############################################################################
###############################################################################
