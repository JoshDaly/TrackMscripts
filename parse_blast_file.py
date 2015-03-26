#!/usr/bin/env python
###############################################################################
#
# __parse_blast_file__.py - Parse blast output file! Type default.
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
import re
import networkx as nx
import matplotlib.pyplot as plt
#import graph_tool.all as gt
import code
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

import numpy as np
np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class AnnotationFileParser(object):
    def __init__(self,l):
        self.readAnnotationData(l)

    def readAnnotationData(self,l):
        tabs            = l.rstrip().split("\t")
        self.pid        = tabs[0].split("_")[0]
        self.sqid       = tabs[0].split("_")[1]
        self.gene       = tabs[0].split("_")[2]
        self.annotation = tabs[1]

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
        self.bodysite_1=                 tabs[5]
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
        self.bodysite_2=                 tabs[22]
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

class Data(object):
    def __init__(self, name):
        self.__name  = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

class BLASTData(object):
    def __init__(self):
        self.blastData      = {}
        self.geneList       = {} # unique list of genes
        self.pidList        = {} # unique list of interPhyla PIDs
        self.phylaTuples    = {}
        self.annotation     = {}
        self.vertexGene     = {}

        # pid_sqid BLAST fmt
        self.pidSqid        = {}

    def getBLASTDataPidSqid(self,file, blast_score_cutoff, evalue_cutoff):
        """parse BLAST output file"""
        with open(file, 'r') as fh:
            for l in fh:
                if l[0:6] == "Query=":
                    query = l.rstrip().split(" ")[-1].split("_")
                    pid1  =  query[0]
                    sqid1 =  query[1]

                    # build node list
                    self.geneList[pid1+"_"+sqid1] = 1

                    for l in fh:
                        if l[0:4] == 'lcl|':
                            whitespace = l.rstrip().split()
                            straightline = whitespace[0].rstrip().split('|')
                            underscore = straightline[1].split("_")
                            pid2 = underscore[0]
                            sqid2 = underscore[1]
                            score = float(whitespace[1])
                            evalue = float(whitespace[2])
                            
                            # set arbitrary cutoff for blast score
                            if score >= blast_score_cutoff and evalue <= evalue_cutoff:
                                if pid1 > pid2:
                                    try:
                                        self.blastData[pid1+"_"+sqid1][pid2+"_"+sqid2] = 1
                                    except KeyError:
                                        self.blastData[pid1+"_"+sqid1] = {pid2+"_"+sqid2:1}
                                else:
                                    try:
                                        self.blastData[pid2+"_"+sqid2][pid1+"_"+sqid1] = 1
                                    except KeyError:
                                        self.blastData[pid2+"_"+sqid2] = {pid1+"_"+sqid1:1}

                                # build node list
                                self.geneList[pid2+"_"+sqid2] = 1
                        elif "Effective" in l:
                            #last_pos = fh.tell()
                            #fh.seek(last_pos)
                            break

    def getSubSetBLASTData(self,file,cutoff):
        """parse BLAST output file"""
        count = 0
        with open(file, 'r') as fh:
            for l in fh:
                if "Query=" in l:

                    if count > cutoff:
                        print "broken yo %d %d" % (count,cutoff)
                        break


                    query = l.rstrip().split(" ")[-1].split("_")
                    pid1 = query[0]
                    sqid1 = query[1]
                    gene1 = query[2]

                    #build node list
                    self.geneList[sqid1+"_"+gene1] = 1

                    for l in fh:
                        bam = re.match('  [0-9]',l)
                        if bam:
                            #print l.rstrip()
                            whitespace = l.rstrip().split()
                            underscore = whitespace[0].split("_")
                            pid2 = underscore[0]
                            sqid2 = underscore[1]
                            gene2 = underscore[2]
                            score = whitespace[1]
                            evalue = whitespace[2]
                            # build interaction data (edges)
                            try:
                                self.blastData[sqid1+"_"+gene1][sqid2+"_"+gene2] = 1
                            except KeyError:
                                self.blastData[sqid1+"_"+gene1] = {sqid2+"_"+gene2:1}
                            try:
                                self.blastData[sqid2+"_"+gene2][sqid1+"_"+gene1] = 1
                            except KeyError:
                                self.blastData[sqid2+"_"+gene2] = {sqid1+"_"+gene1:1}
                            """
                            if pid1 > pid2:
                                try:
                                    self.blastData[sqid1+"_"+gene1][sqid2+"_"+gene2] = 1
                                except KeyError:
                                    self.blastData[sqid1+"_"+gene1] = {sqid2+"_"+gene2:1}
                            else:
                                try:
                                    self.blastData[sqid2+"_"+gene2][sqid1+"_"+gene1] = 1
                                except KeyError:
                                    self.blastData[sqid2+"_"+gene2] = {sqid1+"_"+gene1:1}
                            """
                            # build node list
                            self.geneList[sqid2+"_"+gene2] =1
                        elif "Effective" in l:
                            #last_pos = fh.tell()
                            #fh.seek(last_pos)
                            break
                    count += 1

    def getBLASTData(self,file):
        """parse BLAST output file"""
        with open(file, 'r') as fh:
            #last_pos = fh.tell()
            for l in fh:
                    #print l.rstrip()
                if l[0:6] == "Query=": # in l:
                    query = l.rstrip().split(" ")[-1].split("_")
                    pid1 = query[0]
                    sqid1 = query[1]
                    gene1 = query[2]

                    #build node list
                    self.geneList[sqid1+"_"+gene1] = pid1

                    # add PID to list
                    self.pidList[pid1] = 1

                    for l in fh:
                        bam = re.match('  [0-9]',l)
                        if bam:

                            whitespace = l.rstrip().split()
                            underscore = whitespace[0].split("_")
                            pid2 = underscore[0]
                            sqid2 = underscore[1]
                            gene2 = underscore[2]
                            score = whitespace[1]
                            evalue = whitespace[2]
                            # build interaction data (edges)
                            """
                            try:
                                self.blastData[sqid1+"_"+gene1][sqid2+"_"+gene2] = 1
                            except KeyError:
                                self.blastData[sqid1+"_"+gene1] = {sqid2+"_"+gene2:1}
                            try:
                                self.blastData[sqid2+"_"+gene2][sqid1+"_"+gene1] = 1
                            except KeyError:
                                self.blastData[sqid2+"_"+gene2] = {sqid1+"_"+gene1:1}
                            """


                            if pid1 > pid2:
                                try:
                                    self.blastData[sqid1+"_"+gene1][sqid2+"_"+gene2] = 1
                                except KeyError:
                                    self.blastData[sqid1+"_"+gene1] = {sqid2+"_"+gene2:1}
                            else:
                                try:
                                    self.blastData[sqid2+"_"+gene2][sqid1+"_"+gene1] = 1
                                except KeyError:
                                    self.blastData[sqid2+"_"+gene2] = {sqid1+"_"+gene1:1}

                            # build node list
                            self.geneList[sqid2+"_"+gene2] = pid2

                            # add to PID list
                            self.pidList[pid2] = 1

                        elif "Effective" in l:
                            #last_pos = fh.tell()
                            #fh.seek(last_pos)
                            break

    def getPhylogeneticData(self,file):
        with open(file,'r') as fh:
            for l in fh:
                HDP = HitDataParser(l)
                if l[0:3] != 'hid':
                    HDP.readHitData(l)
                    if HDP.pid in self.pidList: # interPhyla transfer
                        """pid -> (phylum1,genus1,phylum2,genus2)"""
                        self.phylaTuples[HDP.pid] = (HDP.phylum_1,HDP.genus_1,HDP.phylum_2,HDP.genus_2)

    def getAnnotationData(self,file):
        with open(file,'r') as fh:
            for l in fh:
                AFP = AnnotationFileParser(l)
                if AFP.pid in self.pidList:
                    self.annotation[AFP.pid] = AFP.annotation

    def makeNetworkPlot(self):
        """Create network graph"""
        print "Building network plot:"
        print "Creating %d nodes" % (len(self.geneList))
        G=nx.Graph()

        # add genes (nodes)
        G.add_nodes_from(self.geneList.keys())


        # add edges
        for gene1 in self.blastData.keys():
            for gene2 in self.blastData[gene1]:
                if gene1 in self.geneList and gene2 in self.geneList:
                    G.add_edge(gene1,
                               gene2)

        # Build network plot
        pos= nx.spring_layout(G,iterations=500)
        #fig = plt.figure(figsize=(21,10),dpi=300)

        #plt.subplot(1,1,1,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,pos,linewidths=0, alpha=0.7, node_size=10,with_labels=True)
        nx.draw_networkx_edges(G, pos, edge_color = "#373737",width=0.3)
        nx.draw_networkx_labels(G, pos)
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')

        plt.show()

    def getConnectedNodes(self,group_number):
        """Connect dem nodes!"""

        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_2_nodes= {}

        for region in self.blastData.keys():
            if region not in regions_2_nodes:
                D = Data(region)
                regions_2_nodes[region] = D
                nodes |= {D}
            for linked_region in self.blastData[region].keys():
                if linked_region not in regions_2_nodes:
                    D = Data(linked_region)
                    regions_2_nodes[linked_region] = D
                    nodes |= {D}

        for region1 in self.blastData.keys():
            for region2 in self.blastData[region1].keys():
                regions_2_nodes[region1].add_link(regions_2_nodes[region2])

        # Find all the connected components
        number = 1
        # drop to interpretor!
        #code.interact(local=locals())

        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            
            if group_number >0:
                if names_count == group_number:
                    for name in names:
                        print name
                    break
                else:
                    number += 1
            else:
                names = ",".join(names)
                print "Group %i:\t%d\t%s" % (number,names_count, names)
                number += 1
                
    def makeNetworkPlotGraphTool(self,outFile):
        """make network plot using Graph Tool!"""
        # initialise graph
        G = gt.Graph(directed=False)
        geneVertices   = {}
        count = 0

        # convert genes to vertices, numbered # add genes (nodes)
        for gene in self.geneList.keys():
            try:
                if len(self.blastData[gene]) > 30:
                    geneVertices[gene] = G.add_vertex()
                    self.vertexGene[count] = gene
                    count += 1
            except KeyError:
                pass

        # add nodes and edges
        for gene1 in self.blastData.keys():
            if len(self.blastData[gene1]) > 50:
                for gene2 in self.blastData[gene1]:
                    if len(self.blastData[gene2]) > 30:
                        if gene1 != gene2:
                            if gene1 in self.geneList and gene2 in self.geneList:
                                G.add_edge(geneVertices[gene1],
                                           geneVertices[gene2])
        # set pos
        pos = gt.sfdp_layout(G)
        #pos = gt.arf_layout(G, d=0.1)
        #pos = gt.get_hierarchy_control_points(G)
        #pos = gt.fruchterman_reingold_layout(G)

        # draw network plot
        gt.graph_draw(G,
                      pos=pos,
                      vertex_text=G.vertex_index,
                      output_size=(20000,20000),
                      vertex_font_size=100,
                      output=outFile,
                      nodesfirst=True,
                      fmt='png',
                      vertex_size = 1,
                      edge_marker_size = 1
                      )
                            #

    def makeInteractionScatterMap(self,outfile,show_plot=True):
        """make an interaction map using a scatter plot"""
        xs = []
        ys = []
        val = []
        workingIds = self.geneList.keys() # master list of genome tree ids

        #set the number of tick labels
        #for i in range(len(workingIds)):
        #    ticks.append(i)
        #ax.set_xticks(ticks)
        #ax.set_yticks(ticks)

        #Change label names according to taxonomy
        #labels = [item.get_text() for item in ax.get_xticklabels()]
        #for i,v in enumerate(labels):
        #    labels[i] = ordered_tax_string_lowest[i]
        #ax.set_xticklabels(labels,size=5)
        #ax.set_yticklabels(labels,size=5)

        #
        for x in range(len(workingIds)):
            try:
                fred = self.blastData[workingIds[x]]
                for y in range(x+1, len(workingIds)):
                    try:
                        # build upper triangle data
                        val.append(fred[workingIds[y]])
                        xs.append(y)
                        ys.append(x)
                    except KeyError:
                        pass
            except KeyError:
                pass
        #-----
        # build scatter plot
        plt.subplot(1,1,1,autoscale_on=True)
        plt.tight_layout()
        #Plot labels
        plt.xlabel("transferred sequences")
        plt.ylabel("transferred sequences")
        plt.title("Interaction map of interphyla transfers")
        #ax.grid(b=True,which='both')
        plt.scatter(xs,
                    ys,
                    marker='s',
                    alpha = 1,
                    s = 0.1,
                    linewidths = 0
                    )



        #-----

        #rotate x labels 90deg
        #plt.xticks(rotation=90)

        #adjust margin size
        #plt.subplots_adjust(bottom=0.2)


        #plt.xticks(np.arange(min(xs), max(xs)+1,1.0))
        #plt.yticks(np.arange(min(xs), max(xs)+1,1.0))
        # set the plot axes
        #plt.axis([plot_border*-1,
        #          len(workingIds)+plot_border,
        #          plot_border*-1,
        #          len(workingIds)+plot_border])

        if show_plot:
            plt.show()
        else:
            plt.savefig(IPS.outFile,dpi=IPS.dpi,format=IPS.imageFormat)

    def buildUpperMatrix(self):
        workingIds = self.geneList.keys()
        matrix = []



        # build empty matrix
        for x in range(len(workingIds)):
            matrix.append([])
            for y in range(len(workingIds)):
                matrix[x].append(0)
        print "matrix dimensions:"
        print "x = %d" % len(matrix)


        # add data to matrix
        for x in range(len(workingIds)):
            for y in range(x+1,len(workingIds)):
                try:
                    fred = self.blastData[workingIds[x]][workingIds[y]]
                    matrix[x][y] = 1
                except KeyError:
                    pass
        np_matrix = np.array(matrix)


        # single linkage clustering using scipy
        #matrix_collapsed = ssd.squareform(matrix)
        single_linkage = sch.linkage(np_matrix, method='single', metric='euclidean')
        print single_linkage

    def returnGeneID(self,vertexID):

        print self.vertexGene[vertexID]

    def connected_components(self,nodes):
        # List of connected components found. The order is random.
        result = []

        # Make a copy of the set, so we can modify it.
        nodes = set(nodes)

        nn = float(len(nodes))
        remains = nn

        # Iterate while we still have nodes to process.
        while nodes:

            # Get a random node and remove it from the global set.
            n = nodes.pop()

            # This set will contain the next group of nodes connected to each other.
            group = {n}

            # Build a queue with this node in it.
            queue = [n]

            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:

                # Consume the next item from the queue.
                n = queue.pop(0)

                # Fetch the neighbors.
                neighbors = n.links

                # Remove the neighbors we already visited.
                neighbors.difference_update(group)

                # Remove the remaining nodes from the global set.
                remains -= float(len(neighbors))
                #print "%f completed" % (1 - remains/nn)
                nodes.difference_update(neighbors)

                # Add them to the group of connected nodes.
                group.update(neighbors)

                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(group)

        # Return the list of groups.
        return result

    def printHeader(self):
        print "\t".join(['sqid_gene',
                         'pid',
                         'count',
                         'phylum1',
                         'genus1',
                         'phylum2',
                         'genus2'
                         ])

    def printInteractionData(self):
        #print self.blastData
        for gene1 in self.blastData.keys():
            count = 0
            pid =  self.geneList[gene1]
            for gene2 in self.blastData[gene1]:
                count += 1
            try:
                print "\t".join([gene1,
                                 pid,
                                 str(count),
                                 self.phylaTuples[pid][0],
                                 self.phylaTuples[pid][1],
                                 self.phylaTuples[pid][2],
                                 self.phylaTuples[pid][3],
                                 self.annotation[pid]])
            except KeyError:
                pass
                #print "This key %s gone buggered up!" % (pid)



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
    #objects
    BD = BLASTData()

    if args.printData:
        BD.getBLASTData(args.blast_file)
        BD.getAnnotationData(args.anno_file) # need file
        BD.getPhylogeneticData(args.hitData_file) # need file
        BD.printHeader()
        BD.printInteractionData()
    elif args.graphTool:
        if args.pid_sqid_fmt:
            if args.cutoff > 0:
                BD.getSubSetBLASTData(args.blast_file, args.cutoff)
                BD.makeNetworkPlotGraphTool(args.outFile)
            else:
                BD.getBLASTDataPidSqid(args.blast_file)
                BD.makeInteractionScatterMap(args.outFile)
                #BD.makeNetworkPlotGraphTool(args.outFile)
        else:
            if args.cutoff > 0:
                BD.getSubSetBLASTData(args.blast_file, args.cutoff)
                BD.makeNetworkPlotGraphTool(args.outFile)
            else:
                BD.getBLASTData(args.blast_file)
                BD.makeNetworkPlotGraphTool(args.outFile)
    elif args.vertexID:
        if args.cutoff > 0:
            BD.getSubSetBLASTData(args.blast_file, args.cutoff)
            BD.makeNetworkPlotGraphTool(args.outFile,args.vertexID)
            BD.returnGeneID(args.vertexID)
        else:
            BD.getBLASTData(args.blast_file)
            BD.makeNetworkPlotGraphTool(args.outFile,args.vertexID)
            BD.returnGeneID(args.vertexID)

    elif args.group_nodes:
        """Group nodes based on single linkage clustering"""
        BD.getBLASTDataPidSqid(args.blast_file, 
                               args.blast_score_cutoff,
                               args.evalue_cutoff)
        BD.getConnectedNodes(args.group)

    else:
        if args.cutoff > 0:
            BD.getSubSetBLASTData(args.blast_file, args.cutoff)
            BD.printInteractionData()
        else:
            # Parse Blast data, and build data store
            BD.getBLASTData(args.blast_file)

        # build network plot of blast results
        BD.makeNetworkPlot()


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file', help="blast output file")
    parser.add_argument('-c','--cutoff', type=int, default=0, help="Set cutoff for subset analysis")
    parser.add_argument('-p','--printData', default=False, help="Print out interaction data")
    parser.add_argument('-a','--anno_file', default=False, help="Provide table with annotation data")
    parser.add_argument('-hit','--hitData_file', default=False, help="Provide table in hitData format")
    parser.add_argument('-gt','--graphTool', default=False, help="Use graph tool instead of networkx")
    parser.add_argument('-o','--outFile', default='network_test',help="Output file name")
    parser.add_argument('-vid','--vertexID', type=int,help="Search for and return Gene ID using vertex ID")
    parser.add_argument('-pidsqid','--pid_sqid_fmt', default=False,help="pid_sqid format BLAST input file")
    parser.add_argument('-gn','--group_nodes', default=False,help="Group nodes by single linkage clustering")
    parser.add_argument('-g','--group', type=int,default=0,help="Print out group-specific files")
    parser.add_argument('-bsc','--blast_score_cutoff', type=int,default=160,help="Set Blast score cutoff")
    parser.add_argument('-evc','--evalue_cutoff', type=int,default=1,help="Set Evalue score cutoff")
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
