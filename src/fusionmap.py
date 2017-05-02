"""
This module contains classes and functions for handling fusionmap reports.
Reference and fusionmap fusion reads BAM files are optional. 
"""

import bam
import reference
import csv
import collections

Fusion = collections.namedtuple('Fusion',
             [  "ID",
                "UniqueCuttingPositionCount",
                "SeedCount",
                "RescuedCount",
                "Strand",
                "Chromosome1",
                "Position1",
                "Chromosome2",
                "Position2",
                "KnownGene1",
                "KnownGeneStrand1",
                "KnownGene2",
                "KnownGeneStrand2",
                "FusionJunctionSequence",
                "FusionGene",
                "Filter" ])

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])
    return lcs_set


class Report(object):
    """The FusionMap class is the interface for fusionmap report file.
    """

    def __init__(self, report, fusion_reads=None , reference=None ):
        """Initialization of the FusionMap class

        Load a fusionmap report text file. The refernece file and 
        the associated bam file are optional but recommended.

        Args:
            report (str): path to the fusionmap report file
            fusion_reads (str): path to the fusion reads bam file
            reference (str): path to the refernce FASTA file

        Raises:
            Exception(): when loading failed
        """

        self.fusions = self.__parse(report)
        self.final_fusions = self._merge_fusions()
        if fusion_reads:
            self.bam = bam.Bam(fusion_reads)
        else:
            self.bam = None
        
        if reference:
            self.reference = reference.Reference(reference)
        else:
            self.reference = None

    def __parse(self, report_path):
        """Parse the fusionmap report"""
        fusions = []
        with open(report_path,'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            counter = 0
            for row in tsvin:
                if counter == 0:
                    counter = counter + 1
                    continue
                fusion = Fusion._make(row)
                fusions.append(fusion)
        return fusions

    def get_fusions(self):
        return self.final_fusions
        
    def _merge_fusions(self):
        """Merge the fusions returned by fusionmap"""
        raw_fusions = self.fusions
        bucket = {}
        for fusion in raw_fusions:
            if (fusion.Chromosome1, fusion.Chromosome2) not in bucket:
                bucket[(fusion.Chromosome1, fusion.Chromosome2)] = [fusion]
            else:
                bucket[(fusion.Chromosome1, fusion.Chromosome2)].append(fusion)
        final_fusions = []
        for key, fusions in bucket.iteritems():
            final_fusions = final_fusions + self._merge_fusion_same_chroms(fusions)
        return final_fusions
    
    def _merge_fusion_same_chroms(self, fusions, threshold =0.8):
        """Merging fusions which are known from same chromosomes"""
        dimension = len(fusions)
        merging_set = []
        for x in range(dimension):
            for y in range(x,dimension):
                if ( x != y and
                     fusions[x].KnownGene1 == fusions[y].KnownGene1 and
                     fusions[x].KnownGene2 == fusions[y].KnownGene2):
                    longest_common_sequence = lcs(fusions[x].FusionJunctionSequence.lower(),
                                                  fusions[y].FusionJunctionSequence.lower())
                    ratio = float(len(longest_common_sequence)*2)/(len(fusions[x].FusionJunctionSequence) + len(fusions[y].FusionJunctionSequence))
                    if ratio > threshold:
                        merging_set.append((x,y))
        
        duplicates = set()
        for x,y in merging_set:
            duplicates.add(y)
        
        unique_fusions = [fusions[x] for x in range(dimension) if x not in duplicates]
        return unique_fusions
            

