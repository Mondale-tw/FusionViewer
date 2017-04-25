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
            for row in tsvin:
                fusion = Fusion._make(row)
                fusions.append(fusion)
        return fusions

    def get_fusions(self):
        return self.fusions
        
    
