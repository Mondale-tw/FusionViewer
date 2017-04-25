"""reference.py: Reference Module

This module contains classes and functions for handling reference sequences,
which are usuallystored in FASTA and FAI (index) files.
"""

import pyfasta

class Reference(object):
    """The Reference class is the main interface to FASTA data.
    """
    
    def __init__(self, fasta):
        """Initialization of the Reference class 
        
        Load a fasta file.

        Args:
            fasta (str): path to the reference FASTA file

        Raises:
            Exception(): when loading failed
        """
        self.name = fasta
        self.fasta = pyfasta.Fasta(fasta)
        self.length = {}
        for chrom in self.fasta.keys():
            self.length[chrom] = len(self.fasta[chrom])
        
    def get_reference_name(self):
        """Get reference name(path)
        
        Args:
            none

        Return:
            str, the name (path) of the reference
        """
        return self.name

    def get_sequence(self, chr, start, end):
        """Get the raw sequence at "chr" from "start" to "end",
           as in one-based coordinate system.
        
        Args:
            chr (str): name of the chromosome
            start (int): starting position
            end (int): ending position, which will not be included in the  
                       returned sequence

        Return:
            str, the raw sequence
            
        Example:
            ref.get_sequence("chr1",1,9) should return bases 1~9 on chr1.
        """
        return self.fasta.sequence({'chr': chr, 'start': start, 'stop': end})

    def get_reference_length(self):
        """Get the reference length information

        Args:
            None

        Return:
            dict, a dictionary with chrmosome names as keys and lengths as
                  values
        """
        return self.length

