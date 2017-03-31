"""reference.py: Reference Module

This module contains classes and functions for handling reference sequences,
which are usuallystored in FASTA and FAI (index) files.
"""

import pyfasta


class Reference():
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

    def get_sequence(self, chr, start, end):
        """Get the raw sequence at "chr" from "start" to "end"-1
        
        Args:
            chr (str): name of the chromosome
            start (int): starting position
            end (int): ending position, which will not be included in the  
                       returned sequence

        Return:
            str, the raw sequence
            
        Example:
            ref.get_sequence("chr1",1,9) should return the bases 1~8 on chr1.
        """

    def get_info(self):
        """Get the reference information

        Args:
            None

        Return:
            dict, a dictionary containing all chrmosome names and lengths
        """

