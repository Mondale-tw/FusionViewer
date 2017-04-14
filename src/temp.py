import pysam
import os
import sys
"print sys.argv"
class Bam():
    
    
    def __init__ (self, bam):
        """initialization of bam class
        
        Load bam file
        Args:
        bam (str) :is path to a bam file
        
        """
        self.name = bam
        self.bam  = pysam.AlignmentFile( bam , "rb")
        self.List = []
        for read in self.bam.fetch(until_eof=True ):
            self.seq = read.seq
            self.qname = read.qname
            self.pos = read.pos
            self.cigarstring = read.cigarstring
            self.pnext = read.pnext
            self.flag = read.flag
            #self.List.append(self.seq,self.qname,self.pos,self.cigarstring,self.pnext,self.flag)
            self.List.append(read.flag)
            
    def get_names (self):
        return self.List





will = Bam ("N.FusionReads.bam")

print will.get_names()
