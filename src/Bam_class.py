#i make a class to import a SAM file and Parse it"
import pysam

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
			
			self.qname = read.qname
			self.flag = read.flag
			self.pos = read.pos
			self.cigarstring = read.cigarstring
			self.pnext = read.pnext
			self.seq = read.seq
			y = [self.qname,self.flag,self.pos,self.cigarstring,self.pnext,self.seq]
			# y is a list to accumulate all the parsing value
			self.List.append(y)
			

	def get_parsebam (self):
		return self.List