from string import split,lower,upper
from random import random
import os

dir = os.path.dirname(__file__)
filename = os.path.join(dir, '../ref/best_expression_codons_Ecoli')

class CodonData( object ):
	def __init__(self,path=filename):
		self.__AA={}   #dictonary to list of freqency codon_pairs	
		self.__trans={}
		self.read_codon_table(path)

	def __getitem__(self, name):
		if name.upper() in self.__AA.keys():
			return  self.codon( name )
		elif name.lower() in self.__trans.keys():
			return self.translate( name )
		elif ( len( name.lower() ) == 3	and 'n' in name.lower() ):
			return '?'
		else: return len(name)*' '
		
	def codon(self, AA, best=True):
		AA=AA.upper()
		picked_codon=''

		if best: picked_codon = self.__AA[AA][-1][1]
		elif len(self.__AA[AA]) == 1: picked_codon = self.__AA[AA][0][1]
		else:	
			value = random()
			for elem in self.__AA[AA]:
				if value >= elem[0]: picked_codon = elem[1]
			if picked_codon =='': picked_codon =  self.__AA[AA][0][1] 
		return picked_codon 

	def all_codons( self, AA):
		""" Returns a list containing all of the codons for a specified AA 
		"""
		lst=[]
		for one in self.__AA[AA]:
			lst.append(one[1])
		return lst	
	
	def translate(self,codon):
		return self.__trans[codon.lower()]

	def read_codon_table(self,path):
		for line in open(path):
			
			if '#' in line: continue
			codon,AA,freq=split(line)
			codon=codon.lower()
			AA=AA.upper()
			freq=float(freq)
			
			self.__trans[codon]=AA
			
			if AA not in self.__AA.keys():
				self.__AA[AA]=[]
			self.__AA[AA].append( (freq,codon) )	

		for AA in self.__AA.keys(): 
			self.__AA[AA] = sorted(self.__AA[AA])
			sm=0
			for one in self.__AA[AA]:
				sm += one[0]
			# Renormalize the datatable just in case
			t=[]
			for val in self.__AA[AA]:
				t.append( (val[0]/sm, val[1]) )
			
			self.__AA[AA]=t

			sm=0
			for one in self.__AA[AA]:
				sm += one[0]
		


if __name__=='__main__':

	c=CodonData()
	c.codon('X',False)			 
	print c.all_codons('H')
