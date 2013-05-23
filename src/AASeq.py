from re import finditer
from string import upper, split
from sys import path
from CodonData import CodonData

path.insert(0, '../')
import ref.aa_data as ref
import ref.dna_data as dna_ref

def get_ORFs(seq):
	""" Trivial split wrapper. Takes a protein single letter sequence and returns 
	the string split according to stop codons (1-Letter code 'X') 
	"""
	return seq.split('X') 

class InvalidAASeq( BaseException ):
	pass

class AASeq( object ):
	def __init__( self, seq ):				
		self.seq=''
		self.__clean_seq(seq)

	def __len__(self): return len(self.seq)
	def __str__(self): return self.seq	
			
	def __clean_seq( self, in_seq ):
		self.seq,s='',''
		for letter in in_seq.upper():
			if letter not in ref.AA:
				if letter == ' ': continue 
				print 'Warning the sequence member '+letter+' is not a recognized Amino Acid! Skipping'
				raise InvalidAASeq
			s+=letter.upper()

		orfs=get_ORFs(in_seq)
		for one in orfs:
			if len(one) > len(self.seq): self.seq=one
	
	def MW( self ):
		"""Calculates the molecular weight of self.seq in Daltons
		"""
		weight=17.0
		for letter in self.seq:
			if letter == 'X': break
			weight+=ref.AAMW[letter]
		return round(weight,1)

	def E280( self ):
		"""Returns the calculated molar absorptivity of protein fragment self.seq
		calculated at 280nm in water.
		"""
		nTyr=len( [m.start() for m in finditer('Y', self.seq)])
		nTrp=len( [m.start() for m in finditer('W', self.seq)])
		nCys=len( [m.start() for m in finditer('C', self.seq)])
		return round(( (nTyr*1490)+(nTrp*5500)+(nCys*125) ),1)	

	def suggest_dna(self, best=False, CodonObj=None):
		"""Translates the self.seq protein sequence to DNA
		If best == False : will pick a random codon with frequency proportional to frequency table
		If best == True  : will always pick the codon with the highest frequency in frequency table
		"""
		dna=''
		if not CodonObj:  CodonObj = CodonData()	
		for ltr in self.seq:
			dna+=CodonObj.codon( ltr, best )
		
		return dna
			

if __name__ == '__main__':
	import DNASeq as DNAlib
	seq='MEFFGESWKKHLSGEFGKPYFIKLMGFVAEERKHYTVYPPPHQVFTWTQMCDIKDVKVVILGQDPYHGPNQAHGLCFSVQRPVPPPPSLENIYKELSTDIEDFVHPGHGDLSGWAKQGVLLLNAVLTVRAHQANSHKERGWEQFTDAVVSWLNQNSNGLVFLLWGSYAQKKGSAIDRKRHHVLQTAHPSPLSVYRGFFGCRHFSKTNELLQKSGKKPIDWKEL'
	g=AASeq(seq)
	print g.MW()
	print g.E280()
	s= g.suggest_dna()
	g=DNAlib.DNASeq(s)
	g.auto_phase()
	
	ans = DNAlib.translate(g.trips).replace(' ','')
	print ans, ans==seq
	
