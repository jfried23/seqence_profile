from re import finditer
from sys import exit, path
import AASeq as AAlib
from CodonData import CodonData
from string import join, upper, lower, center, rjust, ljust

path.insert(0, '../')
import ref.dna_data as ref

class WrongAA(Exception):
	def __init__(self, AA):
		self.AA = AA

###########UTILITTY FUNCTIONS ##########################################################################

def split_codons( seq, phase=0):
	""" Splits the input sequence into a list of sequence triplets. eq [ 'ATG','CGA','GCN',...]
            If phase != 0 (and/or the sequence is not evenly divisible into units of 3)
            the first ( and/or last ) triplet will be padded with empty chars ' ' to make it
	    a triplet.
	"""
	seq.replace(' ','')
	phase = phase % 3
	i,ii=0,len(seq)
	trips=[]
	if phase != 0:
		i = phase
		trips.append( rjust( seq[0:phase], 3 )) 		
	while ( i < ii ):
		trips.append( seq[i:i+3] )
		i+=3
	trips[-1]= ljust(trips[-1],3)
	return trips

def translate( trips, CodonObj = None ):
	"""
  	Returns translated protein sequence from a list of sequence triplets
        Seq is out put with a space on either side of the 1-letter amino-acid code   
	"""
	if not CodonObj:  CodonObj = CodonData()

	s=''
	for trip in trips:
		s+=center( CodonObj[trip], 3 )
	return s

def count_bases( seq ):
	""" returns a dictionary object mapping each unique char in a sequence to the
	    the number of occurances.
	"""
	count={'a':0,'t':0,'g':0,'c':0,'n':0}
	for ltr in seq:
		if ltr not in count.keys(): count[ltr] =1
 		else: count[ltr]+=1
	if 'n' not in count.keys(): count['n']=0
	return count
 
def GC_content(seq, num=None):
	""" Returns the fractional GC content of self.seq. Must be in range [0.0,1.0]  """
	if not num: num = count_bases(seq)

	seq=seq.replace(' ','').lower()
	return (1.*num['g']+num['c'])/len(seq) #Damn python. Hackey but keeps it real!


def MW( seq, num = None ):
	""" 
	Returns the molecular weight of the DNA sequence 
	If unknown nucleotides 'n' are present in the sequence uses an
        average MW of ~309 g/Mol for that position.
        This function assumes the 5' is phosphorylated.   
	"""
	
	if not num: num = count_bases(seq)
	#if num['n']!=0: print "Caution MW calculation will use an average base MW for the 'N' nucleotides in this sequence!"
        return round((num['a']*313.21)+(num['t']*304.2)+(num['c']*289.18)+(num['g']*329.21)+(num['n']*308.95),1)


def E260( seq, num=None ):
	""" Returns the calculated molar extinction coefficient of DNA sequence
	    for single and double stranded DNA in units of 1/(Mol*L) using nearest neighbor 
            extinction coefficients of bases. If unknown nucleotides 'n' are present in the 
            sequences calculates extinction coefficients using 33 A/ng and 55 A/ng averages
            that are valid for long DNA sequences.	
	"""
	
	if not num: num = count_bases(seq)

	if num['n'] != 0:
		#print "Caution calculating E260 for DNA with unknown \'n\' nucleotides! Will calculate from averages."
		mw = MW(seq,num)
		sigmaSS1= mw/(33.e-3)
		sigmaDS = (2*mw)/(50.e-3)
		return round(sigmaSS1,1), round(sigmaDS,1)
					
	sigmaSS1,sigmaSS2, sigmaDS = 0,0,0
	rc=reverse_complement(seq)
	i,ii=0,len(seq)
	while i < ii:
		sigmaSS1+= ref.NA_E260[ seq[i:i+2] ]
		sigmaSS2+= ref.NA_E260[ rc[i:i+2] ]
		if (i!=0): 
			sigmaSS1-= ref.NA_E260[ seq[i] ]
			sigmaSS2-= ref.NA_E260[ rc[i] ]
		i+=1

	frac_gc = GC_content(seq, num)
	hypo_chrome = 0.059*frac_gc + 0.287*(1.-frac_gc) 
	sigmaDS = (1-hypo_chrome)*(sigmaSS1 + sigmaSS2)

	return round(sigmaSS1,1), round(sigmaDS,1)

def reverse_complement(seq):
	"""
	Returns the complement DNA sequence with same polarty.  5'->3' input maps to 5'->3' output
	nucleotides are assigned complement partners of 'n'
	"""
	s=''
	for letter in seq: s+=ref.DNA_comp[letter]
	return s[::-1]

def complement(seq):
	"""
	Returns the complement DNA sequence with same polarty.  5'->3' input maps to 5'->3' output
	Unknown 'n' nucleotides are assigned complement partners of 'n'
	"""
	s=''
	for letter in seq: s+=ref.DNA_comp[letter]
	return s

def clean_dna_seq( in_seq ):
	seq=''
	for letter in in_seq.lower():
		if letter not in ref.DNA_comp.keys():
			if letter == ' ': continue
			print letter 
			raise InvalidDNASeq()
		else: seq+=letter.lower()		
	return seq

class InvalidDNASeq( BaseException ):
	pass
#########################################################################################################

#############DNASeq Class definition#####################################################################
class DNASeq( object ):
	def __init__( self, seq, phase = 0 ):
		self.number_offset, self.phase = divmod(phase,3)
		self.seq=''
		self.__clean_seq(seq)

		self.trips = split_codons(self.seq, self.phase)
		if len(self.trips[0].replace(' ','')) !=3: self.number_offset+=1

		self.orf=(0,len(self)-1)

	def __len__(self): return len(self.seq)

	def __str__(self): return self.info_str()		
	
	def __clean_seq(self, in_seq):
		seq=''
		for letter in in_seq.lower():
			if letter not in ref.DNA_comp.keys():
				if letter == ' ': continue 
				raise InvalidDNASeq()
			else: seq+=letter.lower()		
		self.seq = seq

	def fetch_position( self, number, expected_AA=None ):
		""" Returns the start position of the first base in the triplet encoding
		    sequence position number. As optional sanity check the char expected_AA will
		    ensure the expected AA is in the designated position
		"""
		CodonObj = CodonData()
		if  expected_AA != None:
			wttrip = self.trips[ (self.number_offset + number)]

			if ( CodonObj[wttrip] != expected_AA ): raise WrongAA(  CodonObj[wttrip] )

		number = ( (self.number_offset + number) * 3 )
		if ' ' in self.trips[0]: number -= self.trips[0].count(' ')
		 
		return number
			

	def digest( self ):
		"""
		Returns a dictionary mapping a restriction enzyme cut site to a restriction enzyme name
		ie {target_site_position:"EnzymeName"}
		"""
		digest={}
		for enzyme in ref.restriction.keys():
			cut_site,offset=ref.restriction[enzyme]
			cuts = [ x+offset for x in [m.start() for m in finditer(cut_site, self.seq)] ]
			for c in cuts: digest[c]=enzyme	
		return digest

	def insert_restrict( self, enzyme_name, pos, phase ):
		"""
		Inserts a restriction site at the specified AA positions with the specified phase.
		Phase can be positive or negative and is not generally limited to the range [0,3)
		"""
		if enzyme_name not in ref.restriction.keys():
			print "Accepted enzymes ", ref.restriction.keys()
			exit()
		pos = self.fetch_position(pos-1) + phase
		new_seq = ref.restriction[ enzyme_name ][0]
		self.seq = self.seq[0:pos] + new_seq + self.seq[pos+len(new_seq):]
		self.set_phase(self.phase) 
		return		

	def info_str(self):
		info_string=''
		es,ed= E260(self.seq)
		info_string += "****************** DNA SEQEUENCE *********************\n"
		info_string += rjust( str(len(self)), 5 ) +" base pairs.\n"
		info_string += "E260 = " + rjust(str(round(es/1.e3, 2)), 9) +" 1/( mM * L) single stranded DNA\n"
		info_string += "E260 = " + rjust(str(round(ed/1.e3,2)), 9) +" 1/( mM * L) double stranded DNA\n"
		info_string +=  "Single stranded molecular weight: " + str(round(MW(self.seq)/1.e3,1))+"kDa.\n"
		
		return info_string

	def translate( self ):
		return translate( self.trips )

	def print_digest(self):
		splitstr=''
		digest=self.digest()
		ii=0
		for site in sorted(digest.keys()):
			interject = ' <--'+digest[site]+'--> '
			splitstr += self.seq[ii:site]+interject
			ii=site
		splitstr += self.seq[ii:]
		return splitstr	
		
	def fancy_print( self, lw=78, extras=False ):
		line=''
		if extras: line+=self.info_str()

		prefix = (self.trips[0].count(' '))*' '   ######Fixed to give the proper specified for 1st line
		seq=prefix+self.seq
		num_lines,remain= divmod( len(seq), lw )
		aaSeq= translate(self.trips)
		
		digest=self.digest()
		cut_line,key_line=' '*len(self),' '*len(self)
		for key in digest.keys():
			cut_line = cut_line[0:key]+'|'+cut_line[key+1:]
			key_line = key_line[0:key-3] + center( digest[key], 6 ) + key_line[key+3:]	

		cut_line = prefix+cut_line 
		key_line = prefix+key_line

		number=''
		n=1
		for i, codon in enumerate(self.trips):
			if ( i < self.number_offset ): 
				number+= ' '*3
				continue
			elif n==1:
				number+= center( str(n), 3)
			elif ( n % 5 == 0 ): 
				number+= center( str(n), 3)
			else: number += ' '*3
			n+=1

		if remain != 0: num_lines+=1
		for i in range(0,num_lines):
			#if i == 0: 
			if i == num_lines: lw=remain
			line += number[i*lw:(i+1)*lw]+'\n'
			line +=  aaSeq[i*lw:(i+1)*lw]+'\n'
			line +=  seq[i*lw:(i+1)*lw]+'\n'
			line +=  cut_line[i*lw:(i+1)*lw]+'\n'
			line +=  key_line[i*lw:(i+1)*lw]+'\n\n'
		return line
			
	def set_phase( self, phase):
		self.number_offset, self.phase = divmod(phase,3)
		self.trips=split_codons(self.seq,self.phase)
		if len(self.trips[0].replace(' ','')) !=3: self.number_offset+=1
			
		
	def auto_phase( self, sub_seq=None ):
		"""
		Finds a reading frame containing the input protein sequence (sub_seq). If no sub_seq is provided or
		the sub_seq is not found in any open reading frame, it uses the phase that generates the longest ORF.
		Returns a tuple (ORF_seq, found_sub) where ORF_seq is translated protein sequence of the longest open
		reading frame and found_sub is a boolean indicating if a specified sub_seq is contained
		within the ORF_seq. If argument sub_seq is default (None) found_sub will default to False
		"""
		best_phi, best_ORF = 0, ''
		found_sub = False
		for phi in range(0,3):
			self.number_offset, self.phase = divmod(phi,3)
			self.trips=split_codons(self.seq, phi)
			if len(self.trips[0].replace(' ','')) !=3: self.number_offset+=1
			this_seq=translate(self.trips).replace(' ','')
			
			orfs = AAlib.get_ORFs(this_seq)			
			for orf in orfs:
				if sub_seq != None:
					sub_seq = sub_seq.upper().replace(' ','')
					if sub_seq in orf: 
						best_phi=phi
						best_ORF=orf
						found_sub = True	
						break
					
				else:
					if len(orf) > len(best_ORF):
						best_ORF = orf
						best_phi = phi
	
				
		self.number_offset, self.phase = divmod(best_phi,3)
		self.trips=split_codons(self.seq, best_phi)
		return (best_ORF, found_sub)

if __name__ == '__main__':
	seq='tatacatatgagcaatgaagacgacatgaaaaaactgtataaacaaatggtgcaggaactggaaaaagcccgtgaccgtatggaaaaactgtataaagaaatggtggaactgattcagaaagcaatcgaactgatgcgtaaaatctttcaggaagtgaaacaagaagttgaaaaagcaatcgaagaaatgaaaaaactgtacgatgaagcgaaaaagaaaattgaacagatgatccagcaaatcaaacagggcggtgacaaacaaaaaatggaagaactgctgaaacgcgcgaaagaagaaatgaaaaaagttaaagataaaatggaaaaactgctggaaaaactgaaacagatcatgcaagaagccaaacaaaagatggagaaactgctgaaacagctgaaagaagaaatgaagaaaatgaaagaaaagatggaaaaactgttaaaagaaatgaaacagcgtatggaagaagtgaaaaagaaaatggatggcgatgacgaactgctggaaaaaattaagaaaaacatcgatgacctgaagaaaattgcggaagacctgatcaaaaaagccgaagaaaacatcaaagaagcgaagaaaattgccgaacagctggtcaaacgcgcgaaacaactgattgaaaaagcaaaacaggtggctgaagaactgatcaagaaaatcctgcaactgatcgaaaaagccaaagaaatcgccgaaaaagtgctgaaaggcactagtgacgg'
	
	g=DNASeq(seq,1)
	g.insert_restrict( "PstI", 12, 0 )
	g.insert_restrict( "HindIII", 5, 0 )
	g.insert_restrict("NheI", 235, 2)
	g.insert_restrict( "XhoI", 243, 0 )		

	print g.fancy_print()
	print g.seq
