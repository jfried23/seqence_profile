from string import split, ljust, rjust, lower, center, join
import sys
from DNASeq import DNASeq, compliment,reverse_compliment,split_codons,translate, WrongAA
from CodonData import CodonData

def print_frag( frag ):
	cod=[]
	keys=[]
	for key in sorted(frag.coding_codon.keys() ):
		cod.append( frag.coding_codon[key] )
		keys.append( center(str(key),3) )
	print keys
	print cod
	print frag.seq()

def compare_codons( c1, c2 ):
	c1 = c1.replace(' ','')
	c2 = c2.replace(' ','')
	if (c1 in c2) or (c2 in c1): return True
	else: return False

def LastIndex(alist, value):
	""" returns the last occuernce of an value in a list """
	return len(alist) - alist[-1::-1].index(value) -1


def read_oligo_file_return_frag_lib(path):
	start_pos=0
	frags=[]
	for line in open(path):
		data=split(line)
		if 'startpos' in line: 
			start_pos=int(data[1])
			continue
		if '****' in line: break
		frags.append( Fragment(data[0],data[1],data[2],data[3])  )

	FL=FragLib( frags, start_pos, path)

	return FL


class Fragment( object ):
	""" Fwd fragments stored in 5'->3' polarity; Rev fragments 3'->5' """
	def __init__(self, name, polarity, number, seq ):

		assert( polarity=='+' or polarity=='-' )	
		self.updated    = False
		self.name       = name
		self.number     = number
		self.__seq      = seq.lower().replace(' ','')
		self.__polarity = polarity	
		self.overlap_data =[ None, None ]
		if self.is_rev(): self.__seq = self.__seq[::-1]
		
		self.coding_codon={}

	def __str__(self):
		s=''
		if self.is_rev(): s+= self.__seq[::-1]
		else: s=self.__seq	
		return ljust(self.name,20)+'  '+ self.__polarity +'  '+self.number +'  '+ s

	def __len__(self): return len(self.seq().replace(' ',''))

	def overlaps(self, frag2, minOver=8 ):
		"""Find if one Fragment overlaps another Fragment object returns string indceces that overlap
		   first value is overlap of self, second value overlap of its binding partner """
		if self.is_fwd() == frag2.is_fwd(): return False 
		pos=-1
		
		least_len=max( len(self), len(frag2) )

		pos = self.seq().find(frag2.seq()[0:minOver] )
		if pos != -1:
			self.overlap_data [1] =  frag2 
			frag2.overlap_data[0] =  self 
			return True
	
		pos = frag2.seq().find( self.seq()[0:minOver] )
		if pos != -1:
			self.overlap_data [0] = frag2 
			frag2.overlap_data[1] = self
			return True

		if pos == -1: return False

	def is_fwd(self): return self.__polarity=='+'
	def is_rev(self): return self.__polarity=='-'

	def contains( self, AAindex): return (int(AAindex) in self.coding_codon.keys() )
	
	def codons( self ):
		vl=[]
		for key in sorted(self.coding_codon.keys()):
			vl.append( self.coding_codon[key] )
		return key, vl

	def mutate(self, index, codon):
		if index not in  self.coding_codon.keys(): return
		cod = self.coding_codon[index]
		if cod[0] == ' ':
			ct = cod.count(' ')
			self.coding_codon[index] = codon[ct:]
		elif cod[-1] == ' ':
			ct = cod.count(' ')
			self.coding_codon[index] = codon[0:(3-ct)]
		else: self.coding_codon[index] = codon
		self.__seq_from_codons()

	def seq(self, coding=True):
		"""Returns the sequence of the fragment in the 'sence' sence """
		if coding:
			if self.is_rev(): return compliment(self.__seq)
			else: return self.__seq
		else:
			return self.__seq

	def __seq_from_codons(self):
		t_seq=''
		for one in sorted(self.coding_codon.keys()): 
			t_seq+=self.coding_codon[one]
		t_seq = t_seq.replace(' ','')
		if self.is_rev():
			t_seq = compliment( t_seq )
		self.__seq = t_seq
		self.updated=True

class FragLib( object ):
	def __init__( self, frags, start_pos, path ):
		self.frags = frags
		self.start_pos = start_pos
		self.__spos, self.__ph = divmod(start_pos,3)
		self.DNAobj = None
		self.CodonData = CodonData()	
		self.min_overlap=8
		
		self.__order_frags()
		self.assemble( self.min_overlap )

		self.__path = path
		self.__made_mods=[]
	

	def __order_frags(self):
		""" Automoatically sorts the fragment list into the proper assembly order
		and ensures all the fragments have a logical place
		"""
		for frag1 in self.frags:
			for frag2 in self.frags:
				frag1.overlaps(frag2, self.min_overlap)	
		head=None
		tail=None
		for frag in self.frags:
			if frag.overlap_data[0] == None:
				assert( head == None )
				head = frag
			if frag.overlap_data[1] == None:
				assert( tail == None )
				tail = frag
		ptr=head
		tmp_frags=[]
		while (ptr != None):
			tmp_frags.append(ptr)
			ptr = ptr.overlap_data[1]

	def algn_file_out(self): #Assumes the fragments were properly shuffled by __order_frags
		pos,pos_n,neg, neg_n='','','',''
		seq = self.DNAobj.seq
		for frag in self.frags:
			st_pos = seq.find( frag.seq() )
			this_seq = frag.seq(False)
			
			if frag.is_fwd():
				diff = st_pos - len(pos)
				name_lbl = '|-> '+frag.name + ' +'+str(frag.number)
				pos_n += diff*' ' + name_lbl + (' '* ( len(frag)- len(name_lbl) ) )
				pos   += diff*' ' + this_seq
			if frag.is_rev():
				diff = st_pos - len(neg)
				name_lbl = '|<- ' + frag.name +' -' + str(frag.number)
				neg += diff*' ' + this_seq
				neg_n += diff*' ' + name_lbl + (' '* ( len(frag)- len(name_lbl) ) )				
		
		mut_info=''
		for i, one in enumerate(self.__made_mods):
			if ( i== 0 ): mut_info += "Mutations:   "
			mut_info += one+' '

		return_data = "vi:set nowrap:\n\n"  +'\n'	
		return_data += self.DNAobj.fancy_print( lw=len(seq)+3 )
		return_data += pos_n+'\n'+pos+'\n'+neg+'\n'+neg_n+'\n'	
		return_data += self.__str__()
		
		h=open( self.__path+'.assembly','w')
		h.writelines( return_data )
		h.close()

	def __str__(self):
		ret_str='startpos '+str( 3*self.__spos + self.__ph ) +'\n'
		altered=[]
		for frag in self.frags:
			if frag.updated: altered.append(frag)
			ret_str += ( frag.__str__() + '\n' )

		if len(altered) != 0: 
			ret_str += '\n\n'+20*'*'+'NEW FRAGMENTS'+20*'*'+'\n'
			for one in altered:
				ret_str += one.__str__() +'\n'
		return ret_str	
			
	def output(self):
		fwd_s, rev_s, tot_s = '','', self.DNAobj.seq
		
		data = self.DNAobj.fancy_print( lw = len(tot_s)+ 5 )

		altered,al_str=[], ''
		for frag in self.frags: 
			align =  tot_s.find( frag.seq() )
			if frag.is_fwd(): 
				pad = ' '*( align - len(fwd_s) )
				fwd_s += ( pad + frag.seq() )
			else:
				pad = ' '*( align - len(rev_s) )
				rev_s += ( pad + frag._Fragment__seq )

			if frag.updated: altered.append(frag)	

		return 'vi:set nowrap:'+'\n'+ data+fwd_s+'\n'+rev_s+'\n\n\n'+self.__str__()	

	def assemble( self, min_overlap = 6 ):
		seq=''
		offset=0
		for i,frag in enumerate(self.frags):
			frag.coding_codon={}
			if i == 0:
				seq += frag.seq().replace(' ','')
				codons = split_codons( frag.seq(), self.__ph )

				##find where to start counting codons:
				offset = self.__spos
				if ( self.__ph != 0 ): offset += 1
				##End counting logic block

				for indx, codon in enumerate(codons):
					frag.coding_codon[ 1 + indx - offset ]=codon
				continue
	
			overlap = seq.find( frag.seq()[0:min_overlap] )
			assert( overlap != -1 )
			if (seq[overlap:] != frag.seq()[0:(len(seq)-overlap)] ): ##CHECK THAT WE HAVE PERFECT ALIGNMENT!
				#get required info
				f1obj = self.frags[i-1]
				f2obj = frag
				f1name = f1obj.name + f1obj._Fragment__polarity +f1obj.number
				f2name = f2obj.name + f2obj._Fragment__polarity +f2obj.number
				print "Mismatch between framgnets "
				print ljust(f1name,20), f1obj.seq(False)
				print ljust(f2name,20), f2obj.seq(False)
				sys.exit()

			seq += frag.seq()[(len(seq) - overlap):]

			num,ph = divmod( overlap-self.start_pos, 3)
			ph=3-ph
   			codons = split_codons( frag.seq(), ph )
						
			for ii, codon in enumerate(codons):
				this_index = num + ii + 1 
				frag.coding_codon[ this_index  ] = codon	
					##########Optional block, double check for consistancy#######
				if self.frags[i-1].contains(this_index):
					other = self.frags[i-1].coding_codon[this_index]
					if (' ' in other) or (' ' in codon): continue
					if other != codon: print this_index, other, codon				

		self.DNAobj = DNASeq( seq, self.start_pos )

	def mutate(self, key):
		AA1,num,AA2 = key.split('.')
		
		if AA2 != '': the_codon = self.CodonData.codon(AA2, False)
		else:   the_codon = ''   #Covers deletions

		try: self.DNAobj.fetch_position( int(num)-1, AA1 )  ##Sanity check, ensure AA1 is indeed at pos num
		except WrongAA as e:
			sys.exit("Resiude number "+str(num)+" is "+e.AA+" not an "+AA1+"! Please check mutation input")

		for frag in self.frags:
			if not frag.contains( num ): continue
			frag.mutate(int(num),the_codon)
			if AA2 !='': frag.name+= '_'+AA1+str(num)+AA2
			else: frag.name+= '_x'+AA1+str(num)

		self.assemble( self.min_overlap )
		
		self.__made_mods.append(key)


	def insert(self, pos, AAseq):
		new_dna_seq, new_frag_seq = '',''
		for aa in AAseq: new_dna_seq+= self.CodonData.codon(aa, False)
		for frag in self.frags:
			if not frag.contains(pos): continue
			if ' ' == frag.coding_codon[int(pos)][:-1]: continue
			
			for cod in sorted(frag.coding_codon.keys()):
				new_frag_seq += frag.coding_codon[cod]
				if cod == int(pos): new_frag_seq += new_dna_seq
			if frag.is_rev(): new_frag_seq  = compliment(new_frag_seq.replace(' ','') )
			frag._Fragment__seq = new_frag_seq 
		self.assemble(self.min_overlap )	
			

if __name__ == '__main__':
	path='/Users/Josh/Documents/Baker Projects/ets_domain/append_to_N_HL3/model1/ETS_ex_repeat.oligos'
	FL=read_oligo_file_return_frag_lib(path)
	
	#print FL.DNAobj.fancy_print()
	#FL.mutate('I.35.A')
	#FL.mutate('G.55.W')
	#FL.mutate('L.5.X')
	FL.mutate('Y.64.S')
	
	#FL.mutate('N.65.X')
	#FL.mutate('T.110.P')
	

	#FL.insert(13,'WGGGGGGW')
	
	print FL.DNAobj.fancy_print()
			


