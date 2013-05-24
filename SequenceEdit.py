import src.DNASeq as libDNA
import src.AASeq  as libAA
from optparse import OptionParser

import ref.dna_data as ref

def get_user_input():
	parser = OptionParser()
	parser.add_option("-s", "--sequence", dest="seq")
	parser.add_option("-p", dest="phase", default=0)
	parser.add_option("-x", "--cutsite", dest="cutsite", default=1)

	(options, args) = parser.parse_args()	

	return options

opt = get_user_input()

def number_differences( nativeS, s1):
	i,count = 0,0
	while i < len(nativeS):
		if nativeS[i] != s1[i]: count+=1
		i+=1
	return count
	


dna=libDNA.DNASeq( opt.seq, int(opt.phase) )

start_pos = dna.fetch_position( int(opt.cutsite) - 3 )
end_pos = dna.fetch_position( int(opt.cutsite) + 3 )

sub=dna.seq[start_pos:end_pos+3]
h=libDNA.DNASeq(sub)
nativeAA=libDNA.translate( h.trips )
print nativeAA 

slots=[]

for a in nativeAA.replace(' ',''): slots.append([])


for enzyme in ref.restriction.keys():
	cut_seq = ref.restriction[ enzyme ][0]	
	for ph in range(0,3):
		tmp=libDNA.DNASeq(dna.seq, int(opt.phase) ) 
		tmp.insert_restrict( enzyme, int(opt.cutsite), int(ph) )
		
		sub=tmp.seq[start_pos:end_pos+3]
		h=libDNA.DNASeq(sub)
		seq = libDNA.translate( h.trips )
		if 'X' in seq: continue 
		diffs = number_differences(nativeAA,seq)
		
		slots[diffs].append( (seq, enzyme, ph) )

for one in slots:
	if one == []: continue
	for site in one:
		print site[0], site[1], site[2] 
