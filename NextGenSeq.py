import src.DNASeq as libDNA
import src.AASeq  as libAA
import src.NW_aln as NW

import string, sys, random, zipfile, os,re
from string import split, ascii_lowercase
from glob import glob
from string import lower
from optparse import OptionParser

optionparser = OptionParser()
optionparser.add_option( '-s', '--sequencing_data', help='path to sequencing zip file' )
optionparser.add_option( '-r', '--ref_data', help='path to reference file' )
optionparser.add_option( '-n', '--num_char', default=50, help='path to reference file' )
optionparser.add_option( '--search', default=True, action="store_false", help='Search all the reference sequences against all the sequencing results.' )



(opt,args) = optionparser.parse_args()

###############Utility Functions###############################################
def longest_sub_seq_without( seq, char):
	'''Returns the longest subseqence in seq without an instance of char
	'''
	ns = [m.start() for m in re.finditer(char, seq)]
	s0,s1=0,0
	i=0

	if len(ns) == 0: return seq
	if len(ns) == 1:
		p1 = seq[0:ns[0]]
		p2 = seq[ns[0]+1:]
		if len(p1) > len(p2): return p1
		else: return p2
 
	while i < len(ns)-1:
		if ( (ns[i+1]-ns[i]+1) > (s1-s0) ): s0,s1=ns[i]+1,ns[i+1]
		i+=1
	
	return seq[s0:s1]

def fixBadZipfile(zipFile):  
     f = open(zipFile, 'r+b')  
     data = f.read()  
     pos = data.find('\x50\x4b\x05\x06') # End of central directory signature  
     if (pos > 0):  
         f.seek(pos + 22)   # size of 'ZIP end of central directory record' 
         f.truncate()  
         f.close() 

def read_assembly(lines):
	name=""
	seq=""
	for line in lines:
		if ">" in line: 
			name=line[1:]
			continue
		seq+=line.rstrip()
	return (name,seq)

def Read_sequencing( path ):
	seq_files=[]
	if os.path.isdir(path):
		seq_files = glob(path+'*.zip')
	else: seq_files.append( path )

	DNAseq_objs={}
	for File in seq_files:
		try:
			fixBadZipfile( File )
			zf = zipfile.ZipFile(File, 'r')
			for f in zf.namelist():
				if '.seq' not in f: continue
				name,seq = read_assembly(zf.open(f,'r').readlines())
				name = name.split('.')[0]
				if 'Term' in f: seq = libDNA.reverse_compliment(seq.lower())
				DNAseq_objs[name]= seq.lower().replace(' ','')					
		except: 
			print "Can not open "+File
			continue

	return DNAseq_objs

def Read_RefDataFASTA( path ):
	RefData = {}
	lines = open(path,'r').readlines()
	
	for i,line in enumerate(lines):
		if '>' not in line: continue
		
		name = lines[i]
		RefData[ name[1:].replace('\n','').replace('\r','') ] = lines[i+1].lower().replace('\n','').replace('\r','').replace(' ','')
	return RefData


	
seq_data = Read_sequencing( opt.seqencing_data )	
ref_data = Read_RefDataFASTA( opt.ref_data )

gbl_out = open('data.txt','w')
gbl_out.write( 'vi:set nowrap:\n' )

files={}

if not os.path.isdir('./data'): os.makedirs('data')
	
dnaScoreMatrix = NW.DNAScore()
prtScoreMatrix = NW.ProtScore()

if opt.search:
	for key, exp_seq in seq_data.iteritems():
		
		for key2, seq in ref_data.iteritems():
				
			i,frag_size=0,40
			frags=[]

			while i < len(seq):	
				fr = seq[i:i+frag_size]
				if len(fr) == frag_size: frags.append(fr)
				i+=frag_size
			
		
			passTest=0
			for frag in frags:
				if frag in exp_seq:
					strt =  exp_seq.find(frag)
					passTest+=1
			
			if len(frags) ==0 : continue
	
			if ( (1.*passTest/len(frags) > 0.5 ) ):
				[align1, align2, malign] = NW.NW( exp_seq, seq, matrix=dnaScoreMatrix )

				if  ':'*opt.num_char not in malign: continue
				print "Match    "+key+" <===> " + key2 				
	
				i=0
				while align2[i] not in ascii_lowercase: i+=1

				seqobj1 = libDNA.DNASeq( exp_seq )
				refobj1 = libDNA.DNASeq( seq )
				
				best_ORF1, found_sub = seqobj1.auto_phase('HHHH')  #Check for his tag
				if not found_sub: best_ORF1, found_sub = seqobj1.auto_phase()    #Else take longest ORF
				
				best_ORF2, found_sub = refobj1.auto_phase('HHHH')  #Check for his tag
				if not found_sub: best_ORF2, found_sub = refobj1.auto_phase()    #Else take longest ORF
				
				[aalign1, aalign2, aamalign] = NW.NW( best_ORF1 , best_ORF2, matrix=prtScoreMatrix )
			
				if key2 not in files.keys(): 
					files[key2] = open('./data/'+key2+'.txt','w')
					files[key2].write( 'vi:set nowrap:\n' ) 

				keyLen = max( len(key), len(key2) )+2
				pkey1=string.ljust(key,keyLen)
				pkey2=string.ljust(key2,keyLen)
				
				data  = ''
				data += seqobj1.fancy_print(len(seqobj1.seq)+5)+'\n\n'

				data += pkey1+aalign1+'\n'
				data += string.ljust(' ',keyLen)+aamalign+'\n' 
				data += pkey2+aalign2+'\n\n'
	
				data += pkey1+align1[i:]+'\n'
				data += string.ljust(' ',keyLen)+malign[i:]+'\n' 
				data += pkey2+align2[i:]+'\n'+'*'*(len(exp_seq)+5) +'\n\n\n\n'
				
		
				files[key2].write( data )
				gbl_out.write( data )

