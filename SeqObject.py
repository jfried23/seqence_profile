import src.DNASeq as libDNA
import src.AASeq  as libAA
from optparse import OptionParser
from os.path import isfile
from os import system
from sys import exit
from string import split, rjust, ljust,center,replace

def get_user_input():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="file", default=None)
	parser.add_option("-s", "--sequence", dest="seq")
	parser.add_option("--out", dest="out", default=True, action="store_false")
	parser.add_option("-p", "--phase", dest="phase", default=None)
	parser.add_option("-a","--analyze", dest="ana", default=True, action="store_false" )
	parser.add_option("--sub", dest="sub",default=None)
	#parser.add_option("-Genewiz", dest="genewiz_file", default=None)
	(options, args) = parser.parse_args()	
	if options.file and not isfile( options.file ):
		exit("Path to sequence file is invalid!")
 
	return options

def analyze_protein( opt ):
	AAobj = libAA.AASeq( opt.seq )
	print '\nProtein molar extenction coefficent: '+str(round(AAobj.E280()/1000.,2))+" 1/mM*L."
	print 'Protein molecular weight: '+ str(round(AAobj.MW()/1000.,2))+" kDa."

def analyze_dna( opt ):
	dna=None
	aaSeq=('',False)
	if opt.sub:
		print 'Doing sub'
		dna= libDNA.DNASeq( opt.seq )
		aaSeq = dna.auto_phase( opt.sub )
	elif opt.phase:
		 dna= libDNA.DNASeq( opt.seq, int(opt.phase) )
	else:
		dna=libDNA.DNASeq( opt.seq )
		aaSeq = dna.auto_phase()
	
	aaseq=dna.translate().replace(' ','')[dna.number_offset:]
	print ">Input sequence"		
	print dna.seq+'\n'
	print ">Digested sequence"
	print dna.print_digest()+'\n'
	print ">Translated to AA"
	print aaseq 
	print '**************************************************************************\n'
	print dna.fancy_print()
	
	aaobjc =  libAA.AASeq( aaSeq[0] )
	print ">Best AA Open Reading Frame"
	print aaobjc
	print '\nProtein molar extenction coefficent: '+str(round(aaobjc.E280()/1000.,2))+" 1/mM*L."
	print 'Protein molecular weight: '+ str(round(aaobjc.MW()/1000.,2))+" kDa."

def read_fasta( opt ):
	pass	

opt=get_user_input()

try:
	analyze_dna(opt)

except libDNA.InvalidDNASeq:
	try:
		analyze_protein(opt)
	except libAA.InvalidAASeq:
		print 'Invalid sequence'
		exit()






"""
need_clear = True

dna=''
try:
	dna=DNASeq( opt.seq, int(opt.phase) )
except InvalidDNASeq: pass

	

if opt.out:
	seq= dna.seq	
	splitstr=''
	digest=dna.digest()
	
	ii=0
	for site in sorted(digest.keys()):
		interject = ' <--'+digest[site]+'--> '
		splitstr += seq[ii:site]+interject
		ii=site
	splitstr += seq[ii:]

	print ">Input sequence"		
	print dna.seq+'\n'
	print ">Digested sequence"
	print splitstr+'\n'
	print ">Translated to AA"
	print dna.translate().replace(' ','')[dna.number_offset:]	
	
if opt.ana:
	if need_clear:
		system('clear')
		need_clear = False
	print '*********************************************************************'
	print '**********              sequence ANALYSIS                *************'
	print '*********************************************************************\n'
	print dna.fancy_print()
	print '*********************************************************************'
	print '*********************************************************************\n'


"""
