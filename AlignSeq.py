import src.NW_aln as NW
import src.DNASeq as libDNA
import src.AASeq  as libAA
import sys
from optparse import OptionParser


def get_user_input():
	parser = OptionParser()

	(options, args) = parser.parse_args()	
	#parser.add_option("--seq", default=args)
	#(options, args) = parser.parse_args()	
	if args != []:
		if len(args) != 2:
			a=''
			for one in args: a+= 2*'\n'+one  
			sys.exit("Invalid sequence arguments. "+ a+'\n')  
		options.seq = args
	if args == []:
		seq1 = raw_input("Enter sequence 1: ")
		seq2 = raw_input("Enter sequence 2: ")
		options.seq = [seq1,seq2]
	return options

def parseline( seq, lw=78 ):
	line=[]
	num_lines,remain= divmod( len(seq), lw )
	if remain != 0: num_lines+=1
	
	for i in range(0,num_lines):
		if i == num_lines: lw=remain
		line.append( seq[i*lw:(i+1)*lw] )
	return line
	


opt = get_user_input()


obj1,obj2, my_matrix = None, None, None 
try:
	obj1=libDNA.DNASeq( opt.seq[0] )
	obj2=libDNA.DNASeq( opt.seq[1] )
	my_matrix = NW.DNAScore()

except libDNA.InvalidDNASeq:
	try:
		obj1 = libAA.AASeq( opt.seq[0] )
		obj2 = libAA.AASeq( opt.seq[1] )
		my_matrix =  NW.ProtScore()
	except libAA.InvalidAASeq:
		print 'Invalid Sequence'
		exit()  

[s1,s2,s3] =  NW.NW( obj1.seq, obj2.seq, matrix = my_matrix )
s1=parseline(s1)
s2=parseline(s2)
s3=parseline(s3)

dataline='\n\n'

if len(s1) < 2: 
	print s1[0]
	print s3[0]
	print s2[0]
else:
	for i, line in enumerate(s1):
		dataline += s1[i] +'\n' + s3[i] +'\n'+ s2[i] +'\n\n'

	print dataline


#if __name__ == '__main__':
#	seq1="ATGFSWEQPLNHIKDFSCRRRTTPQNCY"
#	seq2="ATGFSWEQPFCHINKDLNHIKDFSCKRKTSPQNCY"
#	[s1,s2,s3] =  NW.NW( seq1, seq2, matrix = NW.ProtScore() )
#	print s1,'\n',s3,'\n',s2
