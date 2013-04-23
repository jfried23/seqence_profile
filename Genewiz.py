import lib.DNASeq as libDNA
import lib.AASeq  as libAA
import string, sys, random, zipfile, os
from string import lower
from optparse import OptionParser

optionparser = OptionParser()
optionparser.add_option( '-f', '--file_path', help='path to seqencing zip file' )
optionparser.add_option( '--sub', default=None )


(opt,args) = optionparser.parse_args()

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
	return (name,seq.lower())

fixBadZipfile(opt.file_path)
zf = zipfile.ZipFile(opt.file_path, 'r')
files= zf.namelist()

out_name=opt.file_path.split("/")[-1].split('.')[0]

dna=open("./"+out_name+'_dna.fasta','w')
out=open("./"+out_name+'_aa.fasta','w')

print "./"+out_name+'_aa.fasta'
for f in files:
	if ".seq" not in f: continue
	name,seq=read_assembly(zf.open(f,'r').readlines())
	if 'Term' in f: seq= libDNA.reverse_compliment(seq)
		
	if len(seq) < 20: continue	

	DnaSeqObj = libDNA.DNASeq(seq)
	rf, found = DnaSeqObj.auto_phase(opt.sub)

	if found:
		aaobj = libAA.AASeq( DnaSeqObj.translate().replace(' ','') )
		out.write(">"+name+ aaobj.seq +"\n\n")
	else:
		out.write(">"+name+ DnaSeqObj.translate().replace(' ','') +"\n\n")	
	dna.write(">"+name+DnaSeqObj.seq+"\n\n") 
		
	continue
	
out.close()
dna.close()
 
