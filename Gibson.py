import src.FragLib as libFrag
import sys
from optparse import OptionParser
from string import split

def register_options():
	parser = OptionParser()
	
	parser.add_option("-f", "--file", dest="file", default=None)
	parser.add_option("--auto",action="store_true", default = False)	
	#Little hackey...To read in list of mutations, mutations are first collected in the leftover args list
	#The args list are then added back to the parser as a defult option
	(opt, args) = parser.parse_args()
	parser.add_option("--DONNOTCALLMECOMMANDLINE", dest="muts", default = args  )
	(opt, args) = parser.parse_args()

	
	if opt.file == None:  sys.exit( "A fragment libray is required. '-f /PATH/TO/LIBRARY'")

	return opt

opt = register_options()

FragLib = libFrag.read_oligo_file_return_frag_lib(opt.file)

#if ( opt.auto ): FragLib.DNAobj.auto_phase()

for mut in opt.muts:
	try: 
		FragLib.mutate(mut)
	except ValueError: 
		sys.exit( mut+" is invalid mutation target. (e.g. Proper format A.1.D)")
		
		

print FragLib.DNAobj.fancy_print()
print 10*'*'+'DNA FRAGMENTS TO ORDER' +10*'*'
for frag in FragLib.frags:
	data = split(frag.__str__())
	print data[0]+data[1]+data[2]+';\t'+data[3]
FragLib.algn_file_out()

