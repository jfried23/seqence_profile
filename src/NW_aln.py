from string import lower, strip, split
import os

dir = os.path.dirname(__file__)
filename = os.path.join(dir, '../ref/')

def ReadScoreMatrix( file_path ):
	""" Reads in a PSSM scoreing matrix and returns a dictionary mapping all to all """
	_dict={}
	_key=[]
		
	for line in open(file_path):
		if strip(line)[0] == '#': continue
		data = split(line)
		if _key == []:
			_key = data
			for one in _key:
				_dict[one]={}
			continue

		for i,elem in enumerate(data[1:]):
			_dict[ data[0] ][ _key[i] ] = int(elem) 
	return _dict	

def DNAScore():
	return ReadScoreMatrix(filename+'DNAscore')
def ProtScore():
	return ReadScoreMatrix( filename+'BLOSUM62')

def NW(string1, string2, open_gap=-100, ext_gap = -3, matrix=None):
	"""Performs Needleman-Wunsch alignment of string1 and string2.
    	Prints out the alignment and returns the array of scores and pointers(arrows).

    	Example usage from an interactive shell:
        from NeedlemanWunsch import NW
        Scores, Pointers = NW('PELICAN','COELACANTH')

   	 This is modified from a Perl implementation in the book BLAST by Korf, et al.
   	 """
	string1=string1.replace('\r','').replace('\n','')
	string2=string2.replace('\r','').replace('\n','')

	# initialize scoring and 'arrow' matrices to 0
    	Scores = [[0 for x in range(len(string2)+1)] for y in range(len(string1)+1)]
    	Pointers = [[0 for x in range(len(string2)+1)] for y in range(len(string1)+1)]

    	# initialize borders
    	# for pointers (arrows), use 2 for diagonal, -1 for horizontal, and 1 for vertical moves (an arbitrary system).
    	# I have tried to consistently use i for rows (vertical positions) in the score and pointer tables, and j for columns (horizontal positions).
    	for i in range(len(string1)+1):
        	Scores[i][0] = 0  #gap*i No End Gap penalty
        	Pointers[i][0] = 1
    	for j in range(len(string2)+1):
        	Scores[0][j] = 0  #gap*j No End Gap penalty
        	Pointers[0][j] = -1

    	# fill with scores
   	for i in range(1,len(string1)+1):
        	for j in range(1,len(string2)+1):
            		letter1 = string1[i-1]
            		letter2 = string2[j-1]
			pts = matrix[letter1.upper()][letter2.upper()]
            		DiagonalScore = Scores[i-1][j-1] + pts

	    		if Pointers[0][j-1] != 2: pen = ext_gap
	    		else: pen = open_gap
	
            		HorizontalScore = Scores[i][j-1] + pen
            		UpScore = Scores[i-1][j] + pen
            		# TempScores is list of the three scores and their pointers
            		TempScores = [[DiagonalScore,2],[HorizontalScore,-1],[UpScore,1]]
            		# Now we keep the highest score, and the associated direction (pointer)
           		Scores[i][j], Pointers[i][j] = max(TempScores)

		# backtrace from the last entry.  
   	[i,j] = [len(string1),len(string2)]
    	align1 = ''
    	align2 = ''
	malign = ''
	while [i,j] != [0,0]:
        	if Pointers[i][j] == 2:
           		align1 = align1 + string1[i-1]
            		align2 = align2 + string2[j-1]
			if string1[i-1].upper() == string2[j-1].upper(): malign+=':'
			elif (matrix[string1[i-1].upper()][string2[j-1].upper()] > 0):malign+='.'
			else: malign+=' '
            		i = i - 1
           		j = j - 1
        	elif Pointers[i][j] == -1:
            		align1 = align1 + '-'
            		align2 = align2 + string2[j-1]
			malign+=' '
            		j = j - 1
        	else:
            		align1 = align1 + string1[i-1]
            		align2 = align2 + '-'
			malign+=' '
            		i = i - 1

    	# the alignments have been created backwards, so we need to reverse them:
    	align1 = align1[::-1]
    	align2 = align2[::-1]
	malign = malign[::-1]
    	# in case you want to look at the scores and pointers, the function returns them
	return [align1,align2,malign]

