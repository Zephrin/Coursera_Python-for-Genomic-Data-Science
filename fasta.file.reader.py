#!/usr/bin/python

def seqlength(seqin):
	"will determine seq length"
	
	

try:
	f = open("dna1.fasta")
except IOError:
	print("File not found")
seqs={}
for line in f: #reads entire fasta file, adds all sequences to a dictionary named by the first word in the header
	line=line.rstrip() #removes /n  in line
	if line[0]=='>': #checks if the line begins with > character
		words=line.split() #makes a list of the words in the header as elements
		name=words[0][1:] #adds the first word, minues the > character
		seqs[name]='' #initializes key for the value of "name"
	else: #for seq lines
		seqs[name] = seqs[name]+line #adds new sequence line to the current line
f.close()		#closes the fasta file once you've added all dictionary entries with sequences

#initilize saved sequence variables
longestgene=['no genes tested']
longestseq=''
genenum = 0

print ('the gene lengths are:')
for gene,sequence in seqs.items(): 
	genenum+=1
	#print(name)
	#determine longest seq
	print(len(sequence))
	if len(sequence) > len(longestseq): #will find longest sequence
		del longestgene[0:] #delete all shorter genes
		longestgene.append(gene)		#starts new gene list
		longestseq = sequence #stores longest sequence length
	elif len(sequence)==len(longestseq):
		longestgene.append(gene) #adds tied gene length to list
		
#longest gene printout		
print('total genes =', genenum)
print('the number of longest genes is:', len(longestgene))		
print('the longest gene is ', len(longestseq), 'bp')
print('the longest genes are', longestgene)
#print('the longest seq is', longestseq)

#initialize variables
shortestseqlen=999999999999999999999999999999999
shortestgene=['none']
shortestseq=['none']

for gene,sequence in seqs.items(): 
	#determine shortest seq
	if len(sequence) < shortestseqlen: #will find shortest sequence
		del shortestgene[0:] #delete all loonger genes
		del shortestseq[0:] #delete all loonger genes
		shortestgene.append(gene)		#starts new gene list
		shortestseq.append(sequence) #starts new seq list
		shortestseqlen = len(sequence) #stores shortest sequence length
	elif len(sequence)==shortestseqlen:
		shortestgene.append(gene)		#adds gene
		shortestseq.append(sequence) #adds seq
	
#shortest gene printout		
print('shortest gene is:', shortestseqlen,'bp')		
print('the shortest genes are', shortestgene)
#print('the shortest gene seqs are',shortestseq)

#####ORF identifier

def stopcodon(dna) :
	"""This function checks if given dna sequence has in frame stop codons from position 0."""
	stop_codon=0
	stop_codons=['taa','tag','tga']
	for i in range(0,len(dna),3) :
		codon=dna[i:i+3].lower()
		if codon in stop_codons :
			return i			
	return stop_codon

def startcodon(dna,startsearch):
	"""this fucntion should find start codons in a specified reading frame"""
	for i in range(startsearch,len(dna),3):
		codon=dna[i:i+3]
		if codon == 'atg':
			return i
			
	return startsearch
	
ORFS={} 
ORFpos={}
ORFlen={}
for gene,sequence in seqs.items():
	pos_atg=0 ##CHANGE THIS TO CHANGE FRAMES 
	#0 for frame 1
	#1 for frame 2
	#2 for frame 3
		
	longORFlen=0
	lowseq=sequence.lower() #converts all character to lower case
	num_atg=lowseq.count("atg") #counts num start  in fwd frames 1,2,3
	
	for i in range(num_atg):
		pos_atg=startcodon(lowseq,pos_atg)#finds first/next atg in this frame #REVISED
		#pos_atg=lowseq.find('atg',pos_atg) #finds first/next atg in all frames #ORIGINAL
		newORFlen=(stopcodon(lowseq[pos_atg:]))+3
		if newORFlen>longORFlen: #tests if current ORF is the longest found so far
			#print(pos_atg)
			longORFlen=newORFlen #records longest ORF length
			ORFS[gene]=lowseq[int(pos_atg):int(pos_atg+longORFlen*3)] #records longest ORF seq for each fasta seq
			ORFpos[gene]=pos_atg+1 #records starting pos of longest ORF for each fasta seq
			ORFlen[gene]=longORFlen #dict with length of longest orf for each gene
		
		pos_atg=pos_atg+3 #moves pos_atg to next codon for next search
	
#determines longest gene in fasta file
longestgene=0
longestgenename=''
for gene,length in ORFlen.items():
	if length>longestgene:
		longestgenename=gene
		longestgene=length
		
print('The longest ORF in the fasta file is from fasta sequence', longestgenename)
print('the longest ORF length is', longestgene)		
	
#print(ORFpos)
#print(ORFlen)	
#print(ORFlen["gi|142022655|gb|EQ086233.1|323"])
	
###REPEATS
repeats={} 
repeatn=12  #set repeat search length
for gene,sequence in seqs.items(): #loops through all seq in dictionary "seqs
	for i in range(len(sequence)-(repeatn-1)): #loops through every character of the sequence
		sear_term=sequence[i:i+repeatn] #sets the search term with length n 
		if sear_term in repeats : #tests if search term is in repeat dictionary
			repeats[sear_term]=repeats[sear_term]+1 #if so, adds one to counter
		else:
			repeats[sear_term]=1 #begins counter for new search terms
		
		#for genesear,seqsear in seqs.items(): #loops through all dict again
		#	for i in range(len(sequence)): #scans whole length of sequence

#print(sorted(list(repeats.values())))	#prints ordered list of keys to show how many times each term was found	
#print('cggcgct',repeats['CGGCGCT'])
#print('cggcggc',repeats['CGGCGGC'])
#print('gccgccg',repeats['GCCGCCG'])
#print('tggtggc',repeats['TGGTGGC'])
print ('The most repeated sequence in the FASTA file of length',repeatn,'is :')		
print(max(repeats))
