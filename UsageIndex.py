#! /usr/bin/env python

#A. Bayegan
#This program uses forward and backward partition functions obtained from RNAsampleCDS software to 
#calculate codon usage index for non-overlapping as well as overlapping reading frames.
#The minimum input requirement is a fasta file containg the genes required to be analysed and the coding reading frames of the genes. 

import os,sys,math
from RNAsampleCDS import computeZB,computeZF,computeZF_GC,computeZB_GC,createListOfCompatible5mers,NUCL,bfs,parseConstraintFile
from utility import translateMRNA,GCcont
from aminoAcidAndGeneticCodes import genCode,aa2codon,aaCode

DEBUG=0

def computeCodonProbFromGenesOverlap(seqList,rf):
#compute probability of codon w in the input sequence list 'seqList' for the overlapping reading frames given in rf. For each reading frame of rf
#compute the freq(w)/freq(syn(w)) in the corresponding reading frame. 
#Output is dictionary p_gene where keys are codons and values are sum of the probabilities over all reading frames.
	#-------------------initialize the codon frequency dictionary
	freq = {}
	for cod in genCode.keys():
		freq[cod] = 0.0
	#------------------check the length of the sequences to be 3k+2-----------------------------------
	#~ m = n%3
	#~ l = [] 
	#~ if m!=2:
		#~ for seq in seqList:
			#~ l.append(seq[:n-m-1])
		#~ print "WARNING: the last %d nucleotides ignored in computing codon usage from genes!" % (m+1)
		#~ seqList = l
		#~ if l==[]:
			#~ print 'error: sequences are to short!'
			#~ sys.exit(1)
		#~ n = len(seqList[0])
	#------------------calculate frequency of codons in all reading frames --------------------------------
	for i in rf:
		for seq in seqList:
			n = len(seq)
			if i[0]=='+': #positive strands
				r=int(i)
				for idx in range(r,n-2+r,3): #seq[0:n-2], seq[1:n-1], seq[2,n] code in reading frames +0, +1 and +2 respectively
					codon = seq[idx:idx+3]
					if len(codon)==3:
						freq[codon] += 1
			elif i[0]=='-': #negative strands
				r=abs(int(i))
				for idx in range((n-2+r),r,-3): #seq[0:n-2], seq[1:n-1], seq[2,n] code in reading frames +0, +1 and +2 respectively
					codon = seq[idx-3:idx][::-1]
					if len(codon)==3:
						freq[codon] += 1
	if DEBUG:
		print "codon freq in the fasta file:"
		for k in sorted(freq.keys()):
			if freq[k]!=0:
				print k,freq[k]
	#----------------calculate frequency of synonymous codons in all reading frames-------------------------------
	freq_syn={}
	for codon in genCode.keys(): 
		freq_syn[codon] = 0.0
	for codon in freq.keys():
		for syncod in aa2codon[aaCode[genCode[codon].upper()]]:
			freq_syn[codon] += freq[syncod]
	if DEBUG:
		print "total synonymous codon freq in the fasta file: "
		for k in sorted(freq_syn.keys()):
			if freq_syn[k]!=0:
				print k,aa2codon[aaCode[genCode[k].upper()]],freq[k],freq_syn[k]
	#---------------calculate the codon usage index-----------------------------------------------------------------
	p = {}
	for codon in genCode:
		p[codon]=0.0
	for codon in freq.keys():
		if freq_syn[codon]!= 0:
			p[codon] = freq[codon]/freq_syn[codon]
	return p,freq,freq_syn

def computeExpectedCodonProbForSeqs(seqDict,rf,gcDiff,constraintfile=None,filename=None,threshold=1):
#compute probability of codon w for a set of nucleotide sequences
	totalfreq={}
	totalfreqSyn={}
	for codon in genCode.keys():
		totalfreq[codon]=0.0
		totalfreqSyn[codon]=0.0
	p_hyp = {};totCDS=0;
	for codon in genCode.keys():
		p_hyp[codon] = 0.0
	for ids,seq in seqDict.iteritems():
		print ids
		peptides={}
		if(constraintfile):
			constraint = parseConstraintFile(constraintfile,len(seq))
		else:
			constraint = None
		#---------------generate peptides for the given rf ------------
		#translate sequences
		for r in rf:
			if r[0]=='+':
				peptides[r] = translateMRNA(seq[int(r):])
			elif r[0]=='-':
				peptides[r] = translateMRNA(seq[::-1][int(r)+2:])
		n = len(peptides[r])
		
		#if the length of the peptides in all reading frames are not identical make them equal
		N = []
		N.append(n)
		for pep in peptides.values():
			if len(pep)!=n:
				N.append(len(pep))
		n=min(N)
		if len(N)!= 1:
			for k in peptides.keys():
				peptides[k] = peptides[k][:n]
		
		#add the 'O' peptides in reading frames that are not defined in the input
		for i in ['+0','+1','+2','-0','-1','-2']:
			if i not in peptides.keys():
				peptides[i] = 'O'*n
		if DEBUG: print peptides
		print peptides
		print constraint
		#--------------compute pf and frequency dictionary-------------
		if gcDiff==-1: #no gc-content defined
			ZF = computeZF(peptides,constraint,filename,threshold)
			ZB = computeZB(peptides,constraint,filename,threshold)
			ZFval  = 0.0
			for ch3 in NUCL:
				for ch4 in NUCL:
					ZFval += ZF[(n,ch3,ch4)]
			print ZFval
			totCDS += ZFval
			gclow = -1
			gcup = -1
		else: #gc-content defined
			ZF = computeZF_GC(peptides,constraint,filename,threshold)
			ZB = computeZB_GC(peptides,constraint,filename,threshold)
			cnt = round((gcDiff/100.0)*len(seq) )
			gc = GCcont(seq)
			gclow = gc-cnt  
			gcup = gc+cnt
			if (gclow<0): gclow=0
			if (gcup>3*n+2) : gcup=3*n+2
			#~ print len(seq),gc,gclow,gcup
			#~ ZBval_gc  = 0.0
			#~ for ch1 in NUCL:
				#~ for ch2 in NUCL:
					#~ for x in range(3*n+3):
						#~ ZBval_gc += ZF[(n,x,ch1,ch2)]
			#~ #print "gc,len",gc,len(seq)
			#~ print "total:",ZBval_gc
			#~ 
			#~ ZBval  = 0.0
			#~ for ch1 in NUCL:
				#~ for ch2 in NUCL:
					#~ for x in range(3*n+3):
						#~ ZBval += ZB[(0,x,ch1,ch2)]
			#~ print "totalZB:",ZBval
			#~ 
			ZFval  = 0.0
			for ch3 in NUCL:
				for ch4 in NUCL:
					for x in range(int(gclow) , int(gcup+1)):
						ZFval += ZF[(n,x,ch3,ch4)]
			#print ZFval
			totCDS += ZFval
		freq,freq_syn = computeExpectedCodonFreqForPeptide(rf,peptides,ZF,ZB,gclow,gcup,filename,threshold)
		for codon in genCode.keys():
			totalfreq[codon] += freq[codon]
			totalfreqSyn[codon] += freq_syn[codon]
	print "#total number of sequences that code exactly the same peptides as the given fasta file: ",totCDS
	for codon in genCode.keys():
		if totalfreqSyn[codon] != 0:
			p_hyp[codon] = totalfreq[codon] / totalfreqSyn[codon]
	return p_hyp,totalfreq, totalfreqSyn

def computeExpectedCodonFreqForPeptide(rf,peptides,ZF,ZB,gclow,gcup,filename=None,threshold=1):
#compute probability of codon w in the set of all sequences that translate 'peptides' in the corresponding reading frames
	n  = len(peptides['+0'])
	seqLen = 3*n+2
	#----------Initialize freq ---------------
	freq = {}
	for codon in genCode.keys():
		freq[codon] = 0.0
	#-----------calculate frequency of codons in all reading frames of interest -------------------
	###
	if gclow != -1: #if gc-content is defined
		for k in range(n):
			hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			for s in L:
				for x1 in range(0,3*k+3):
					for x2 in range(0,3*(n-k)+3):	
						#print s,k,x1,x2,ZF[(k,x1,s[0],s[1])],ZB[(k+1,x2,s[3],s[4])]
						if( gclow <= x1 + x2 + GCcont(s[2]) <= gcup) : 
							#print '*',k,s,GCcont(s[2]),x1,x2, ZF[(k,x1,s[0],s[1])],ZB[(k+1,x2,s[3],s[4])]
							if '+0' in rf : freq[s[0:3]]+= (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							if '+1' in rf : freq[s[1:4]]+= (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							if '+2' in rf : freq[s[2:5]] += (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							if '-0' in rf : freq[s[0:3][::-1]]+= (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							if '-1' in rf : freq[s[1:4][::-1]]+= (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							if '-2' in rf : freq[s[2:5][::-1]] += (ZF[(k,x1,s[0],s[1])]) * ZB[(k+1,x2,s[3],s[4])]
							#~ #if ZF[(k,x1,s[0],s[1])] * ZB[(k+1,x2,s[3],s[4])]!= 0:
								#~ #print 's:%s\tgc2:%d\tk:%d\tx1:%d\tx2:%d\tZF[k,x1,s[0],s[1]]:%f\tZB[(k+1,x2,s[3],s[4])]:%d\tfreq:%d'%(s,GCcont(s[2]),k,x1,x2, ZF[(k,x1,s[0],s[1])],ZB[(k+1,x2,s[3],s[4])],ZF[(k,x1,s[0],s[1])]*ZB[(k+1,x2,s[3],s[4])])
	else: #if gc-content is not defined
		for k in range(n):
			hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			for s in L:
				if '+0' in rf : freq[s[0:3]]+= (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				if '+1' in rf : freq[s[1:4]]+= (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				if '+2' in rf : freq[s[2:5]] += (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				if '-0' in rf : freq[s[0:3][::-1]]+= (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				if '-1' in rf : freq[s[1:4][::-1]]+= (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				if '-2' in rf : freq[s[2:5][::-1]] += (ZF[(k,s[0],s[1])]) * ZB[(k+1,s[3],s[4])]
				
	###			
				#~ if ZF[(k,s[0],s[1])] * ZB[(k+1,s[3],s[4])]!= 0:
					#~ print k,s,ZF[(k,s[0],s[1])],ZB[(k+1,s[3],s[4])],freq[s[0:3]]
	#----------check frequency matrix (only works without blosum)--------------------------------------------------------
	sum0 = 0.0
	ZFval  = 0.0
	for k in freq.keys():
		sum0 += freq[k]
	if gclow!=-1:
		for ch3 in NUCL:
			for ch4 in NUCL:
				for x in range(int(gclow) , int(gcup)+1):
					ZFval += ZF[(n,x,ch3,ch4)]
		print "####",ZFval
	else:
		for ch1 in NUCL:
			for ch2 in NUCL:
				ZFval += ZF[(n,ch1,ch2)]	
	print (sum0 , ZFval*len(rf)*n)	
	#assert(sum0 == ZFval*len(rf)*n)
	#~ #print (sum0 , ZFval*len(rf)*n)		
	#----------calculate frequency of synonymous codons-----------------------------------------
	freq_syn={}
	for codon in genCode.keys():
		freq_syn[codon] = 0.0
	for codon in freq.keys():
		for syncod in aa2codon[aaCode[genCode[codon].upper()]]:
			freq_syn[codon] += freq[syncod]
	if DEBUG:
		print 'expected codon freq for peptides:', peptides
		for k in sorted(freq_syn.keys()):
			if freq_syn[k]!=0:
				#print k,aa2codon[aaCode[genCode[k].upper()]],freq[k],freq_syn[k]
				print k,freq[k],freq_syn[k]
	#---------------calculate the codon usage index-----------------------------------------------------------------
	#~ p = {}
	#~ for codon in genCode.keys():
		#~ p[codon]=0.0
	#~ for codon in freq.keys():
		#~ if freq_syn[codon]!=0:
			#~ assert((freq[codon]/freq_syn[codon])<=1)
	return freq,freq_syn

def parseFasta(fastafile):
	d={}
	fin = open(fastafile,'r')
	line = fin.readline().strip()
	if not line:
		print "Error in fasta file!"
		sys.exit(1)
	seqid = line[1:]
	line = fin.readline().strip()
	seq=''
	while line:
		if line[0]!='>':
			seq += line
		elif line[0]=='>':
			if seqid in d:
				print 'sequences in the fasta file must have unique ids!'
				sys.exit(1)
			d[seqid]=seq.upper().replace('T','U')
			seq=''
			seqid = line[1:]
		line = fin.readline().strip()
	d[seqid] = seq.upper().replace('T','U')
	return d
 	
def checkCodonUsage(seqList,rf,gcDiff,constraint=None):
	#codon usage computed directly from the sequences must be identical to the codon usage computed from partition functions
	l=0
	RNA0=[]
	for seq in seqList:
		if(constraintfile):
			constraint = parseConstraintFile(constraintfile,len(seq))
		else:
			constraint = None
		p1={}
		for codon in genCode.keys():
			p1[codon]=0.0
		peptides={}
		#translate sequences
		for r in rf:
			if r[0]=='+':
				peptides[r] = translateMRNA(seq[int(r):])
			elif r[0]=='-':
				peptides[r] = translateMRNA(seq[::-1][int(r)+2:])
		
		n = len(peptides[r])
		#if the length of the peptides in all reading frames are not identical make them equal
		N = []
		N.append(n)
		for pep in peptides.values():
			if len(pep)!=n:
				N.append(len(pep))
		n=min(N)
		if len(N)!= 1:
			for k in peptides.keys():
				peptides[k] = peptides[k][:n]
		
		#add the 'O' peptides in reading frames that are not defined in the input
		for i in ['+0','+1','+2','-0','-1','-2']:
			if i not in peptides.keys():
				peptides[i] = 'O'*n
		if DEBUG: print peptides
		#print peptides
		RNA = bfs(peptides,filename,threshold,constraint)
		if gcDiff==-1:
			for rna in RNA:
				RNA0.append(rna)
		#~ for rna in RNA:
			#~ print rna,GCcont(rna)
		else: #gc-content defined
			cnt = round((gcDiff/100.0)*len(seq) )
			gc = GCcont(seq)
			gclow = gc-cnt  
			gcup = gc+cnt
			if (gclow<0): gclow=0
			if (gcup>3*n+2) : gcup=3*n+2
			print gc,gclow,gcup
			seqLen = len(seqList[0])
			for rna in RNA:
				if int(gclow <= GCcont(rna) <=  gcup):
					print rna,GCcont(rna)
					RNA0.append(rna)
	l =  len(RNA0)
	p1,f1,f2 = computeCodonProbFromGenesOverlap(RNA0,rf)
	print '#number of seq from bruteforce: %d' % l
	p2,f3,f4 = computeExpectedCodonProbForSeqs(seqList,rf,gcDiff)
	for codon in genCode.keys():
		print codon,p1[codon],p2[codon]
	assert(p1==p2)

def printusageIndex(exe):
	print "Usage: %s -f fastaFile -r readingFrames -c constraint -m BlosumMatrix -t threshold -gc GC_content" % name
	sys.exit(1)

if __name__=='__main__':
	#------------------default setting------------------
	gcDiff = -1
	threshold=1
	filename=None
	constraintfile=None
	rf=[]
	fastafile=None
	#------------------parse input arguments-------------
	if len(sys.argv) < 3:
		printusageIndex(sys.argv[0])
	args     = sys.argv[1:]
	#fasta file indicated by '-f' flag
	for i in range(len(args)):
		if args[i]=='-f':
			if i+1!=len(args) and args[i+1][0]!='-':
				fastafile = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: peptide file must be secified after -f flag!\n "
				printusageIndex(sys.argv[0])
	#coding reading frames indicated by '-r' flag
	for i in range(len(args)):
		if args[i]=='-r':
			if i+1!=len(args) and args[i+1][0]!='-':
				read = args[i+1]
				for r in read.split(','):
					rf.append(r)
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR:coding reading frames must be specified after -r flag!\n"
				printusageIndex(sys.argv[0])
	#sequence constraint indicated by -c flag
	for i in range(len(args)):
		if args[i]=='-c':
			if i+1!=len(args) and args[i+1][0]!='-':
				constraintfile = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: constraint file must be secified after -c flag!\n "
				printusageIndex(sys.argv[0])
	#similarity matrix filename, with -m flag
	for i in range(len(args)):
		if args[i]=='-m':
			if i+1!=len(args) and args[i+1][0]!='-':
				filename = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else: 
				print "ERROR: matrix file must be specified after -m flag!\n"
				printusageIndex(sys.argv[0])
	#blosum threshold, with -t flag
	for i in range(len(args)):
		if args[i]=='-t':
			if i+1!=len(args) and args[i+1][0]!='-':
				threshold = int(args[i+1])
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: similarity threshold must be specified after -t flag!\n"
				printusageIndex(sys.argv[0])
	#GC content with -gc flag
	for i in range(len(args)):
		if args[i]=='-gc':
			gcDiff = int(args[i+1])
			del[args[i+1]]; del[args[i]]
			break
	#verbose output, with -v flag
	for i in range(len(args)):
		if args[i]=='-v':
			VERBOSE = 1
			del[args[i]]
			break
	
	if rf==[]: #the above check does not work when only -f is defined
		print "coding reading frames must be specified!\n"
		printusageIndex(sys.argv[0])
	if fastafile==None:
		print "fasta file name must be specified!\n"
		printusageIndex(sys.argv[0])
	seqDict = parseFasta(fastafile)
	#checkCodonUsage(seqList,rf,gcDiff,constraintfile)
	p_gene,f_gene,fsyn_gene = computeCodonProbFromGenesOverlap(seqDict.values(),rf)
	p_hyp,f_hyp,fsyn_hyp = computeExpectedCodonProbForSeqs(seqDict,rf,gcDiff,constraintfile,filename,threshold)
	for codon in genCode.keys():
		if p_hyp[codon]!=0:
			print codon,f_gene[codon],fsyn_gene[codon],f_hyp[codon],fsyn_hyp[codon],p_gene[codon],p_hyp[codon],p_gene[codon]/p_hyp[codon]
		else:
			print codon,f_gene[codon],fsyn_gene[codon],f_hyp[codon],fsyn_hyp[codon],p_gene[codon],p_hyp[codon],'NA'
