#! /usr/bin/env python
import os,sys
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
 	
def GCcont(seq):
	sum=0
	for ch in seq:
		if ch=='G' or ch=='C':
			sum += 1
	return sum


if __name__=='__main__':
	if len(sys.argv)!=2:
		print "Usage: %s fastafile" % sys.argv[0]
		sys.exit(1)
	seqdict = parseFasta(sys.argv[1])
	gclist=[]
	for acc,seq in seqdict.iteritems():
		gclist.append(GCcont(seq)/float(len(seq)) * 100)
	print "GC:%.1f - %.1f"%(min(gclist),max(gclist))
		
		
