#! /usr/bin/env python
# aminoAcidAndGeneticCodes.py
# P.Clote

# This file can be imported into another as a module,
# since it simply performs many assignments.

import string

NUCL   = ['A','C','G','U']

stopCodon = ['UGA', 'UAA', 'UAG']

iupacNuc = {'A':['A'],'C':['C'],'G':['G'],'U':['U'],'R':['A','G'],'Y':['C','U'],'S':['G','C'],
			'W':['A','U'],'K':['G','U'],'M':['A','C'],'B':['C','G','U'],'D':['A','G','U'],
			'H':['A','C','U'],'V':['A','C','G'],'N':['A','C','G','U']}

aa2codon = {'*': ['UGA', 'UAA', 'UAG'], 'A': ['GCG', 'GCA', 'GCC', 'GCU'], 'C': ['UGC', 'UGU'],
			'E': ['GAA', 'GAG'], 'D': ['GAU', 'GAC'], 'G': ['GGU', 'GGG', 'GGA', 'GGC'], 'F': ['UUU', 'UUC'],
			'I': ['AUC', 'AUA', 'AUU'], 'H': ['CAU', 'CAC'], 'K': ['AAA', 'AAG'], 'M': ['AUG'],
			'L': ['CUU', 'CUG', 'CUC', 'CUA', 'UUG', 'UUA'], 'N': ['AAC', 'AAU'], 'Q': ['CAA', 'CAG'],
			'P': ['CCU', 'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGU', 'UCU', 'UCG', 'UCA'], 
			'R': ['CGG', 'AGG', 'CGU', 'AGA', 'CGA', 'CGC'], 'T': ['ACC', 'ACA', 'ACG', 'ACU'], 'W': ['UGG'], 
			'V': ['GUC', 'GUU', 'GUG', 'GUA'], 'Y': ['UAU', 'UAC'], 's': ['UCC']}

codon2aa = {}
codon2aa[ 'UUU' ] = 'F'
codon2aa[ 'UCU' ] = 'S'
codon2aa[ 'UAU' ] = 'Y'
codon2aa[ 'UGU' ] = 'C'
codon2aa[ 'UUC' ] = 'F'
codon2aa[ 'UCC' ] = 'S'
codon2aa[ 'UAC' ] = 'Y'
codon2aa[ 'UGC' ] = 'C'
codon2aa[ 'UUA' ] = 'L'
codon2aa[ 'UCA' ] = 'S'
codon2aa[ 'UAA' ] = '*'
codon2aa[ 'UGA' ] = '*'
codon2aa[ 'UUG' ] = 'L'
codon2aa[ 'UCG' ] = 'S'
codon2aa[ 'UAG' ] = '*'
codon2aa[ 'UGG' ] = 'W'
codon2aa[ 'CUU' ] = 'L'
codon2aa[ 'CCU' ] = 'P'
codon2aa[ 'CAU' ] = 'H'
codon2aa[ 'CGU' ] = 'R'
codon2aa[ 'CUC' ] = 'L'
codon2aa[ 'CCC' ] = 'P'
codon2aa[ 'CAC' ] = 'H'
codon2aa[ 'CGC' ] = 'R'
codon2aa[ 'CUA' ] = 'L'
codon2aa[ 'CCA' ] = 'P'
codon2aa[ 'CAA' ] = 'Q'
codon2aa[ 'CGA' ] = 'R'
codon2aa[ 'CUG' ] = 'L'
codon2aa[ 'CCG' ] = 'P'
codon2aa[ 'CAG' ] = 'Q'
codon2aa[ 'CGG' ] = 'R'
codon2aa[ 'AUU' ] = 'I'
codon2aa[ 'ACU' ] = 'T'
codon2aa[ 'AAU' ] = 'N'
codon2aa[ 'AGU' ] = 'S'
codon2aa[ 'AUC' ] = 'I'
codon2aa[ 'ACC' ] = 'T'
codon2aa[ 'AAC' ] = 'N'
codon2aa[ 'AGC' ] = 'S'
codon2aa[ 'AUA' ] = 'I'
codon2aa[ 'ACA' ] = 'T'
codon2aa[ 'AAA' ] = 'K'
codon2aa[ 'AGA' ] = 'R'
codon2aa[ 'AUG' ] = 'M'
codon2aa[ 'ACG' ] = 'T'
codon2aa[ 'AAG' ] = 'K'
codon2aa[ 'AGG' ] = 'R'
codon2aa[ 'GUU' ] = 'V'
codon2aa[ 'GCU' ] = 'A'
codon2aa[ 'GAU' ] = 'D'
codon2aa[ 'GGU' ] = 'G'
codon2aa[ 'GUC' ] = 'V'
codon2aa[ 'GCC' ] = 'A'
codon2aa[ 'GAC' ] = 'D'
codon2aa[ 'GGC' ] = 'G'
codon2aa[ 'GUA' ] = 'V'
codon2aa[ 'GCA' ] = 'A'
codon2aa[ 'GAA' ] = 'E'
codon2aa[ 'GGA' ] = 'G'
codon2aa[ 'GUG' ] = 'V'
codon2aa[ 'GCG' ] = 'A'
codon2aa[ 'GAG' ] = 'E'
codon2aa[ 'GGG' ] = 'G'

#WARNING: one letters added to regular AA codes: letter 'O' represents any aa or stop
IupacAA = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','O']

genCode = { }
genCode[ 'UUU' ] = 'Phe'
genCode[ 'UCU' ] = 'Ser'
genCode[ 'UAU' ] = 'Tyr'
genCode[ 'UGU' ] = 'Cys'
genCode[ 'UUC' ] = 'Phe'
genCode[ 'UCC' ] = 'Ser'
genCode[ 'UAC' ] = 'Tyr'
genCode[ 'UGC' ] = 'Cys'
genCode[ 'UUA' ] = 'Leu'
genCode[ 'UCA' ] = 'Ser'
genCode[ 'UAA' ] = 'Stop'
genCode[ 'UGA' ] = 'Stop'
genCode[ 'UUG' ] = 'Leu'
genCode[ 'UCG' ] = 'Ser'
genCode[ 'UAG' ] = 'Stop'
genCode[ 'UGG' ] = 'Trp'
genCode[ 'CUU' ] = 'Leu'
genCode[ 'CCU' ] = 'Pro'
genCode[ 'CAU' ] = 'His'
genCode[ 'CGU' ] = 'Arg'
genCode[ 'CUC' ] = 'Leu'
genCode[ 'CCC' ] = 'Pro'
genCode[ 'CAC' ] = 'His'
genCode[ 'CGC' ] = 'Arg'
genCode[ 'CUA' ] = 'Leu'
genCode[ 'CCA' ] = 'Pro'
genCode[ 'CAA' ] = 'Gln'
genCode[ 'CGA' ] = 'Arg'
genCode[ 'CUG' ] = 'Leu'
genCode[ 'CCG' ] = 'Pro'
genCode[ 'CAG' ] = 'Gln'
genCode[ 'CGG' ] = 'Arg'
genCode[ 'AUU' ] = 'Ile'
genCode[ 'ACU' ] = 'Thr'
genCode[ 'AAU' ] = 'Asn'
genCode[ 'AGU' ] = 'Ser'
genCode[ 'AUC' ] = 'Ile'
genCode[ 'ACC' ] = 'Thr'
genCode[ 'AAC' ] = 'Asn'
genCode[ 'AGC' ] = 'Ser'
genCode[ 'AUA' ] = 'Ile'
genCode[ 'ACA' ] = 'Thr'
genCode[ 'AAA' ] = 'Lys'
genCode[ 'AGA' ] = 'Arg'
genCode[ 'AUG' ] = 'Met'
genCode[ 'ACG' ] = 'Thr'
genCode[ 'AAG' ] = 'Lys'
genCode[ 'AGG' ] = 'Arg'
genCode[ 'GUU' ] = 'Val'
genCode[ 'GCU' ] = 'Ala'
genCode[ 'GAU' ] = 'Asp'
genCode[ 'GGU' ] = 'Gly'
genCode[ 'GUC' ] = 'Val'
genCode[ 'GCC' ] = 'Ala'
genCode[ 'GAC' ] = 'Asp'
genCode[ 'GGC' ] = 'Gly'
genCode[ 'GUA' ] = 'Val'
genCode[ 'GCA' ] = 'Ala'
genCode[ 'GAA' ] = 'Glu'
genCode[ 'GGA' ] = 'Gly'
genCode[ 'GUG' ] = 'Val'
genCode[ 'GCG' ] = 'Ala'
genCode[ 'GAG' ] = 'Glu'
genCode[ 'GGG' ] = 'Gly'

aaCode = {}
aaCode[ 'ALA' ] =  'A'
aaCode[ 'ARG' ] =  'R'
aaCode[ 'ASN' ] =  'N' 
aaCode[ 'ASP' ] =  'D' 
aaCode[ 'CYS' ] =  'C' 
aaCode[ 'GLN' ] =  'Q' 
aaCode[ 'GLU' ] =  'E' 
aaCode[ 'GLY' ] =  'G' 
aaCode[ 'HIS' ] =  'H' 
aaCode[ 'ILE' ] =  'I' 
aaCode[ 'LEU' ] =  'L' 
aaCode[ 'LYS' ] =  'K' 
aaCode[ 'MET' ] =  'M' 
aaCode[ 'PHE' ] =  'F' 
aaCode[ 'PRO' ] =  'P' 
aaCode[ 'SER' ] =  'S' 
aaCode[ 'THR' ] =  'T' 
aaCode[ 'TRP' ] =  'W' 
aaCode[ 'TYR' ] =  'Y' 
aaCode[ 'VAL' ] =  'V'
aaCode[ 'STOP' ] = '*' 

singleLetterAAcodes = ['A','R','N','D','C','Q','E','G','H','I',
 'L','K','M','F','P','S','T','W','Y','V']

nuclCodes = ['A','C','G','T','U']


threeLetterAAcodes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]


def complement(s):
  #return the complement of a single nucleotide
  #uncertain data is correctly treated and the same case is returned
  t = s.upper()
  if   t=="A" : t = "T"
  elif t=="C" : t = "G"
  elif t=="G" : t = "C"
  elif t in ["T","U"] :  t = "A"
  elif t=="R" : t = "Y"
  elif t=="Y" : t = "R"
  elif t=="M" : t = "K"
  elif t=="K" : t = "M"
  elif t=="S" : t = "W"
  elif t=="W" : t = "S"
  elif t=="H" : t = "D"
  elif t=="D" : t = "H"
  elif t=="B" : t = "V"
  elif t=="V" : t = "B"
  #if t=="N" then complement is "N"
  if s.isupper() :
    return t
  else:
    return t.lower()

def reverseComplement(s):
  #s is string; this function returns the reverse complement
  L = [] 
  for ch in s: L.append(complement(ch))
  L.reverse()   #note: NOT L = L.reverse()
  return( string.join(L,"") ) 
    #return the string corresponding to L 
    #Since string is immutable, string appending creates new strings
    #Doing this for an entire genome is extremely wasteful of space
    #thus we use the mutable list and when finished convert it to a string

def complementStr(s):
  #s is string; this function returns the complement
  L = [] 
  for ch in s: L.append(complement(ch))
  return( string.join(L,"") ) 


def reverse(s):
  #s is string; this function returns the reverse complement
  L = [] 
  for ch in s: L.append(ch)
  L.reverse()   #note: NOT L = L.reverse()
  return( string.join(L,"") ) 


geneticCode3letter = {
'UUU' : 'Phe',
'UCU' : 'Ser',
'UAU' : 'Tyr',
'UGU' : 'Cys',
'UUC' : 'Phe',
'UCC' : 'Ser',
'UAC' : 'Tyr',
'UGC' : 'Cys',
'UUA' : 'Leu',
'UCA' : 'Ser',
'UAA' : 'Stop',
'UGA' : 'Stop',
'UUG' : 'Leu',
'UCG' : 'Ser',
'UAG' : 'Stop',
'UGG' : 'Trp',
'CUU' : 'Leu',
'CCU' : 'Pro',
'CAU' : 'His',
'CGU' : 'Arg',
'CUC' : 'Leu',
'CCC' : 'Pro',
'CAC' : 'His',
'CGC' : 'Arg',
'CUA' : 'Leu',
'CCA' : 'Pro',
'CAA' : 'Gln',
'CGA' : 'Arg',
'CUG' : 'Leu',
'CCG' : 'Pro',
'CAG' : 'Gln',
'CGG' : 'Arg',
'AUU' : 'Ile',
'ACU' : 'Thr',
'AAU' : 'Asn',
'AGU' : 'Ser',
'AUC' : 'Ile',
'ACC' : 'Thr',
'AAC' : 'Asn',
'AGC' : 'Ser',
'AUA' : 'Ile',
'ACA' : 'Thr',
'AAA' : 'Lys',
'AGA' : 'Arg',
'AUG' : 'Met',
'ACG' : 'Thr',
'AAG' : 'Lys',
'AGG' : 'Arg',
'GUU' : 'Val',
'GCU' : 'Ala',
'GAU' : 'Asp',
'GGU' : 'Gly',
'GUC' : 'Val',
'GCC' : 'Ala',
'GAC' : 'Asp',
'GGC' : 'Gly',
'GUA' : 'Val',
'GCA' : 'Ala',
'GAA' : 'Glu',
'GGA' : 'Gly',
'GUG' : 'Val',
'GCG' : 'Ala',
'GAG' : 'Glu',
'GGG' : 'Gly'
}

aaCodes3to1 = {
'ALA' : 'A', 
'ARG' : 'R', 
'ASN' : 'N', 
'ASP' : 'D', 
'CYS' : 'C', 
'GLN' : 'Q', 
'GLU' : 'E', 
'GLY' : 'G', 
'HIS' : 'H', 
'ILE' : 'I', 
'LEU' : 'L', 
'LYS' : 'K', 
'MET' : 'M', 
'PHE' : 'F', 
'PRO' : 'P', 
'SER' : 'S', 
'THR' : 'T', 
'TRP' : 'W', 
'TYR' : 'Y', 
'VAL' : 'V' 
}


aaCodes1to3 = {
  'A' : 'ALA',
  'R' : 'ARG',
  'N' : 'ASN',
  'D' : 'ASP',
  'C' : 'CYS',
  'Q' : 'GLN',
  'E' : 'GLU',
  'G' : 'GLY',
  'H' : 'HIS',
  'I' : 'ILE',
  'L' : 'LEU',
  'K' : 'LYS',
  'M' : 'MET',
  'F' : 'PHE',
  'P' : 'PRO',
  'S' : 'SER',
  'T' : 'THR',
  'W' : 'TRP',
  'Y' : 'TYR',
  'V' : 'VAL'
}

    

kyteDoolittle3Letter = {
  'ALA' : 1.8,
  'ARG' : -4.5,
  'ASN' : -3.5,
  'ASP' : -3.5,
  'CYS' : 2.5,
  'GLN' : -3.5,
  'GLU' : -3.5,
  'GLY' : -0.4,
  'HIS' : -3.2,
  'ILE' : 4.5,
  'LEU' : 3.8,
  'LYS' : -3.9,
  'MET' : 1.9,
  'PHE' : 2.8,
  'PRO' : -1.6,
  'SER' : -0.8,
  'THR' : -0.7,
  'TRP' : -0.9,
  'TYR' : -1.3,
  'VAL' : 4.2
}

def kyteDoolittle(s):
  return kyteDoolittle3Letter[aaCodes1to3[s]]

