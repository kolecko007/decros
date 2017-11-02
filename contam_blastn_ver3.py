import glob
import sys
import os
import re
import pickle
from Bio.Blast import NCBIXML

'''
python contam_blastn_ver3.py total ratio hit_threshold hit_length Own_name blastn_xml_file fasta_file min_hit_len_threshold

total = min kmer frequency
hit_treshold = percent identity
hit_length = len_aln/len_query
ratio - kmer =  ratio of own/hit to be non-contaminat
own_name = name of the organism (e.i.: in file where names are: >gnl|Contaminated_dataset|0_a1 12345667 K; 25... - thr name will be read as Contaminated_dataset)
blastn_xml_file - blast against all data from the sequencing project
fasta_file = reformated fasta file names have to be like this: >gnl|yourname|contig_name etc... contig_name etc is the name of the sequence in assembler ouput
numbers of sequences shoudl grow up byt one
all analyzed dataset have to be reformatted like this before the blast database is made.
'''


def make_coverage_dict(filename):
  infile = open(filename)
  lines = infile.readlines()
  cov_dict = {}
  for line in lines:
    cov_dict[line.split()[0]] = line.split()[1]
 # print cov_dict
  return cov_dict

def load_coverage_dict(filename):
  infile = open(filename)
  cov_dict = pickle.load(infile)
  return cov_dict
  
def analyze_hit(hit):
  print '*******************************%s********************************' % (hit.query)
  mskey = 0    #counts number of missin gocverage values - hopefully will be 0.
  own_kmer = cov_dict[hit.query]   #gets coverage values
  if own_kmer == 0:
    own_kmer = 0.000000000000000000000000000000000000000000000000000001
  #own_kmer = own_kmer.split(';')[1]
  status = 'good'
  for aln in hit.alignments:      #Goes through every hit
    #print aln.title
    hitname = aln.title.split()[1]
    hitname = hitname.split('|')[1]    #builds name
    print '***************', hitname
    if hitname != sys.argv[5] and float(aln.hsps[0].identities)/float(len(aln.hsps[0].query)) >= float(hit_treshold) and (float(len(aln.hsps[0].query.replace('-','')))/float(hit.query_letters) >= float(sys.argv[8]) or len(aln.hsps[0].query.replace('-','')) > sys.argv[4]):
      print hitname
      hit_kmer = cov_dict[aln.title.split()[1]]
      if hit_kmer == 0:
        hit_kmer = 0.00000000000000000000000000000000000000000000000001
        
      #hit_kmer = hit_kmer.split(';')[1]
      zaznam.write('%s %s %s %s\n' % (hit.query,  hitname,  own_kmer, hit_kmer))
      outprdel.write(aln.title)
      #print float(own_kmer)/float(hit_kmer)
      if float(own_kmer)/float(hit_kmer) < float(ratio):
        status = 'bad'
        list_contam.append(hitname)
      elif int(own_kmer) <int(total):
        status = 'bad'
  if int(own_kmer) < int(total):
    status = 'bad'
  print status
  return (status, mskey)


def make_seq_dict(filename):
  infile = open(filename)
  line = infile.read()
  infile.close()
  seqs = line.split('>')
  seqs = seqs[1:]
  d = {}
  for seq in seqs:  
    d[seq.split('\n')[0].strip()] = ''.join(seq.split('\n')[1:])
  
  return d



'''*******************Actual Script*******************'''



cov_dict = load_coverage_dict('RPKM_dit.dat')
total = sys.argv[1]
ratio = sys.argv[2]
hit_treshold = sys.argv[3]
good_names = [sys.argv[5]]
own = sys.argv[5]
PRDEL = open('PRDEL.txt','w')
out = open('%s_clean3.fasta' % (own),'w')
deleted = open('%s_deleted3.fasta' % (own),'w')
short = open('%s_short3.fasta' % (own),'w')
outprdel = open('prdel3.txt','w')
PRDEL2 = open('PRDEL23.fas','w')
zaznam = open('%s_records.txt' % (own),'w')
  
  
print 'making sequence dictionary'
d = make_seq_dict(sys.argv[7])
print 'dictionary done'
  
list_contam = []
  
records = NCBIXML.parse(open(sys.argv[6]))
  
counter = range(10000, 1000100, 10000)

lowkmerc = 0
shortc = 0
contamc = 0
good = 0
unrecogn = 0

c = 0
contamc = 0
for record in records:
  c = c + 1
  status = analyze_hit(record)
  stat = status[0]
  if len(d[record.query]) < 150:
    short.write('>%s\n%s\n' % (record.query, d[record.query]))
  elif stat == 'good':
    out.write('>%s\n%s\n' % (record.query, d[record.query]))
  elif stat == 'bad':
    deleted.write('>%s\n%s\n' % (record.query, d[record.query]))
  elif stat == 'PRDEL':
    PRDEL2.write('>%s\n%s\n' % (record.query, d[record.query]))
  else:
    print 'What the fuck, email Martin!!!'
  
  
print '\n\n\n\n\n\n\n\n\n\n***********************************************************************************************************'
print lowkmerc
print shortc
print contamc
print unrecogn
  
  
out.close()
deleted.close()
short.close()
  
uniq = []
for i in list_contam:
  if i not in uniq:
    uniq.append(i)    
  
d = {}
  
for i in uniq:
  d[i] = 0
  

for i in list_contam:
  d[i] = d[i] + 1
  
print d
zaznam2 = open('%s_recordsdict.txt' % (own),'w')
for i in d:
  zaznam2.write('%s %s\n' % (i,d[i]))
