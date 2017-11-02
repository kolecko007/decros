import pickle
import glob
import sys
import os


final_d = {}

files = [fname for fname in glob.glob('*_names.dict')]

for f in files:
  d = pickle.load(open(f))
  infile = open('%s_stats.txt' % (f.split('_')[0][2:]))
  lines = infile.readlines()
  infile.close()

  total_mapped = 0
  for line in lines:
    total_mapped = total_mapped + int(line.split()[2])

  for line in lines:
    try:
      RPKM = float(line.split()[2])/((float(line.split()[1])/1000)*(float(total_mapped)/1000000))
    except ZeroDivisionError:
      RPKM = 0
  
    if line.split()[0] != '*':
      final_d[d[line.split()[0]]] = RPKM

print final_d



pickle.dump(final_d, open('RPKM_dit.dat','w'))






