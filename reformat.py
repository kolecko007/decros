import sys
import pickle

infile = open(sys.argv[1])
line = infile.read()
infile.close()

seqs = line.split('>')[1:]

out = open(sys.argv[2],'w')


d = {}

c = 0
for seq in seqs:
  if len(''.join(seq.split('\n')[1:])) >= int(sys.argv[4]):
    out.write('>gnl|%s|%s_%s\n%s\n' % (sys.argv[3], c, seq.split()[0], ''.join(seq.split('\n')[1:])))
    d[seq.split()[0]] = 'gnl|%s|%s_%s' % (sys.argv[3], c, seq.split()[0])
    c = c + 1
out.close()

out = open('%s_names.dict' % (sys.argv[2].split('.')[0]),'wb')

pickle.dump(d, out)
out.close()
