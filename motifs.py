import re
from Bio import SeqIO
fin = open('full.fasta', 'r')
fout = open('full_cut.csv', 'w')
line=''
for record in SeqIO.parse(fin,'fasta'):
    with open('Cutmotifs.txt') as motifs:
        for line in motifs:
            line=line.strip()
            liner = record.id+'\t'+str(record.seq)+'\t'+str(len(record.seq))+'\t'+str(line)
            fout.write(liner)
            for m in re.finditer(line,str(record.seq)):
                diner = '\t'+str(m.start())+'\t'+ str(m.end())+'\t'+str(m.group())
                fout.write(diner)
            fout.write('\n')
fin.close()
fout.close()
