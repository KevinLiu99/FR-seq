sample = "BIR-A"

from Bio import SeqIO
import os
import pysam
import subprocess
import time
from multiprocessing.dummy import Pool
import sys
import gzip
import numpy as np
import random
import pandas as pd
import csv
np.set_printoptions(suppress=True)

#R1,R2比对程序

def R1_R2_alignment(R1,R2):
	Ori = []
	for seq1,seq2 in zip(R1,R2):
		if seq1 == 'T' and seq2 == 'C':
			Ori.append('C')
		elif seq1 == 'G' and seq2 == 'A':
			Ori.append('G')
		elif seq1 == 'C' and seq2 == 'T':
			Ori.append('T')
		elif seq1 == 'A' and seq2 == 'G':
			Ori.append('A')
		elif seq1 == 'A' and seq2 == 'A':
			Ori.append('A')
		elif seq1 == 'C' and seq2 == 'C':
			Ori.append('C')
		elif seq1 == 'G' and seq2 == 'G':
			Ori.append('G')
		elif seq1 == 'T' and seq2 == 'T':
			Ori.append('T')
		else:
			Ori.append('N')

	return Ori

#得到反向互补序列
def get_complement_reverse(seq):
	complement_table = str.maketrans("ATCGU", "TAGCA")
	return seq.translate(complement_table)[::-1]

#检测其他突变频率

def get_mutation_nums(R1,R2):
	M = 0
	K = 0
	S = 0
	W = 0
	for seq1,seq2 in zip(R1,R2):
		if seq1 == 'A' and seq2 == 'C':
			M += 1
		elif seq1 == 'C' and seq2 == 'A':
			M += 1
		elif seq1 == 'G' and seq2 == 'T':
			K += 1
		elif seq1 == 'T' and seq2 == 'G':
			K += 1
		elif seq1 == 'C' and seq2 == 'G':
			S += 1
		elif seq1 == 'G' and seq2 == 'C':
			S += 1
		elif seq1 == 'A' and seq2 == 'T':
			W += 1
		elif seq1 == 'T' and seq2 == 'A':
			W += 1
	
	
	return M, K, S, W;

#推导original sequence
os.chdir("./")
start = time.time()
Original = {}
for R1,R2 in zip(SeqIO.parse(f'../{sample}_1.fq.gz.processed',"fastq"),SeqIO.parse(f'../{sample}_2.fq.gz.processed',"fastq")):
	if(R1.id == R2.id):
		Ori_id = R1.id
		Original[Ori_id] = R1_R2_alignment(R1.seq,R2.seq)
	else:
		pass

end = time.time()
print("seconds:", end-start)

#检测reads mutation number

start = time.time()
read_num = 0
mutation_num = 0
AC_mut = GT_mut = CG_mut = AT_mut = 0

for R1,R2 in zip(SeqIO.parse(f'../{sample}_1.fq.gz.processed',"fastq"),SeqIO.parse(f'../{sample}_2.fq.gz.processed',"fastq")):
	read_num += 1
	if(R1.id in Original.keys() and R1.id == R2.id):
		AC,GT,CG,AT = get_mutation_nums(R1.seq,R2.seq)
		mutation_num = AC + GT + CG + AT
	Original.setdefault(R1.id, ).append(mutation_num)
	mutation_num = 0

print(read_num)
end = time.time()
print("seconds:", end-start)

#list mutation num in each read

start = time.time()
mutation_list = []
for key,value in Original.items():
	mutation_list.append(value[-1])
temp = np.array(mutation_list)
unique = np.unique(temp)

Mut = dict.fromkeys(unique,0)
for key,value in Original.items():
	if(value[-1] in Mut.keys()):
		Mut[value[-1]] = int(Mut[value[-1]]) + 1

end = time.time()
print("second:", end-start)

#print results

datas = []
header = list(Mut.keys())
datas.append(Mut)

with open(f'../{sample}_mutation.csv','a') as f:
	writer = csv.DictWriter(f,fieldnames=header)
	writer.writeheader()
	writer.writerows(datas)
