from Bio import SeqIO
import os
import pysam
from itertools import islice
import subprocess
import time
from multiprocessing.dummy import Pool
import sys
import gzip
import numpy as np
import random

sample = "M4sT_dedup.trim.CFS"

List_4sT = {}

with open(f"../{sample}.sam","r") as samfile:
	lines = samfile.readlines()
	for r_lines in islice(lines,98,None):
		r_lines = r_lines.split("\t")
		Name = r_lines[0]
		Call_4sT = r_lines[10]
		List_4sT[Name] = Call_4sT



N = 0
count = 0

for key,value in List_4sT.items():
	if "Z" in value or "F" in value:
		Z_num = value.count("Z")
		F_num = value.count("F")
		N = N + Z_num + F_num
		count += 1
print("4st number after processed:" + str(N))
print("reads include at least one 4st:" + str(count))

				

table_4sT = [[0 for _ in range(10)] for _ in range(10)]
n = 0
Table = []
BIR = 0
BIR_read = 0
S_replication = 0
S_read = 0
for key,value in List_4sT.items():
	n += 1
	Z_count = value.count("Z")
	F_count = value.count("F")
	if Z_count == 0 and F_count == 0:
		continue;
	elif Z_count >=1 and F_count >= 1:
		BIR += Z_count + F_count
		BIR_read += 1
	else:
		S_replication += Z_count + F_count
		S_read += 1
	Z_count = min(9,max(0,Z_count))
	F_count = min(9,max(0,F_count))
	table_4sT[F_count][Z_count] += 1
for i in range(10):
	for j in range(10):
		Table.append(table_4sT[i][j])
print(Table)
print(n)
print("4sT in BIR_region:" + " " + str(BIR))
print("4sT in cReplication:" + " " + str(S_replication))
print("4sT read in BIR:" + " " + str(BIR_read))
print("4sT read in cReplication:" + " " + str(S_read))


