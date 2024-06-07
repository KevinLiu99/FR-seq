sample = "FR-HU_dedup"

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

def get_mutation_rates(R1,R2):
	M = 0
	K = 0
	S = 0
	W = 0
	num = 0
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
		num += 1
	
	M_mut = M/(num*2)
	K_mut = K/(num*2)
	S_mut = S/(num*2)
	W_mut = W/(num*2)
	
	return M_mut, K_mut, S_mut, W_mut;

# get_information函数

def get_information(r_lines):
	ref = {}
	r_lines = r_lines.split('\t')
	flag = r_lines[1]
	if (flag == '4'):
		sys.exit()
	id = r_lines[0]
	Chr = r_lines[2]
	Start = r_lines[3]
	End = str(int(r_lines[3]) + len(r_lines[9]))
	pos = Chr + ":" + Start + "-" + End
	ref_seq = os.popen('samtools faidx /lustre/user/liclab/lisky/lkw/proj/genome/Gallus_gallus/galGal6/galGal6.fa %s | awk "NR != 1{print}" | xargs echo -n | sed s/[[:space:]]//g' % pos)
	ref_seq = ref_seq.read()
	ref_seq = ref_seq.upper()
	if (flag == '16'):
		ref_seq = get_complement_reverse(ref_seq)
	key = id
	ref[key] = [flag,pos]
	ref.setdefault(key, ).append(ref_seq)
	return ref

#Step1:数据质控

os.system(f"/lustre/user/liclab/lisky/lkw/proj/HiC_met/meth_tools/hammer-seq-master/hbs_tools.v1.1/hbs_process path/to/R1 path/to/R2 path/to/hairpin_adapter path/to/sequencing_adapter")

#Step7:去重

os.system(f"fastp -i R1_input -I R2_input -o R1_output -O R2_output -D")

#Step2:还原original sequence
os.chdir("./")
start = time.time()
Original = {}
for R1,R2 in zip(SeqIO.parse(f"../{sample}_1.fq.gz.processed","fastq"),SeqIO.parse(f"../{sample}_2.fq.gz.processed","fastq")):
	if(R1.id == R2.id):
		Ori_id = R1.id
		Original[Ori_id] = R1_R2_alignment(R1.seq,R2.seq)
	else:
		pass

end = time.time()
print("seconds:", end-start)
'''
#4ST calling测试

ST = []
for Seq1,Seq2,Seq3 in zip(test_11, test_12, Ori_2):
	if (Seq1 == 'A' and Seq2 == 'G' and Seq3 == 'A'):
		ST.append('F')
	elif (Seq1 == 'C' and Seq2 == 'T' and Seq3 == 'T'):
		ST.append('Z')
	else:
		ST.append('-')

print (ST)
'''
#Step3:去除错误建库结构及低质量hairpin片段

for k in list(Original.keys()):
	S = Original[k]
	num =0
	for sequence in S:
		if (sequence == 'N'):
			num = num + 1
		else:
			num = num
	#print (num)
	if (num >= 5):
		del Original[k]

#突变率检测

start = time.time()
read_num = 0
AC_mut = GT_mut = CG_mut = AT_mut = 0

for R1,R2 in zip(SeqIO.parse(f"../{sample}_1.fq.gz.processed","fastq"),SeqIO.parse(f"../{sample}_2.fq.gz.processed","fastq")):
	read_num += 1
	if(R1.id in Original.keys() and R1.id == R2.id):
		AC,GT,CG,AT = get_mutation_rates(R1.seq,R2.seq)
	else:
		continue
	AC_mut = AC_mut + AC
	GT_mut = GT_mut + GT
	CG_mut = CG_mut + CG
	AT_mut = AT_mut + AT

AC_mut_rate = AC_mut/read_num
GT_mut_rate = GT_mut/read_num
CG_mut_rate = CG_mut/read_num
AT_mut_rate = AT_mut/read_num
print (AC_mut_rate)
print (GT_mut_rate)
print (CG_mut_rate)
print (AT_mut_rate)

end = time.time()
print("seconds:", end-start)


#Step2:4ST calling
number = 0
N = 0
start = time.time()
for R1,R2 in zip(SeqIO.parse(f"../{sample}_1.fq.gz.processed","fastq"), SeqIO.parse(f"../{sample}_2.fq.gz.processed","fastq")):
	if ( R1.id == R2.id and R1.id in Original.keys()):
		call_4ST = []
		key = R1.id
		for Seq1,Seq2,Qual1,Qual2 in zip(R1.seq, R2.seq,R1.letter_annotations["phred_quality"],R2.letter_annotations["phred_quality"]):
			if (Seq1 == 'A' and Seq2 == 'G'):
               			call_4ST.append('F')
			elif (Seq1 == 'C' and Seq2 == 'T'):
				call_4ST.append('Z')
			elif ((Seq1 == 'G' and Seq2 == 'A') or (Seq1 == 'T' and Seq2 == 'C')):
				call_4ST.append('C')
			else:
				call_4ST.append('-')
		Original.setdefault(key, ).append(call_4ST)
		if 'Z' in call_4ST or 'F' in call_4ST:
			number += 1
			N = call_4ST.count('Z') + call_4ST.count('F') + N
	else:
		continue

print("4st number original:" + str(N))
print("reads include at least one 4st original:" + str(number))
end = time.time()
print("seconds:", end-start)


#Step4:将Original序列含有4ST的读为fastq文件
start = time.time()
with open('../original.fastq','w') as Origi:
	for key,value in Original.items():
		ID = "@" + key
		Value = value[-1]
		Seq = str(''.join(value[0:-1]))
               # if ("Z" in Value or "F" in Value):
		Origi.write(str(ID) + "\n" + Seq + "\n" + "+" + "\n" + str(''.join(Value)) + "\n")
               # else:
               #         continue

#for key in list(Original.keys()):
#	if ("Z" in Original[key][-1] or "F" in Original[key][-1]):
#		continue
#	else:
#		del Original[key]

end = time.time()
print("seconds:", end-start)
'''
# 将Original序列不含有4ST的读为另一个fastq文件

with open('../Original_None.fastq','w') as zero:
        for key,value in Original.items():
                ID = "@" + key
                Value = value[-1]
                Seq = str(''.join(value[0:-1]))
                if ("Z" in Value or "F" in Value):
                        continue
                else:
                        zero.write(str(ID) + "\n" + Seq + "\n" + "+" + "\n" + str(''.join(Value)) + "\n")

'''
#Step5: Genome Alignment
start = time.time()
os.system(f"bowtie2 -p 8 -x /lustre/user/liclab/lisky/lkw/proj/genome/hg19index_bt2/hg19 -U ../original.fastq --al ../original_align.fastq --un ../original_unalign.fastq -S ../tmp1.sam")
os.system(f"bowtie2 -p 8 -x /lustre/user/liclab/lisky/lkw/proj/genome/Gallus_gallus/galGal6/galGal6.fa -U ../original_align.fastq --al ../original_align_both.fastq --un ../original_human.fastq -S ../tmp2.sam")
os.system(f"bowtie2 -p 8 -x /lustre/user/liclab/lisky/lkw/proj/genome/Gallus_gallus/galGal6/galGal6.fa -U ../original_unalign.fastq --al ../original_gallus.fastq --un ../original_unalign_both.fastq -S ../tmp3.sam")
os.system(f"rm ../tmp*")
os.system(f"bowtie2 -p 8 -x /lustre/user/liclab/lisky/lkw/proj/genome/hg19index_bt2/hg19 -U ../original_human.fastq -S ../{sample}.sam")
os.system(f"bowtie2 -p 8 -x /lustre/user/liclab/lisky/lkw/proj/genome/Gallus_gallus/galGal6/galGal6.fa -U ../original_gallus.fastq -S ../{sample}_gallus.sam")
os.system(f"samtools view -@ 8 -q 20 -h ../{sample}.sam > ../{sample}.trim.sam")
os.system(f"samtools view -hbS -@ 16 ../{sample}.trim.sam > ../{sample}.trim.bam")
#os.system(f"samtools view -hSf 16 -@ 16 ../{sample}.trim.sam > ../{sample}.rev.trim.sam")
#os.system(f"samtools view -hSF 16 -@ 16 ../{sample}.trim.sam > ../{sample}.for.trim.sam")
end = time.time()
print("seconds:", end-start)

## extract 4sT pattern from samfile

List_4sT = {}

with open(f"../{sample}.trim.sam","r") as samfile:
	lines = samfile.readlines()
	for r_lines in islice(lines,96,None):
		r_lines = r_lines.split("\t")
		Name = r_lines[0]
		Call_4sT = r_lines[10]
		flag = r_lines[1]
		if flag == "16":
			Call_4sT = Call_4sT[::-1]
		List_4sT[Name] = Call_4sT


## count 4sT & 脱氨C 
N = 0
Z = 0
F = 0
C = 0
count = 0

for key,value in List_4sT.items():
	if "Z" in value or "F" in value:
		Z_num = value.count("Z")
		F_num = value.count("F")
		C_num = value.count("C")
		N = N + Z_num + F_num
		Z = Z + Z_num
		F = F + F_num
		C = C + C_num
		count += 1
print("R1 4st number after processed:" + str(Z))
print("R2 4st number after processed:" + str(F))
print("4st number after processed:" + str(N))
print("reads include at least one 4st:" + str(count))
print("脱氨的C的个数为：" + str(C))

##Step6: 根据4sT分布绘制表格(人) 			

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
		pass;
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

'''
#选出正负链均含有4st的reads（个数>=2）
start = time.time()

R = 0
LIST = []
for key,value in List_4sT.items():
	Z_count = value.count('Z')
	F_count = value.count('F')
	if Z_count >= 2 and F_count >= 2:
		R = R + 1			
		LIST.append(key)
	else:
		continue

end = time.time()
print("seconds:", end-start)

#打印reads清单

with open("Reads_list.txt","w") as name:
	for i in LIST:
		name.write(i + "\n")

#从sam文件里提取目标reads

os.system(f"samtools sort -@16 -o  ../{sample}.trim.bam.sorted  ../{sample}.trim.bam")

os.system(f"samtools index ../{sample}.trim.bam.sorted")

target_read_ids_set = set(LIST)

samfile = pysam.AlignmentFile(f"../{sample}.trim.bam.sorted", "rb")

outfile = pysam.AlignmentFile(f"../{sample}.valid.bam", "wb", header=samfile.header)

for read in samfile.fetch():
    if read.query_name in target_read_ids_set:
        outfile.write(read)

samfile.close()
outfile.close()
'''
##Test1:4sT在All, Broken fork & Reversal fork上的分布

# 4st在reads上位置_all
start = time.time()
Z1 = [0]*150
F1 = [0]*150
for key in list(List_4sT.keys()):
	for i in range(0,len(List_4sT[key])):
		
		P = List_4sT[key][i]
		
		if P == "Z":
			Z1[i] = Z1[i] + 1
		elif P == "F":
			F1[i] = F1[i] + 1
		else:
			continue 

end = time.time()
print("seconds:", end-start)

# 4st在reads上位置_one4sT
start = time.time()
Z2 = [0]*150
F2 = [0]*150
N1 = 0
for key in list(List_4sT.keys()):
	if ("Z" in List_4sT[key] and "F" not in List_4sT[key]) or ("F" in List_4sT[key] and "Z" not in List_4sT[key]):
		N1 += 1
		for i in range(0,len(List_4sT[key])):
		
			P = List_4sT[key][i]
		
			if P == "Z":
				Z2[i] = Z2[i] + 1
			elif P == "F":
				F2[i] = F2[i] + 1
			else:
				continue 

end = time.time()
print("seconds:", end-start)
print(N1)
# 4st在reads上位置_two4sT
start = time.time()
Z3 = [0]*150
F3 = [0]*150
N2 = 0
for key in list(List_4sT.keys()):
	if "Z" in List_4sT[key] and "F" in List_4sT[key]:
		N2 += 1
		for i in range(0,len(List_4sT[key])):
		
			P = List_4sT[key][i]
		
			if P == "Z":
				Z3[i] = Z3[i] + 1
			elif P == "F":
				F3[i] = F3[i] + 1
			else:
				continue 

end = time.time()
print("seconds:", end-start)
print(N2)
# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z1),len(F1))+1)]

A = [Position,Z1,F1]

output = open(f"{sample}_position_all.xls","w",encoding='gbk')
for i in range(len(A)):
	for j in range(len(A[i])):
		output.write(str(A[i][j]))
		output.write('\t')
	output.write('\n')
output.close()

# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z2),len(F2))+1)]

B = [Position,Z2,F2]

output = open(f"{sample}_position_one4sT.xls","w",encoding='gbk')
for i in range(len(B)):
	for j in range(len(B[i])):
		output.write(str(B[i][j]))
		output.write('\t')
	output.write('\n')
output.close()

# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z3),len(F3))+1)]

C = [Position,Z3,F3]

output = open(f"{sample}_position_two4sT.xls","w",encoding='gbk')
for i in range(len(C)):
	for j in range(len(C[i])):
		output.write(str(C[i][j]))
		output.write('\t')
	output.write('\n')
output.close()

#Step6:根据4sT分布绘制表格(鸡)
##步骤同人

List_4sT_gallus = {}

with open(f"../{sample}_gallus.sam","r") as samfile:
	lines = samfile.readlines()
	for r_lines in islice(lines,466,None):
		r_lines = r_lines.split("\t")
		Name = r_lines[0]
		Call_4sT = r_lines[10]
		flag = r_lines[1]
		if flag == "16":
			Call_4sT = Call_4sT[::-1]
		List_4sT_gallus[Name] = Call_4sT



N = 0
count = 0

for key,value in List_4sT_gallus.items():
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
for key,value in List_4sT_gallus.items():
	n += 1
	Z_count = value.count("Z")
	F_count = value.count("F")
	if Z_count == 0 and F_count == 0:
		pass;
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

##Test1(鸡)

# 4st在reads上位置_all
start = time.time()
Z1 = [0]*150
F1 = [0]*150
for key in list(List_4sT_gallus.keys()):
	for i in range(0,len(List_4sT_gallus[key])):
		
		P = List_4sT_gallus[key][i]
		
		if P == "Z":
			Z1[i] = Z1[i] + 1
		elif P == "F":
			F1[i] = F1[i] + 1
		else:
			continue 

end = time.time()
print("seconds:", end-start)

# 4st在reads上位置_one4sT
start = time.time()
Z2 = [0]*150
F2 = [0]*150
N1 = 0
for key in list(List_4sT_gallus.keys()):
	if ("Z" in List_4sT_gallus[key] and "F" not in List_4sT_gallus[key]) or ("F" in List_4sT_gallus[key] and "Z" not in List_4sT_gallus[key]):
		N1 += 1
		for i in range(0,len(List_4sT_gallus[key])):
		
			P = List_4sT_gallus[key][i]
		
			if P == "Z":
				Z2[i] = Z2[i] + 1
			elif P == "F":
				F2[i] = F2[i] + 1
			else:
				continue 

end = time.time()
print("seconds:", end-start)
print(N1)
# 4st在reads上位置_two4sT
start = time.time()
Z3 = [0]*150
F3 = [0]*150
N2 = 0
for key in list(List_4sT_gallus.keys()):
	if "Z" in List_4sT_gallus[key] and "F" in List_4sT_gallus[key]:
		N2 += 1
		for i in range(0,len(List_4sT_gallus[key])):
		
			P = List_4sT_gallus[key][i]
		
			if P == "Z":
				Z3[i] = Z3[i] + 1
			elif P == "F":
				F3[i] = F3[i] + 1
			else:
				continue 

end = time.time()
print("seconds:", end-start)
print(N2)
# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z1),len(F1))+1)]

A = [Position,Z1,F1]

output = open(f"{sample}_position_all_gallus.xls","w",encoding='gbk')
for i in range(len(A)):
	for j in range(len(A[i])):
		output.write(str(A[i][j]))
		output.write('\t')
	output.write('\n')
output.close()

# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z2),len(F2))+1)]

B = [Position,Z2,F2]

output = open(f"{sample}_position_one4sT_gallus.xls","w",encoding='gbk')
for i in range(len(B)):
	for j in range(len(B[i])):
		output.write(str(B[i][j]))
		output.write('\t')
	output.write('\n')
output.close()

# 将4st位置分布表格输出

Position = [y for y in range(1,max(len(Z3),len(F3))+1)]

C = [Position,Z3,F3]

output = open(f"{sample}_position_two4sT_gallus.xls","w",encoding='gbk')
for i in range(len(C)):
	for j in range(len(C[i])):
		output.write(str(C[i][j]))
		output.write('\t')
	output.write('\n')
output.close()


