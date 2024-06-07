import pysam
import array
import time
import os
from collections import ChainMap
from collections import defaultdict
bam_file_path = "../FR-HU_dedup.trim.markdup.bam"

def extract_DI_value(data_list):
	for i, item in enumerate(data_list):
		if item[0] == 'DI':
			return i

def find_all_occurrences(input_str, char):
    positions = []
    position = input_str.find(char)
    while position != -1:
        positions.append(position)
        position = input_str.find(char, position + 1)
    return positions

start = time.time()

## 根据比对位置分出大的Duplication sets
N = 0
Chr = 0
Pos = 0
bamfile = pysam.AlignmentFile(bam_file_path, "rb")
outfile = pysam.AlignmentFile(f"{bam_file_path}.modified", "wb", template = bamfile)
for read in bamfile.fetch():
	if(read.reference_name == Chr and read.reference_start == Pos):	
		new_tags = read.tags + [('DI', N)]
	else:
		Chr = read.reference_name
		Pos = read.reference_start
		N += 1
		new_tags = read.tags + [('DI', N)]
	read.tags = new_tags
	#print(read.tags)
	outfile.write(read)

bamfile.close()
outfile.close()

## 根据4sT位置得到sub duplication sets

os.system(f"samtools index {bam_file_path}.modified")

tag = 0
unique_molecular = 0
modified_bamfile = pysam.AlignmentFile(f"{bam_file_path}.modified","rb")
modified_outfile = pysam.AlignmentFile(f"{bam_file_path}.modified.filterd","wb", template = modified_bamfile)
for read in modified_bamfile.fetch():
	quality_string = ''.join(chr(q + 33) for q in read.query_qualities)
	i = extract_DI_value(read.tags)
	Z_pos = find_all_occurrences(quality_string, 'Z')
	F_pos = find_all_occurrences(quality_string, 'F')
	pos = Z_pos + F_pos
	if(read.tags[i][1] != tag):
		tag = read.tags[i][1]
		Query = []
		Query.append(pos)
		new_tags = read.tags + [('Di', 0)]
		read.tags = new_tags
	elif(read.tags[i][1] == tag):
		L = len(Query)
		if pos not in Query:
			Query.append(pos)
			new_tags = read.tags + [('Di', L)]
		else:
			for i,query in enumerate(Query):
				if(pos == query):
					new_tags = read.tags + [('Di', i)]
		read.tags = new_tags
	print(read.tags)
	modified_outfile.write(read)
end = time.time()
print(end-start)
modified_bamfile.close()
modified_outfile.close()

## 根据标签去重

os.system(f"samtools index {bam_file_path}.modified.filterd")

input_file = pysam.AlignmentFile(f"{bam_file_path}.modified.filterd", "rb")
Dedup_file = pysam.AlignmentFile(f"{bam_file_path}.deduped", "wb", header=input_file.header)

reads_to_keep = []

for read in input_file:
	key = str(read.tags[-2][1]) + ":" + str(read.tags[-1][1])
	print(key)
	break
