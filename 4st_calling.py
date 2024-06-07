#4ST calling
number = 0
N = 0
n = 0
with open(f"../test_1.fa", "r") as R1, open(f"../test_2.fa", "r") as R2:
	while True:
		R1_id = R1.readline().strip()
		R2_id = R2.readline().strip()
		if not R1_id or not R2_id:
			break
		elif R1_id == R2_id:
			n += 1
			R1_seq = R1.readline().strip()
			R2_seq = R2.readline().strip()
			R1_seq = list(R1_seq)
			R2_seq = list(R2_seq)
			call_4ST = []
			key = R1_id
			for Seq1,Seq2 in zip(R1_seq, R2_seq):
				if (Seq1 == 'A' and Seq2 == 'G'):
               				call_4ST.append('F')
				elif (Seq1 == 'C' and Seq2 == 'T'):
					call_4ST.append('Z')
				else:
					call_4ST.append('-')
			if 'Z' in call_4ST or 'F' in call_4ST:
				number += 1
				N = call_4ST.count('Z') + call_4ST.count('F') + N
			else:
				continue
		else:
			continue
print(n)
print("4st number original:" + str(N))
print("reads include at least one 4st original:" + str(number))

		
