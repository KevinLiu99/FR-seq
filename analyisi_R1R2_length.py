from Bio import SeqIO

def compare_read_lengths(file_path_r1, file_path_r2):
    r1_gt_r2_count = 0
    r1_lt_r2_count = 0
    r1_eq_r2_count = 0

    with open(file_path_r1, "r") as handle_r1, open(file_path_r2, "r") as handle_r2:
        for record_r1, record_r2 in zip(SeqIO.parse(handle_r1, "fastq"), SeqIO.parse(handle_r2, "fastq")):
            read_length_r1 = len(record_r1.seq)
            read_length_r2 = len(record_r2.seq)

            # 判断R1>R2、R1<R2和R1=R2
            if read_length_r1 > read_length_r2:
                r1_gt_r2_count += 1
            elif read_length_r1 < read_length_r2:
                r1_lt_r2_count += 1
            else:
                r1_eq_r2_count += 1

    return r1_gt_r2_count, r1_lt_r2_count, r1_eq_r2_count

# 假设R1和R2序列保存在两个FASTQ文件中，分别是"R1.fastq"和"R2.fastq"
file_path_r1 = "../BIR-2_test_dedup_1.fq.gz.processed"
file_path_r2 = "../BIR-2_test_dedup_2.fq.gz.processed"

# 获取三个集合的read对数
count_r1_gt_r2, count_r1_lt_r2, count_r1_eq_r2 = compare_read_lengths(file_path_r1, file_path_r2)
All = count_r1_gt_r2 + count_r1_lt_r2 + count_r1_eq_r2

# 打印结果
print("R1>R2集合的read:", count_r1_gt_r2/All)
print("R1<R2集合的read:", count_r1_lt_r2/All)
print("R1=R2集合的read:", count_r1_eq_r2/All)

