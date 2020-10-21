from sys import argv
from collections import Counter
from collections import defaultdict
dic_pos = defaultdict(list)
#base = ["A","T","C","G"]
with open(argv[1]) as f:
	for line in f:
		line = line.strip().split('\t')
		pos = str(line[2]) + '_' + str(line[6]) + '_' + str(line[7])
		mu = str(line[4])
		dic_pos[pos].append(mu)
f.close()
dic_total_pos = defaultdict(int)
with open(argv[2]) as file2:
	for row in file2:
		row = row.strip().split('\t')
		total_pos = str(row[2]) + '_' + str(row[6]) + '_' + str(row[7])
		dic_total_pos[total_pos] +=1
file2.close()
for k, v in dic_pos.items():
	r = Counter(v)
	con = k.split('_')
	cov = dic_total_pos[k]
	print(con[0],con[1],con[2],cov,r['A'],r['T'],r['C'],r['G'])
