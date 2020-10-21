from sys import argv
num_eng = {}
with open(argv[1]) as f1:
	for line in f1:
		if not line.startswith('#'):
			line = line.strip()
			tmp = line.split('\t')
			k = str(tmp[0])+"_"+str(tmp[1])
		#	print(k)
			num_eng[k] = line

with open(argv[2]) as f2:
	for line in f2:
		if not line.startswith('#'):
			line = line.strip()
			tmp = line.split('\t')
			k = str(tmp[0])+"_"+str(tmp[1])
#			print(k,num_eng.keys())
			if k in num_eng.keys():
				print(num_eng[k])
