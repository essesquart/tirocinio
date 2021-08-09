import sys

CD_HIT_file = open(sys.argv[2])
lines = CD_HIT_file.readlines()
CD_HIT_file.close()
for a in range(0,len(lines)):
	lines[a]=lines[a].split('\n')[0]

fasta_file = open(sys.argv[1])
fasta = fasta_file.readlines()
fasta_file.close()
fasta_dict={}
for a in range(0,len(fasta),2):
	fasta_dict[fasta[a].split('\n')[0]]=fasta[a+1].split('\n')[0]

cluster_counter = 0
for line in lines:
	if line[0] == '>':
		current_output = []
		current_cluster = line.split(' ')[1]
		print(current_cluster)
		cluster_counter = cluster_counter+1
		size = 0
		while lines[cluster_counter][0] != '>':
			size = size+1
			seq = lines[cluster_counter].split(', ')[1].split('...')[0]	
			current_output.append(seq)
			current_output.append(fasta_dict[seq])
			cluster_counter = cluster_counter+1
			if cluster_counter == len(lines):
				break
		if size >= int(sys.argv[3]):
			print("Stampato:"+ current_cluster)
			with open(sys.argv[4]+'/'+ current_cluster+'.fa','w') as fasta_file:
				for x in current_output:
					fasta_file.writelines(x + "\n")

