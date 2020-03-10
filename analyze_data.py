#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
import os, shutil, sys
import threading
import main_script as ms

def get_hosts_list(folder):
	records = list(SeqIO.parse('dataset/' + folder +'/sequences.fasta', 'fasta'))

	hosts = []
	for seq in records:
		host = seq.id.split('_')[0]
		if host not in hosts:
			hosts.append(host)

	return hosts

def get_infectors_list(folder):
	trans_file = 'dataset/' + folder +'/transmission_network.txt'
	f = open(trans_file)
	f.readline()
	infectors = []

	for line in f.readlines():
		parts = line.split('\t')
		if not parts[0] in infectors:
			infectors.append(parts[0])

	f.close()
	return infectors

def get_favites_data_list():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	data_list = []
	for folder in folders:
		idx = folder.rindex('_')
		name = folder[:idx]
		if name not in data_list:
			data_list.append(name)

	# print(data_list, len(data_list))
	data_list.sort()
	return data_list

def get_infector_ratio(folder):
	infector_count = len(get_infectors_list(folder))
	host_count = len(get_hosts_list(folder))
	return host_count/infector_count

def get_favites_infector_ratio_summary():
	for group in get_favites_data_list():
		sum = 0
		for i in range(1,21):
			folder = group + '_' + str(i)
			sum += get_infector_ratio(folder)
		avg = sum/20
		print(group, round(avg, 2))

def get_agv_csv_list(csv_list, start_index):
	eff_len = len(csv_list[0].strip().split(',')) - start_index
	sums = []
	for i in range(eff_len):
		sums.append(0.0)
	for line in csv_list:
		parts = line.strip().split(',')
		for i in range(start_index, len(parts)):
			sums[i - start_index] += float(parts[i])

	for i in range(eff_len):
		sums[i] /= len(csv_list)
		sums[i] = round(sums[i], 3)
	
	return sums

def get_favites_result_summary(csv_file):
	file_name = csv_file.split('/')[-1]
	file_name = 'grouped.' + file_name
	result = open('results/analyze_data/' + file_name, 'w+')

	f = open(csv_file)
	result.write(f.readline())
	lines = f.readlines()

	for group in get_favites_data_list():
		grouplines = [idx for idx in lines if idx.startswith(group)]
		avg = get_agv_csv_list(grouplines, 1)
		result.write('{},{},{},{},{},{},{},{},{},{}\n'.format(group,avg[0],avg[1],avg[2],avg[3],avg[4],avg[5],avg[6],avg[7],avg[8]))

	f.close()

def print_data_summary_():
	data_dir = 'dataset/'
	min_leaves_count = 999999
	max_leaves_count = 0
	total_leaves_count = 0

	min_hosts_count = 999999
	max_hosts_count = 0
	total_hosts_count = 0

	folders = next(os.walk(data_dir))[1]
	# folders = ['SEIR01_sl250_mr025_nv10_1']
	for folder in folders:
		print(folder)
		records = list(SeqIO.parse(data_dir + folder +'/sequences.fasta', 'fasta'))

		leaves = []
		for seq in records:
			leaves.append(seq.id.split('_')[0])

		min_leaves_count = min(min_leaves_count, len(leaves))
		max_leaves_count = max(max_leaves_count, len(leaves))
		total_leaves_count += len(leaves)

		hosts = list(set(leaves))
		min_hosts_count = min(min_hosts_count, len(hosts))
		max_hosts_count = max(max_hosts_count, len(hosts))
		total_hosts_count += len(hosts)

	print('min_leaves_count', min_leaves_count)
	print('max_leaves_count', max_leaves_count)
	print('avg_leaves_count', total_leaves_count/len(folders))

	print('min_hosts_count', min_hosts_count)
	print('max_hosts_count', max_hosts_count)
	print('avg_hosts_count', total_hosts_count/len(folders))

def main():
	# get_favites_infector_ratio_summary()
	# get_favites_result_summary('results/favites_directed_comparison/bootstrap.100.phyloscanner.tnet.new.tnet.bias.th.50.csv')
	

if __name__ == "__main__": main()
