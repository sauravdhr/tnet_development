#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
# import get_edges as ge
import operator
import os, shutil, sys
import threading

def get_sequences_and_network():
	favites_data_main = '/home/saurav/research/Favites_data_from_sam/'
	datasets = next(os.walk(favites_data_main))[1]
	# print(datasets)

	for dataset in datasets:
		cur_dir = favites_data_main + dataset
		folders = next(os.walk(cur_dir))[1]
		# print(folders)

		for folder in folders:
			fasta_file = cur_dir+ '/' +folder+ '/FAVITES_output/error_free_files/sequence_data.fasta'
			net_file = cur_dir+ '/' +folder+ '/FAVITES_output/error_free_files/transmission_network.txt'
			new_dir = 'dataset/' +folder
			if not os.path.exists(new_dir):
				os.mkdir(new_dir)
			fasta_copy = new_dir+ '/sequence_data.fasta'
			net_copy = new_dir+ '/transmission_network.txt'
			if os.path.exists(fasta_file):
				shutil.copy(fasta_file, fasta_copy)
			if os.path.exists(net_file):
				shutil.copy(net_file, net_copy)
			# print(os.path.exists(fasta_file))

def rename_and_clean_sequences():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	# print(folders)
	for folder in folders:
		print(folder)
		favites_fasta = data_dir + folder +'/sequence_data.fasta'
		records = list(SeqIO.parse(favites_fasta, 'fasta'))
		for record in records:
			parts = record.id.split('|')
			record.id = parts[1] + '_' + parts[0]

		SeqIO.write(records, data_dir + folder +'/sequences.fasta', 'fasta')
		os.remove(favites_fasta)

def multithreadings(input_file, output_dir, bootstrap):
	cmd = 'raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n favites -k'.format(input_file, output_dir, bootstrap)
	# print(cmd)
	os.system(cmd)

def run_raxml_with_threading(bootstrap):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	t = []

	for folder in folders:
		RAxML_folder = os.path.abspath(data_dir + folder + '/RAxML_output')
		fasta_file = os.path.abspath(data_dir + folder + '/sequences.fasta')
		RAxML_info = RAxML_folder + '/RAxML_info.favites'
		if not os.path.exists(RAxML_folder):
			os.mkdir(RAxML_folder)
		if os.path.exists(RAxML_info):
			continue

		t.append(threading.Thread(target=multithreadings, args=(fasta_file, RAxML_folder, bootstrap)))

	# print('len', len(t))
	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def create_bootstrap_trees():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		bootstrap_file = data_dir + folder + '/RAxML_output/RAxML_bootstrap.favites'
		bootstrap_folder = data_dir + folder + '/RAxML_output/bootstrap_trees'
		if not os.path.exists(bootstrap_file):
			continue
		if not os.path.exists(bootstrap_folder):
			os.mkdir(bootstrap_folder)

		f = open(bootstrap_file)
		tree_list = f.readlines()

		for i in range(len(tree_list)):
			file = open(bootstrap_folder + '/' + str(i) + '.bootstrap.tree', 'w')
			file.write(tree_list[i])

		# break

def root_bootstrap_trees():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		bootstrap_folder = data_dir + folder + '/RAxML_output/bootstrap_trees'
		rooted_bootstrap_folder = data_dir + folder + '/rooted_bootstrap_trees'
		bootstrap_trees = next(os.walk(bootstrap_folder))[2]
		output_folder = os.path.abspath(rooted_bootstrap_folder)
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		for tree in bootstrap_trees:
			input_tree = bootstrap_folder + '/' + tree
			i = int(tree.split('.')[0])
			cmd = 'raxmlHPC -f I -m GTRGAMMA -t {} -n {} -w {}'.format(input_tree, str(i), output_folder)
			# print(cmd)
			os.system(cmd)
			try:
				os.remove(output_folder + '/RAxML_info.' + str(i))
			except:
				print('RAxML_info does not exist')
		# break

def create_phyloscanner_input():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		input_folder = data_dir + folder + '/rooted_bootstrap_trees'
		output_folder = data_dir + folder + '/phyloscanner_input'
		tree_list = next(os.walk(input_folder))[2]
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		for tree in tree_list:
			i = int(tree.rstrip().split('.')[1])
			rooted_tree = input_folder + '/' + tree
			rename_tree = output_folder + '/bootstrap.InWindow_'+ str(1000+i*100) +'_to_'+ str(1099+i*100) +'.tree'
			shutil.copy(rooted_tree, rename_tree)

def create_phyloscanner_input_with_tree_size(size):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		input_folder = data_dir + folder + '/phyloscanner_input'
		output_folder = data_dir + folder + '/phyloscanner_input_' + str(size)
		tree_list = next(os.walk(input_folder))[2]
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		for i in range(size):
			rooted_tree = input_folder + '/' + tree_list[i]
			rename_tree = output_folder + '/' + tree_list[i]
			shutil.copy(rooted_tree, rename_tree)

def run_phyloscanner(bootstrap = 0):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		if bootstrap == 0:
			input_folder = data_dir + folder + '/phyloscanner_input'
			output_folder = 'outputs/' + folder + '/phyloscanner_output_100_bootstrap'
			input_file = input_folder + '/bootstrap.InWindow_'
		else:
			input_folder = data_dir + folder + '/phyloscanner_input_' + str(bootstrap)
			output_folder = 'outputs/' + folder + '/phyloscanner_output_' + str(bootstrap) + '_bootstrap'
			input_file = input_folder + '/bootstrap.InWindow_'

		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
			cmd = 'PhyloScanner/phyloscanner_analyse_trees_old.R {} favites -od {} s,0 --overwrite --tipRegex="^(.*)_(.*)$"'.format(input_file, output_folder)
			os.system(cmd)

		# cmd = 'PhyloScanner/phyloscanner_analyse_trees_old.R {} favites -ct -od {} s,0 --overwrite --tipRegex="^(.*)_(.*)$"'.format(input_file, output_folder)
		# print(cmd)
		# os.system(cmd)

		# break

def cmd_phyloscanner(input_file, output_folder):
	cmd = 'PhyloScanner/phyloscanner_analyse_trees_old.R {} favites -od {} s,0 --overwrite --tipRegex="^(.*)_(.*)$"'.format(input_file, output_folder)
	print(cmd)
	# os.system(cmd)

def run_phyloscanner_multithreaded(bootstrap = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	t = []

	for folder in folders:
		if bootstrap == 100:
			input_folder = data_dir + folder + '/phyloscanner_input'
			output_folder = 'outputs/' + folder + '/phyloscanner_output_100_bootstrap'
			input_file = input_folder + '/bootstrap.InWindow_'
		else:
			input_folder = data_dir + folder + '/phyloscanner_input_' + str(bootstrap)
			output_folder = 'outputs/' + folder + '/phyloscanner_output_' + str(bootstrap) + '_bootstrap'
			input_file = input_folder + '/bootstrap.InWindow_'

		if not os.path.exists(output_folder):
			print('Inside',folder)
			os.mkdir(output_folder)
			t.append(threading.Thread(target=cmd_phyloscanner, args=(input_file, output_folder)))

	print('len_t', len(t))
	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def run_tnet_old_multiple_times(input_file, output_file, time = 100):
	temp_out_file = output_file + '.temp'
	edge_dict = {}
	result = open(output_file, 'w+')

	for t in range(time):
		cmd = 'python3 tnet_old/tnet.py {} {}'.format(input_file, temp_out_file)
		# print(cmd)
		os.system(cmd)
		e_list = []
		print('Run', t)

		# Read result from temp_out_file and save to edge_dict
		f = open(temp_out_file)
		f.readline()
		for line in f.readlines():
			parts = line.rstrip().split('\t')
			edge = parts[0]+'->'+parts[1]

			if edge not in e_list:
				if edge in edge_dict:
					edge_dict[edge] += 1
				else:
					edge_dict[edge] = 1
				e_list.append(edge)

		f.close()
		os.remove(temp_out_file)
		os.remove(input_file + '.temp')
		os.remove(input_file + '.tnet.log')

	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	# print(edge_dict)

	for x, y in edge_dict.items():
		result.write('{}\t{}\n'.format(x, y))

	result.close()

def run_tnet_new_multiple_times(input_file, output_file, time = 100):
	temp_out_file = output_file + '.temp'
	edge_dict = {}
	source_count = {}
	result = open(output_file, 'w+')

	for t in range(time):
		cmd = 'python3 tnet_dev.py {} {}'.format(input_file, temp_out_file)
		# print(cmd)
		os.system(cmd)
		e_list = []
		print('Run', t)

		# Read result from temp_out_file and save to edge_dict
		f = open(temp_out_file)
		f.readline()
		for line in f.readlines():
			parts = line.rstrip().split('\t')
			edge = parts[0]+'->'+parts[1]

			if edge not in e_list:
				if edge in edge_dict:
					edge_dict[edge] += 1
				else:
					edge_dict[edge] = 1
				e_list.append(edge)

		f.close()
		os.remove(temp_out_file)

	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	# print(edge_dict)

	for x, y in edge_dict.items():
		result.write('{}\t{}\n'.format(x, y))

	result.close()

def run_tnet_old(times = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		tree_file = data_dir + folder + '/RAxML_rooted_tree.tree'
		out_file = data_dir + folder + '/tnet_old_' + str(times) + '.tnet'
		# print(tree_file, out_file, times)
		run_tnet_old_multiple_times(tree_file, out_file, times)
		# break

def run_tnet_old_multithreaded(times = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	t = []

	for folder in folders:
		print('Inside',folder)
		input_dir = data_dir + folder + '/rooted_bootstrap_trees'
		output_dir = 'outputs/' + folder + '/tnet_old_' + str(times) + '_times'
		if not os.path.exists('outputs/' + folder):
			os.mkdir('outputs/' + folder)
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
		tree_list = next(os.walk(input_dir))[2]

		for tree in tree_list:
			tree_file = input_dir + '/' + tree
			name = tree.split('.')[1]
			out_file = output_dir + '/' + name +'.tnet_old'
			t.append(threading.Thread(target=run_tnet_old_multiple_times, args=(tree_file, out_file, times)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def run_tnet_new_single_folder(input_dir, output_dir, times = 100):
	t = []
	tree_list = next(os.walk(input_dir))[2]
	for tree in tree_list:
		tree_file = input_dir + '/' + tree
		name = tree.split('.')[1]
		out_file = output_dir + '/' + name +'.tnet'
		if not os.path.exists(out_file):
			t.append(threading.Thread(target=run_tnet_new_multiple_times, args=(tree_file, out_file, times)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()


def run_tnet_new_multithreaded(times = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		input_dir = data_dir + folder + '/rooted_bootstrap_trees'
		output_dir = 'outputs/' + folder + '/tnet_new_' + str(times) + '_bootstrap'
		if not os.path.exists('outputs/' + folder):
			os.mkdir('outputs/' + folder)
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
			run_tnet_new_single_folder(input_dir, output_dir, times)

def create_tnet_bootstrap_output(bootstrap):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		phylo_bootstrap_folder = data_dir + folder + '/phyloscanner_input_' + str(bootstrap)
		if not os.path.exists(phylo_bootstrap_folder):
			print('Invalid bootstrap value')
			return

		input_folder = 'outputs/' + folder + '/tnet_new_100_bootstrap/'
		output_folder = 'outputs/' + folder + '/tnet_new_' + str(bootstrap) + '_bootstrap'
		file_list = next(os.walk(phylo_bootstrap_folder))[2]
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
			for file in file_list:
				window = int(file.split('_')[1])
				index = (window - 1000)//100
				source_file = input_folder + str(index) + '.tnet_new'
				copy_file = output_folder + '/' + str(index) + '.tnet_new'
				print(source_file, copy_file)
				shutil.copy(source_file, copy_file)

def create_directed_tnet_bootstrap_summary(tree_folder, threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		edge_dict = {}
		bootstrap_folder = data_dir + folder + '/' + tree_folder
		output_folder = data_dir + folder + '/tnet_new_bootstrap_summary_directed/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		if not os.path.exists(output_folder + tree_folder + '_th_' + str(threshold) + '_summary.csv'):
			result = open(output_folder + tree_folder + '_th_' + str(threshold) + '_summary.csv', 'w+')
			file_list = next(os.walk(bootstrap_folder))[2]

			for file in file_list:
				tnet_file = bootstrap_folder + '/' + file
				tnet_edges = ge.get_mul_tnet_edges(tnet_file, threshold)
				for edge in tnet_edges:
					if edge in edge_dict:
						edge_dict[edge] += 1
					else:
						edge_dict[edge] = 1

			edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
			for x, y in edge_dict.items():
				result.write('{},{}\n'.format(x, y))

def create_undirected_tnet_bootstrap_summary(tree_folder, threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		edge_dict = {}
		bootstrap_folder = data_dir + folder + '/' + tree_folder
		output_folder = data_dir + folder + '/tnet_new_bootstrap_summary_undirected/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		if not os.path.exists(output_folder + tree_folder + '_th_' + str(threshold) + '_summary.csv'):
			result = open(output_folder + tree_folder + '_th_' + str(threshold) + '_summary.csv', 'w+')
			file_list = next(os.walk(bootstrap_folder))[2]

			for file in file_list:
				tnet_file = bootstrap_folder + '/' + file
				tnet_edges = ge.get_mul_tnet_undirected_edges(tnet_file, threshold)
				for edge in tnet_edges:
					parts_edge = edge.rstrip().split('->')
					rev_edge = parts_edge[1]+ '->' +parts_edge[0]
					if edge in edge_dict:
						edge_dict[edge] += 1
					elif rev_edge in edge_dict:
						edge_dict[rev_edge] += 1
					else:
						edge_dict[edge] = 1

			edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
			for x, y in edge_dict.items():
				result.write('{},{}\n'.format(x, y))

def check_and_clean():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	count = 0

	for folder in folders:
		# print('Inside',folder)
		check_folder = 'outputs/' + folder + '/tnet_best_tree'
		if os.path.exists(check_folder):
			# original = '/home/saurav/research/FAVITES_compare_TNet_v2/outputs/'+ folder +'/tnet_best_tree/bestTree.1.tnet_new'
			# copy = 'outputs/'+ folder +'/tnet_best_tree/bestTree.1.tnet_new'
			# if os.path.exists(original):
			# 	shutil.copy(original, copy)
			# count += 1
			file_list = next(os.walk(check_folder))[2]
			count += len(file_list)
			if len(file_list) == 10: print(folder)
			for file in file_list:
				if "tnet_new_rand_mod_bug_fixed" in file:
					check_file = check_folder + '/' + file
					print(check_file)
					os.remove(check_file)
			# 	print(folder)
			# 	# os.remove(check_file)
			# 	shutil.rmtree(check_folder + '/tnet_new_10_bootstrap')

	print('Done',count, 'out of:', len(folders)*9)

def main():
	# get_sequences_and_network()
	# rename_and_clean_sequences()
	# run_raxml_with_threading(100)
	# create_bootstrap_trees()
	# root_bootstrap_trees()
	# create_phyloscanner_input()
	# create_phyloscanner_input_with_tree_size(50)
	# run_phyloscanner(50)
	# run_phyloscanner_multithreaded(50)
	# run_tnet_old_multithreaded()
	# run_tnet_new_multithreaded()
	# create_tnet_bootstrap_output(10)
	# create_tnet_bootstrap_output(50)
	# create_directed_tnet_bootstrap_summary('tnet_new_10_bootstrap', 30)
	# create_undirected_tnet_bootstrap_summary('tnet_new_10_bootstrap', 30)
	check_and_clean()




if __name__ == "__main__": main()