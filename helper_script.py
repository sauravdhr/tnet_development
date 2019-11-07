#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
import os, shutil, sys
import threading
import main_script as ms

def run_script(folder, script):
	cmd = './' + folder +'/'+ script
	print(cmd)
	os.system(cmd)

def run_raxml_scripts_with_threading(folder):
	t=[]
	scripts = next(os.walk(folder))[2]
	for script in scripts:
		t.append(threading.Thread(target=run_script, args=(folder, script)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def create_raxml_scripts_with_bootstrap(bootstrap, script_folder):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	cmd_list = []

	for folder in folders:
		RAxML_folder = os.path.abspath(data_dir + folder + '/RAxML_output')
		fasta_file = os.path.abspath(data_dir + folder + '/sequences.fasta')
		if not os.path.exists(RAxML_folder):
			os.mkdir(RAxML_folder)

		cmd = 'raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n favites -k'.format(fasta_file, RAxML_folder, bootstrap)
		# print(cmd)
		cmd_list.append(cmd)

	print(len(cmd_list))
	if not os.path.exists(script_folder):
		os.mkdir(script_folder)

	for i in range(len(cmd_list)):
		script = open(script_folder + '/' + str(i%60) + '.raxml_cmd.sh', 'a+')
		script.write(cmd_list[i] + '\n')
		script.close()


def multithreadings(input_file, output_dir, bootstrap):
	cmd = 'raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n favites -k'.format(input_file, output_dir, bootstrap)
	# print(cmd)
	os.system(cmd)


def run_raxml_without_scripts_threading(bootstrap):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	t = []

	for folder in folders:
		RAxML_folder = os.path.abspath(data_dir + folder + '/RAxML_output')
		fasta_file = os.path.abspath(data_dir + folder + '/sequences.fasta')
		RAxML_info = RAxML_folder + '/RAxML_info.favites'
		if os.path.exists(RAxML_info):
			continue

		t.append(threading.Thread(target=multithreadings, args=(fasta_file, RAxML_folder, bootstrap)))

	# print('len', len(t))
	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def root_raxml_best_tree():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside :', folder)
		best_tree = data_dir + folder + '/RAxML_output/RAxML_bestTree.favites'
		output_folder = os.path.abspath(data_dir + folder + '/RAxML_output')
		cmd = 'raxmlHPC -f I -m GTRGAMMA -t {} -n bestTree.favites -w {}'.format(best_tree, output_folder)
		os.system(cmd)

def run_tnet_new_besttree_multithreaded(times = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]
	t = []

	for folder in folders:
		print('Inside',folder)
		input_dir = data_dir + folder + '/RAxML_output'
		output_dir = 'outputs/' + folder + '/tnet_best_tree/'
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

		tree_file = input_dir + '/RAxML_rootedTree.bestTree.favites'
		out_file = output_dir + '/bestTree.' + str(times) +'.tnet_new'
		t.append(threading.Thread(target=ms.run_tnet_new_multiple_times, args=(tree_file, out_file, times)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def run_tnet_old_besttree(times = 100):
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		input_dir = data_dir + folder + '/RAxML_output'
		output_dir = 'outputs/' + folder + '/tnet_best_tree/'
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

		tree_file = input_dir + '/RAxML_rootedTree.bestTree.favites'
		out_file = output_dir + '/bestTree.' + str(times) +'.tnet_old'
		ms.run_tnet_old_multiple_times(tree_file, out_file, times)

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



def run_phyloscanner_besttree():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print('Inside',folder)
		input_folder = data_dir + folder + '/RAxML_output'
		output_folder = 'outputs/' + folder + '/phyloscanner_best_tree'
		input_file = input_folder + '/RAxML_rootedTree.bestTree.favites'

		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
			cmd = 'PhyloScanner/phyloscanner_analyse_trees_old.R {} favites -ct -od {} s,0 --overwrite --tipRegex="^(.*)_(.*)$"'.format(input_file, output_folder)
			os.system(cmd)

def main():
	# get_sequences_and_network()
	# rename_and_clean_sequences()
	# create_raxml_scripts_with_bootstrap(100, 'raxml_scripts')
	# run_raxml_scripts_with_threading('raxml_scripts')
	# root_raxml_best_tree()
	# run_tnet_new_besttree_multithreaded(1)
	# run_tnet_old_besttree(1)
	# run_phyloscanner_besttree()
	print_data_summary_()

	



if __name__ == "__main__": main()