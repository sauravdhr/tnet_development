# This scripts compares Tnet and Phyloscanner on CDC detaset.
# Each outbreak is handeled separately in this script.
# We know the source of only 10 out of these outbreaks.

# Library Imports
from Bio import SeqIO
import operator
import shutil, os
import get_edges as ge
import main_script as ms

# Global Variables
known_outbreaks = ['AA', 'AC', 'AI', 'AJ', 'AQ', 'AW', 'BA', 'BB', 'BC', 'BJ']
sources = ['AA45','AC124','AI4','AJ199','AQ89','AW2','BA3','BB45','BC46','BJ28']

def get_true_transmission_edges(outbreak):
	seq_list = list(SeqIO.parse('CDC/'+ outbreak +'/sequences.fasta', 'fasta'))

	true_edges = []
	hosts = []
	for seq in seq_list:
		hosts.append(seq.id.split('_')[0])

	hosts = list(set(hosts))
	source = sources[known_outbreaks.index(outbreak)]
	source = source.replace(outbreak, '')
	hosts.remove(source)

	for host in hosts:
		true_edges.append(source +'->'+ host)

	return true_edges

def run_new_tnet_cdc_multithreaded(times = 100):
	for outbreak in known_outbreaks:
		input_folder = 'CDC/' + outbreak + '/tnet_input'
		output_folder = 'CDC/' + outbreak + '/tnet_new_mod_rand_bootstrap'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
			ms.run_tnet_new_single_folder(input_folder, output_folder, times)

def run_old_tnet_cdc(times = 100):
	for outbreak in known_outbreaks:
		input_folder = 'CDC/' + outbreak + '/tnet_input/'
		output_folder = 'CDC/' + outbreak + '/tnet_old_fixed_bootstrap/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		file_list = next(os.walk(input_folder))[2]
		for file in file_list:
			input_file = input_folder + file
			parts = file.split('.')
			output_file = output_folder + parts[1] + '.tnet_old_fixed'
			# print(input_file, output_file)
			ms.run_tnet_old_multiple_times(input_file, output_file, times)

def run_old_tnet_cdc_single_tree(times = 100):
	for outbreak in known_outbreaks:
		input_file = 'CDC/'+outbreak+'/tnet_input/RAxML_rootedTree.25'
		output_folder = 'CDC/'+outbreak+'/tnet_single_tree/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		
		output_file = output_folder + 'single_tree.' + str(times) + '.tnet_old_fixed'
		if not os.path.exists(output_file):
			# print(input_file, output_file)
			ms.run_tnet_old_multiple_times(input_file, output_file, times)

def run_new_tnet_cdc_single_tree(times = 100):
	for outbreak in known_outbreaks:
		input_file = 'CDC/'+outbreak+'/tnet_input/RAxML_rootedTree.25'
		output_folder = 'CDC/'+outbreak+'/tnet_single_tree/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		
		output_file = output_folder + 'single_tree.' + str(times) + '.tnet_new_min'
		if not os.path.exists(output_file):
			# print(input_file, output_file)
			ms.run_tnet_new_multiple_times(input_file, output_file, times)

def create_cdc_tnet_summary_directed(threshold):
	for outbreak in known_outbreaks:
		print('Inside', outbreak)

		input_folder = 'CDC/' + outbreak + '/tnet_old_bootstrap'
		output_folder = 'CDC/' + outbreak + '/tnet_old_bootstrap_summary_directed'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		edge_dict = {}
		result = open(output_folder + '/tnet_old_bootstrap' + '_th_' + str(threshold) + '_summary.csv', 'w+')
		file_list = next(os.walk(input_folder))[2]

		for file in file_list:
			tnet_file = input_folder + '/' + file
			tnet_edges = ge.get_mul_tnet_edges(tnet_file, threshold)
			for edge in tnet_edges:
				if edge in edge_dict:
					edge_dict[edge] += 1
				else:
					edge_dict[edge] = 1

		edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
		for x, y in edge_dict.items():
			result.write('{},{}\n'.format(x, y))

def create_cdc_tnet_summary_undirected(threshold):
	for outbreak in known_outbreaks:
		print('Inside', outbreak)
		input_folder = 'CDC/' + outbreak + '/tnet_new_bootstrap'
		output_folder = 'CDC/' + outbreak + '/tnet_new_bootstrap_summary_undirected'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		edge_dict = {}
		result = open(output_folder + '/tnet_new_bootstrap'+ '_th_' + str(threshold) + '_summary.csv', 'w+')
		file_list = next(os.walk(input_folder))[2]

		for file in file_list:
			tnet_file = input_folder + '/' + file
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
	count = 0
	total = len(known_outbreaks)
	for outbreak in known_outbreaks:
		check_folder = 'CDC/' + outbreak + '/tnet_single_tree/'
		if os.path.exists(check_folder):
			# shutil.rmtree(check_folder)
			check_file = check_folder + 'single_tree.1.tnet_new_min'
			if os.path.exists(check_file):
				os.remove(check_file)
			file_list = next(os.walk(check_folder))[2]
			count += len(file_list)
			# check_file = check_folder + file_list[0]
			# if os.stat(check_file).st_size == 0:
			# 	print(folder)

	print('Progress:', count, 'out of', total*6)


def main():
	# run_new_tnet_cdc_multithreaded(100)
	# run_new_tnet_cdc_single_tree(100)
	create_cdc_tnet_summary_directed(50)
	# create_cdc_tnet_summary_undirected(40)
	# check_and_clean()
	# get_true_transmission_edges('BJ')
	# run_old_tnet_cdc()
	# run_old_tnet_cdc_single_tree(100)


if __name__ == "__main__": main()
