#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
import get_edges as ge
import operator
import os, shutil, sys
import threading
import xlrd

def create_clean_sequences_gisaid(input_fasta, output_fasta):
	data_dir = 'covid_19/GISAID/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	wb = xlrd.open_workbook(data_dir + 'gisaid_cov2020_acknowledgement_table.xls')
	sheet = wb.sheet_by_index(0)

	id_location_dict = {}
	for i in range(4, sheet.nrows):
		id_location_dict[sheet.cell_value(i, 0)] = sheet.cell_value(i, 2)

	list_temp = []
	for record in records:
		accession_id = record.id.split('|')[1]
		# parts = id_location_dict[accession_id].split('/')
		# print(accession_id, parts)
		# if len(parts) < 3:
		# 	print('Unknown')
		# else:
		# 	print(parts[2])
		list_temp.append(len(record.seq))
		record.id = accession_id
		record.name = ''
		record.description = ''

	print(min(list_temp), max(list_temp))

	# SeqIO.write(records, data_dir + output_fasta, 'fasta')

def create_clean_sequences_ncbi(input_fasta, output_fasta):
	data_dir = 'covid_19/NCBI/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	# f = open(data_dir + 'data_table.csv')
	# f.readline()
	# id_location_dict = {}
	# for line in f.readlines():
	# 	parts = line.split(',')
	# 	print(parts[2])
	# 	id_location_dict[sheet.cell_value(i, 0)] = sheet.cell_value(i, 2)

	temp = []
	for record in records:
		accession_id = record.id.split('|')[1].split('.')[0]
		# parts = id_location_dict[accession_id].split('/')
		# print(accession_id, parts)
		# if len(parts) < 3:
		# 	print('Unknown')
		# else:
		# 	print(parts[2])
		# temp.append(len(record.seq))
		print(accession_id)
		record.id = accession_id
		record.name = ''
		record.description = ''

	# print(len(temp), min(temp), max(temp))

	SeqIO.write(records, data_dir + output_fasta, 'fasta')

def run_raxml_with_pthreads(fasta_file, bootstrap, threads):
	data_dir = 'covid_19/NCBI/'
	RAxML_folder = os.path.abspath(data_dir + '/RAxML_output_complete')
	input_file = os.path.abspath(data_dir + fasta_file)

	if not os.path.exists(RAxML_folder):
		os.mkdir(RAxML_folder)

	cmd = 'raxmlHPC-PTHREADS -T {} -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n ncbi -k'.format(threads, input_file, RAxML_folder, bootstrap)
	# print(cmd)
	os.system(cmd)

def main():
	# create_clean_sequences_gisaid('gisaid_cov2020_sequences_world_complete_high_coverage.fasta', 'clean_sequences_world_nations.fasta')
	# create_clean_sequences_ncbi('ncbi_sars-cov-2_complete_sequences_align.fasta', 'clean_complete_align_sequences.fasta')
	run_raxml_with_pthreads('clean_complete_align_sequences.fasta', 100, 50)


if __name__ == "__main__": main()