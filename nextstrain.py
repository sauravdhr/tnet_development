#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
from collections import Counter
import csv
import get_edges as ge
import json
import main_script as ms
import operator
import random
import os, shutil, sys
import threading
import xlrd
import tnet_treetime as tnet

def analize_nextstrain_metadata(metadata):
	tsv = csv.reader(open(metadata), delimiter="\t")
	country = []
	for row in tsv:
		if row[8] == 'Human':
			country.append(row[3])

	counter = Counter(country)
	print('Total country', len(counter), 'Seq', sum(counter.values()))
	country.clear()

	for name, count in counter.items():
		if count >= 14:
			country.append(name)

	country.remove('Timor-Leste')
	print('Total country', len(country))

	filtered_sequences = []
	tsv = csv.reader(open(metadata), delimiter="\t")
	for row in tsv:
		if row[8] == 'Human' and row[3] in country:
			filtered_sequences.append(row[7])

	print('filtered_sequences', len(filtered_sequences))
	source_fasta = 'covid_19/GISAID/gisaid_06_12_clean.fasta'
	output_fasta = 'covid_19/nextstrain/nextstrain_sequences_06_12.fasta'
	records = list(SeqIO.parse(source_fasta, 'fasta'))
	new_records = []
	print('records', len(records))
	for record in records:
		if record.id in filtered_sequences:
			new_records.append(record)
		# else:
		# 	print(record.id)

	print('new_records', len(new_records))
	SeqIO.write(new_records, output_fasta, 'fasta')

def main():
	analize_nextstrain_metadata('covid_19/nextstrain/nextstrain_metadata_06_12.tsv')

if __name__ == "__main__": main()
