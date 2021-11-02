import os
from csv import writer

from utils.utils import *


def load_data():

	file = 'data.csv'

	if 'outfile.csv' in file == True:

		os.remove('outfile.csv') 

	data = open(file, 'r')

	name = []
	seq = [] 

	for line in data:
		
		sequence = line.split(',')[1].strip('\n').upper()
		strand_name = line.split(',')[0]

		seq.append(sequence)
		name.append(strand_name)

	return name, seq 

def write_data(name, seq):

	rows = []

	for i in range(len(name)):

		rows.append([name[i],
					 seq[i],
					 complement(seq[i]),
					 len(seq[i]),
					 # cvrt_to_binary(seq[i]),
					 gc_content(seq[i]),
					 molecular_weight(seq[i]),
					 extinction_coeff(seq[i])
					 melting_temp(seq[i],2)])


	with open('outfile.csv', 'w') as f:
		csv_writer = writer(f)
		csv_writer.writerow(['name', 'sequence', 'complement', 'length', 'gc content (%)', 'molecular weight (g/mol)', 'extinction (M-1 cm-1)', 'Tm (C)'])
		for row in rows:
			csv_writer.writerow(row)


name, seq = load_data()

dna = 'GTGAGTAGGTAGAGA'

print(f'DNA sequence {dna} \nHas {len(dna)} bases.')
melting_temp(dna, 2, 0, 0, 100, 50)
