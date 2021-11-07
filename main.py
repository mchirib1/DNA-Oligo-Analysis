import os
from csv import writer

from utils.utils import *


def load_data(infile):

	if 'outfile.csv' in  os.listdir(os.getcwd()) == True:

		os.remove('outfile.csv') 

	data = open(infile, 'r')

	seq_names = []
	seq_pool = [] 

	for line in data:
		
		sequence = line.split(',')[1].strip('\n').upper()
		strand_name = line.split(',')[0]

		seq_pool.append(sequence)
		seq_names.append(strand_name)

	return seq_names, seq_pool 

def write_data(seq_names, seq_pool, dnac, Na, K, Tris, Mg):
	# DNA concentration in uM and Ion concentrations in mM 

	rows = []

	for i in range(len(seq_names)):

		rows.append([seq_names[i],
					 seq_pool[i],
					 complement(seq_pool[i]),
					 len(seq_pool[i]),
					 cvrt_to_binary(seq_pool[i]),
					 gc_content(seq_pool[i]),
					 molecular_weight(seq_pool[i]),
					 extinction_coeff(seq_pool[i]),
					 melting_temp(seq_pool[i], dnac, Na, K, Tris, Mg)[0],
					 nn_thermodynamics(seq_pool[i])[0],
					 nn_thermodynamics(seq_pool[i])[1]
					 ])


	with open('outfile.csv', 'w') as f:
		csv_writer = writer(f)
		csv_writer.writerow(['name',
							 'sequence',
							 'complement',
							 'length',
							 'binary',
							 'gc content (%)',
							 'molecular weight (g/mol)',
							 'extinction (M-1 cm-1)',
							 'salt adjusted melting temp (C)',
							 'delta H kcal/mol',
							 'delta S cal/(mol*k)'
							 ])
		for row in rows:
			csv_writer.writerow(row)


name, seq_pool = load_data('data.csv')
write_data(name, seq_pool, 2.5, 5, 0, 10, 1.2)
check_repeats(seq_pool)

