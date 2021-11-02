import math

def flip_seq(seq):

	rev = seq[::-1]

	return rev

def complement(seq):

	comp = '' 

	for char in seq:

		if char == 'A':
			comp += 'T'

		elif char == 'T':
			comp += 'A'

		elif char == 'G':
			comp += 'C'

		elif char == 'C':
			comp += 'G'						

	return comp

def base_frequency(seq):
	wA = seq.count('A')
	xT = seq.count('T')
	yC = seq.count('C')
	zG = seq.count('G')

	return wA, xT, yC, zG

def gc_content(seq):

	_, _, yC, zG = base_frequency(seq)

	gc_content = round(((zG + yC) / len(seq)) * 100)

	return gc_content

def molecular_weight(seq):

	wA, xT, yC, zG = base_frequency(seq)

	mol = (wA * 313.21) + (xT * 304.2) + (yC * 289.18) + (zG * 329.21) - 61.96

	return mol

def extinction_coeff(seq):

	dada = 27400
	dadc = 21200
	dadg = 25000
	dadt = 22800

	dcda = 21200
	dcdc = 14600
	dcdg = 18000
	dcdt = 15200

	dgda = 25200
	dgdc = 17600
	dgdg = 21600
	dgdt = 20000

	dtda = 23400
	dtdc = 16200
	dtdg = 19000
	dtdt = 16800

	da = 15400
	dc = 7400
	dg = 11500
	dt = 8700

	individualSum = 0
	individualBases = []
	for i, _ in enumerate(seq):


		if i == 0:
			pass

		elif i == len(seq)-1:
			pass

		else:
			if seq[i] == 'A':
				individualSum += da
				individualBases.append('A')
			if seq[i] == 'C':
				individualSum += dc
				individualBases.append('C')
			if seq[i] == 'G':
				individualSum += dg
				individualBases.append('G')
			if seq[i] == 'T':
				individualSum += dt
				individualBases.append('T')

	neighborSum = 0
	neighborBases = []
	for i, _ in enumerate(seq):	

		try:
			if seq[i] == 'A' and seq[i+1] == 'A':
				neighborSum += dada
				neighborBases.append('AA')
			elif seq[i] == 'A' and seq[i+1] == 'C':
				neighborSum += dadc
				neighborBases.append('AC')
			elif seq[i] == 'A' and seq[i+1] == 'G':
				neighborSum += dadg
				neighborBases.append('AG')
			elif seq[i] == 'A' and seq[i+1] == 'T':
				neighborSum += dadt
				neighborBases.append('AT')

			elif seq[i] == 'C' and seq[i+1] == 'A':
				neighborSum += dcda
				neighborBases.append('CA')
			elif seq[i] == 'C' and seq[i+1] == 'C':
				neighborSum += dcdc
				neighborBases.append('CC')
			elif seq[i] == 'C' and seq[i+1] == 'G':
				neighborSum += dcdg
				neighborBases.append('CG')
			elif seq[i] == 'C' and seq[i+1] == 'T':
				neighborSum += dcdt
				neighborBases.append('CT')

			elif seq[i] == 'G' and seq[i+1] == 'A':
				neighborSum += dgda
				neighborBases.append('GA')
			elif seq[i] == 'G' and seq[i+1] == 'C':
				neighborSum += dgdc
				neighborBases.append('GC')
			elif seq[i] == 'G' and seq[i+1] == 'G':
				neighborSum += dgdg
				neighborBases.append('GG')
			elif seq[i] == 'G' and seq[i+1] == 'T':
				neighborSum += dgdt
				neighborBases.append('GT')

			elif seq[i] == 'T' and seq[i+1] == 'A':
				neighborSum += dtda
				neighborBases.append('TA')
			elif seq[i] == 'T' and seq[i+1] == 'C':
				neighborSum += dtdc
				neighborBases.append('TC')
			elif seq[i] == 'T' and seq[i+1] == 'G':
				neighborSum += dtdg
				neighborBases.append('TG')
			elif seq[i] == 'T' and seq[i+1] == 'T':
				neighborSum += dtdt
				neighborBases.append('TT')

		except IndexError:
			pass

	ext = (neighborSum - individualSum)
	
	return ext

def melting_temp(seq, dnac, na, k, tris, mg):
	
	dh, ds = nn_thermodynamics(seq)

	k = (dnac/2) * 1e-6
	R = 1.987
	Tm = (1000 * dh)/(ds + (R * math.log(k))) - 273.15

	corr = salt_correction(seq, na, k, tris, mg)
	melting_tm = 1 / (1 / (Tm + 273.15) + corr) - 273.15

	return round(melting_tm, 2), round(Tm, 2)

def nn_thermodynamics(seq):

	# Code adapted from biopython. Open-source but check copyright code before use
	seq = seq
	comp = complement(seq)
	dh_sum = 0
	ds_sum = 0
	d_h = 0  # Names for indexes
	d_s = 1  # 0 and 1

	# Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594
	nn_table =	SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
	{"init": (0.2, -5.7),
	"init_A/T": (2.2, 6.9),
	"init_G/C": (0, 0),
	"init_oneG/C": (0, 0),
	"init_allA/T": (0, 0),
	"init_5T/A": (0, 0),
	"sym": (0, -1.4),
	"AA/TT": (-7.6, -21.3),
	"AT/TA": (-7.2, -20.4),
	"TA/AT": (-7.2, -20.4),
	"CA/GT": (-8.5, -22.7),
	"GT/CA": (-8.4, -22.4),
	"CT/GA": (-7.8, -21.0),
	"GA/CT": (-8.2, -22.2),
	"CG/GC": (-10.6, -27.2),
	"GC/CG": (-9.8, -24.4),
	"GG/CC": (-8.0, -19.0)}

	# {
 #    "init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
 #    "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
 #    "sym": (0, -1.4),
 #    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
 #    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
 #    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
 #    "GG/CC": (-8.0, -19.9)} 

	# Type: General initiation value
	dh_sum += nn_table["init"][d_h]
	ds_sum += nn_table["init"][d_s]

	# Type: Penalty if 5' end is T
	if seq[0] == ("T"):
		dh_sum += nn_table["init_5T/A"][d_h]
		ds_sum += nn_table["init_5T/A"][d_s]
	if seq[-1] == ("A"):
		dh_sum += nn_table["init_5T/A"][d_h]
		ds_sum += nn_table["init_5T/A"][d_s]

	# Type: Different values for G/C or A/T terminal basepairs
	ends = seq[0] + seq[-1]
	AT = ends.count("A") + ends.count("T")
	GC = ends.count("G") + ends.count("C")
	dh_sum += nn_table["init_A/T"][d_h] * AT
	ds_sum += nn_table["init_A/T"][d_s] * AT
	dh_sum += nn_table["init_G/C"][d_h] * GC
	ds_sum += nn_table["init_G/C"][d_s] * GC

	for basenumber in range(len(seq) - 1):
		neighbors = (
			seq[basenumber : basenumber + 2]
			+ '/'
			+ comp[basenumber : basenumber +2]
		)

		if neighbors in nn_table:
			dh_sum += nn_table[neighbors][d_h]
			ds_sum += nn_table[neighbors][d_s]
		elif neighbors[::-1] in nn_table:
			dh_sum += nn_table[neighbors[::-1]][d_h]
			ds_sum += nn_table[neighbors[::-1]][d_s]	
		
		print(f'With NN pair {neighbors}: Delta H = {round(dh_sum,2)} kcal/mol -- Delta S = {round(ds_sum,2)} cal/mol*k')				

	dh_sum = round(dh_sum,2)
	ds_sum = round(ds_sum,2)

	print(f'The final Delta H = {dh_sum} and Delta S = {ds_sum}')

	return dh_sum, ds_sum

def salt_correction(seq, na, k, tris, mg):

	# assuming pH 7 and milimolar units
	# correction used is Owczarzy 2008

	mon = (na + k + (tris / 2.0)) * 1e-3
	di = mg * 1e-3

	a = 3.92e-5
	b = -9.11e-6
	c = 6.26e-5
	d = 1.42e-5
	e = -4.82e-4
	f = 5.25e-4
	g = 8.31e-5

	if mon > 0:
		R = math.sqrt(di) / mon
		if R < 0.22:
			corr = (4.29 * (gc_content(seq) / 100) - 3.95) * 1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2
			return corr
		
		elif R < 6.0:
			a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
			d = 1.42e-5 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
			g = 8.31e-5 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)

	corr = (a + b * math.log(mg) + (gc_content(seq) / 100) * (c + d * math.log(mg)) + 1 / (2.0 * (len(seq) - 1)) * (e + f * math.log(mg) + g * (math.log(mg) ** 2)))
	return corr

def cvrt_to_binary(seq):

	binary_seq = ''

	for char in seq:

		if char == 'A':
			binary_seq += '00'
		elif char == 'C':
			binary_seq += '01'
		elif char == 'G':
			binary_seq += '10'
		elif char == 'T':
			binary_seq += '11'

	return binary_seq

def cvrt_to_nucleo(seq):

	end = len(seq)

	nucleo_seq = ''

	for i in range(0, end-1, 2):

		if seq[i] + seq[i+1] == '00':
			nucleo_seq += 'A'
		elif seq[i] + seq[i+1] == '01':
			nucleo_seq += 'C'
		elif seq[i] + seq[i+1] == '10':
			nucleo_seq += 'G'
		elif seq[i] + seq[i+1] == '11':
			nucleo_seq += 'T'

	return nucleo_seq

def check_repeats(strand_pool):

	unique = []
	
	for strand in seq:

		if (strand in unique) == False:
			unique.append(strand)
			
		else:
			print(f'{strand} is repeated')
