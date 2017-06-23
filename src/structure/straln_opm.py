# Name: straln_opm.py
# Language: python3
# Libraries: 
# Description: Runs the structural alignments needed for HOMEP
# Author: Edoardo Sarti
# Date: Aug 17 2016

import os
import sys
import multiprocessing
import subprocess
import re
import time
from support_opm import *
from QuickPylignMe import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo


def fast_blosum(blosum62_address, seq_1, seq_2):
	this_name = 'fast_blosum'
	seqname_1 = seq_1[0]
	seqname_2 = seq_2[0]
	sequence_1 = seq_1[1]
	sequence_2 = seq_2[1]

	RE = RegExp()
	f = open(blosum62_address, 'r')
	text = f.read().split('\n')
	matrix = {}
	for nl in range(len(text)):
		fields = text[nl].split()
		if RE.emptyline.match(text[nl]):
			continue
		elif (len(fields) == 20):
			col_names = []
			for nf in range(len(fields)):
				col_names.append(fields[nf])
		elif (len(fields) == 21):
			row_name = str(fields[0])
			for nf in range(1, len(fields)):
				matrix[(row_name, col_names[nf-1])] = float(fields[nf])

	potential = np.zeros((len(sequence_1), len(sequence_2)), dtype=[('D', 'f'),('H', 'f'), ('V', 'f')])
	for na1 in range(len(sequence_1)):
		for na2 in range(len(sequence_2)):
			potential[na1, na2]['D'] += matrix[(sequence_1[na1], sequence_2[na2])]

	profile = []
	for i in range(len(sequence_1)):
		profile.append([])
		for j in range(len(sequence_2)):
			profile[i].append('')
	
	gap = Gaps()
	NWM = needleman_wunsch(sequence_1, sequence_2, potential, profile, gap)
	alignments, scores = traceback(sequence_1, sequence_2, NWM)
	# Choose only the alignments that reach (0,0)
	full_alignments = []
	print(sequence_1, sequence_2)
	for p in alignments:
		if alignments[p][-1] == (0,0):
			full_alignments.append((list(reversed(alignments[p])), scores[p]))
#			aligned_seqs, matches, signatures = to_fasta((sequence_1, sequence_2), list(reversed(alignments[p])))
#			print("Score {0}:\n> {1}:\n{2}\n{3}\n> {4}:\n{5}\n{6}\n".format(scores[p], 'name1', ''.join(str(p) for p in aligned_seqs[0]), signatures[0], 'name2', ''.join(str(p) for p in aligned_seqs[1]), signatures[1]))

	full_alignments = sorted(full_alignments, key=lambda fa: -fa[1])
	best_aln, best_score = full_alignments[0]

	aligned_seqs, matches, signatures = to_fasta((sequence_1, sequence_2), best_aln)

	print_log(this_name, "Score {0}:\n> {1}:\n{2}\n{3}\n> {4}:\n{5}\n{6}\n".format(best_score, seqname_1, ''.join(str(p) for p in aligned_seqs[0]), signatures[0], seqname_2, ''.join(str(p) for p in aligned_seqs[1]), signatures[1]))

	return (aligned_seqs, best_score)


def FrTMjob(data):
	tic = time.time()
	on_cluster = False
	if len(data) == 2:
		on_cluster = True
		dataname, task = data
		data = read_data(dataname)

	locations, target, exelist = data
	topologytype, topology, chain_1 = target
	this_name = 'FrTMjob_{0}'.format(chain_1)
	straln_path = locations['OPT']['stralnexe']
	seqaln_path = locations['OPT']['seqalnexe']
	
	topology_path = locations['FSYSPATH'][topologytype] + topology + '/'
	strfasta_filename = topology_path + locations['TREE']['straln'] + 'str_' + chain_1 + '_fasta.dat'
	strpdb_filename = topology_path + locations['TREE']['straln'] + 'str_' + chain_1 + '_pdb.dat'
	seqfasta_filename = topology_path + locations['TREE']['seqaln'] + 'seq_' + chain_1 + '_fasta.dat'
	str_repo_path = locations['FSYSPATH']['repocstraln']
	seq_repo_path = locations['FSYSPATH']['repocseqaln']

	fasta_tmpfolder_path = topology_path + locations['TREE']['straln'] + 'tmp_' + chain_1 + '_fasta/'
	pdb_tmpfolder_path = topology_path + locations['TREE']['straln'] + 'tmp_' + chain_1 + '_pdb/'
	if (on_cluster and task == 'from_scratch' and (os.path.exists(strfasta_filename) or os.path.exists(strpdb_filename) or os.path.exists(seqfasta_filename) or
	    os.path.exists(fasta_tmpfolder_path) or os.path.exists(pdb_tmpfolder_path))):
		# Someone else is taking care of this, so shut down
		print_log(this_name, "Start {0}".format(chain_1))
		print_log(this_name, "Someone else is taking care of this, so shut down")
		print_log(this_name, "Stop {0}".format(chain_1))
		return

	if on_cluster and task == 'continue':
		print_log(this_name, "Calculation to be continued for {0}".format(chain_1))
		fasta_tmpfolder_path_occ_filename = fasta_tmpfolder_path + 'occupied'
		pdb_tmpfolder_path_occ_filename = pdb_tmpfolder_path + 'occupied'
		if os.path.exists(fasta_tmpfolder_path_occ_filename) or os.path.exists(pdb_tmpfolder_path_occ_filename):
			print_log(this_name, "Found occupied tag in folder {0} or {1}".format(fasta_tmpfolder_path_occ_filename, pdb_tmpfolder_path_occ_filename))
			# Someone else is taking care of this, so shut down
			print_log(this_name, "Start {0}".format(chain_1))
			print_log(this_name, "Someone else is taking care of this, so shut down")
			print_log(this_name, "Stop {0}".format(chain_1))
			return
		else:
			fasta_tmpfolder_path_occ_file = open(fasta_tmpfolder_path_occ_filename, 'w')
			subprocess.call(["echo", "{0}".format(this_name)], stdout=fasta_tmpfolder_path_occ_file)
			fasta_tmpfolder_path_occ_file.close()
			pdb_tmpfolder_path_occ_file = open(pdb_tmpfolder_path_occ_filename, 'w')
			subprocess.call(["echo", "{0}".format(this_name)], stdout=pdb_tmpfolder_path_occ_file)
			pdb_tmpfolder_path_occ_file.close()

	# Checks for needed locations:
	# Superfamily path
	if not os.path.exists(topology_path):
		raise NameError("ERROR: Superfamily {0} not found in path {1}".format(topologytype+' '+topology, topology_path))
	# structure/ folder
	if not os.path.exists(topology_path + locations['TREE']['str']):
		raise NameError("ERROR: {0} folder not found in path {1}".format(locations['TREE']['str'], topology_path))
	# Main pdb file
	if not os.path.exists(topology_path + locations['TREE']['str'] + chain_1 + '.pdb'):
		raise NameError("ERROR: File {0} not found in {1}".format(chain_1 + '.pdb', topology_path + locations['TREE']['str']))
	# Secondary pdb files
	for chain_2 in exelist:
		if not os.path.exists(topology_path + locations['TREE']['str'] + chain_2 + '.pdb'):
			raise NameError("ERROR: File {0} not found in {1}".format(chain_2 + '.pdb', topology_path + locations['TREE']['str']))
	
	print_log(this_name, "Start {0}".format(chain_1))
	sys.stdout.flush()

	if (on_cluster and task == 'from_scratch') or not on_cluster:
		# Creates the temporary folder for sequence alignments
		os.mkdir(fasta_tmpfolder_path)
		# Creates the temporary folder for structure alignments
		os.mkdir(pdb_tmpfolder_path)

	repo_info_str_pdb = {chain_1 : repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_pdb.dat')[1]}
	repo_info_str_fasta = {chain_1 : repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_fasta.dat')[1]}
	repo_info_seq_fasta = {chain_1 : repo_inspector(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat')[1]}

	pdb1_filename = topology_path + locations['TREE']['str'] + chain_1 + '.pdb'
	for chain_2 in exelist:
		pdb2_filename = topology_path + locations['TREE']['str'] + chain_2 + '.pdb'
		if chain_2 in repo_info_str_fasta[chain_1] and chain_2 in repo_info_str_pdb[chain_1]:
			print_log(this_name, "Alignment {0} {1} already present in repository".format(pdb1_filename[-10:], pdb2_filename[-10:]))
			continue

		# Defines filenames
		fasta_output_filename = fasta_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_fasta.tmp'
		pdb_output_filename = 'straln_' + chain_1 + '_' + chain_2 + '_pdb.tmp' # This filename is without path for a constraint on flags length of frtmalign (see below)
		stdout_filename = fasta_tmpfolder_path + 'output_' + chain_1 + '_' + chain_2 + '.tmp'


		# If sequence and structure temporary files are present, jump
		if os.path.exists(fasta_output_filename) and os.path.exists(pdb_output_filename):
			print_log(this_name, "Alignment {0} {1} already present in temporary folder".format(pdb1_filename[-10:], pdb2_filename[-10:]))
			continue

		# Fr-TM-align call
		# The Fr-TM-align command cannot exceed a certain length. Thus, we are compelled to run it locally into structures/
		# The stdout file from Fr-TM-align can go directly to the right temporary folder (but needs to be cleaned)
		stdout_file = open(stdout_filename, 'w')
		fnull = open(os.devnull, 'w')
		print_log(this_name, "{0} {1} {2} -o {3}".format(straln_path, pdb1_filename[-10:], pdb2_filename[-10:], pdb_output_filename))
		p = subprocess.Popen([straln_path, pdb1_filename[-10:], pdb2_filename[-10:], '-o', pdb_output_filename], stdout=stdout_file, stderr=fnull, cwd=topology_path+locations['TREE']['str'])
		p.wait()
#		print_log(this_name, "Alignment done")
		fnull.close()
		stdout_file.close()

		# Moves the Fr-TM-align output file into the structure temporary folder
#		print_log(this_name, "Rename output file")
		os.rename(topology_path + locations['TREE']['str'] + pdb_output_filename, pdb_tmpfolder_path + pdb_output_filename)

		# Reads and then removes the stdout file from Fr-TM-align
#		print_log(this_name, "Reads and removes output file")
		stdout_file = open(stdout_filename, 'r')
		text = stdout_file.read().split('\n')
		stdout_file.close()
		os.remove(stdout_filename)
		
		# From the stdout file from Fr-TM-align, it takes the RMSD, TM-score and the two aligned sequences
#		print_log(this_name, "Parses output")
		chkaln = -1000
		for nl in range(len(text)):
			if "Aligned length" in text[nl]:
				fields = re.split('=|,|\s',text[nl])
				fields = list(filter(None, fields))
				RMSD = float(fields[4])
				TMscore = float(fields[6])
			elif chkaln+1 == nl:
				seq_1 = text[nl]
			elif chkaln+3 == nl:
				seq_2 = text[nl]
			elif "denotes the residue pairs of distance" in text[nl]:
				chkaln = nl
		print_log(this_name, "FrTMAlign output (RMDS, TMscore): {0}  {1}".format(RMSD, TMscore))
		
		# Creates a sequence temporary file already correctly formatted
#		print_log(this_name, "Creates temporary file for sequence")
		fasta_output_file = open(fasta_output_filename, 'w')
		fasta_output_file.write(">" + chain_1 + "\n" + seq_1.replace('\x00', '') + "\n>" + chain_2 + "\n" + seq_2.replace('\x00', '') + "\n\nRMSD\t{0:.2f}\nTM-score\t{1:.5f}\nstr_SEQID\t{2:.5f}\n\n".format(RMSD, TMscore, calculate_seqid((seq_1, seq_2))))
		fasta_output_file.close()

	fasta1_filename = topology_path + locations['TREE']['seq'] + chain_1 + '.fa'
	fasta1_file = open(fasta1_filename, 'r')
	text_f1 = fasta1_file.read().split('\n')
	for line in text_f1:
		if line and line[0] != '>':
			sequence_1 = line.strip()
			break
	fasta1_file.close()

	print_log(this_name, "Performing sequence alignments from sequence")
	seqfasta_file = open(seqfasta_filename, 'w')

	for chain_2 in exelist:
		seqfasta_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not chain_2 in repo_info_seq_fasta[chain_1]:
			fasta2_filename = topology_path + locations['TREE']['seq'] + chain_2 + '.fa'
			fasta2_file = open(fasta2_filename, 'r')
			text = text_f1 + fasta2_file.read().split('\n')
			fasta2_file.close()
			tmpmerge_filename = topology_path + locations['TREE']['seqaln'] + 'in_' + chain_1 + '_' + chain_2 + '.fa'
			tmpmerge_file = open(tmpmerge_filename, 'w')
			for line in text:
				tmpmerge_file.write(line+'\n')
			tmpmerge_file.close()
	
#			result = fast_blosum(blosum_address, (chain_1, sequence_1), (chain_2, sequence_2))
			tries = 0
			tmpout_filename = topology_path + locations['TREE']['seqaln'] + 'out_' + chain_1 + '_' + chain_2 + '.fa'
			while (not os.path.exists(tmpout_filename)):
				if tries >= 30:
					raise NameError("There is something wrong with this: {0}".format(tmpmerge_filename))
				fnull = open(os.devnull, 'w')
				p = subprocess.Popen([seqaln_path, '-in', tmpmerge_filename, '-out', tmpout_filename], stdout=fnull, stderr=fnull)
				p.wait()
				fnull.close()
				tries += 1

			tmpout_file = open(tmpout_filename, 'r')
			text = tmpout_file.read().split('\n')
			tmpout_file.close()

			os.remove(tmpmerge_filename)
			os.remove(tmpout_filename)
			cc = -1
			seqaln = {}
			for line in text:
				if not line:
					continue
				if line[0] == '>':
					cc += 1
					seqaln[cc] = ''
				else:
					seqaln[cc] += line.strip()
			if cc > 1:
				print_log(this_name, "Error: more than two sequences in the seq_aln between {0} and {1}:\n{2}".format(chain_1, chain_2, text))

			text = ">{0}\n{1}\n>{2}\n{3}\n\nseq_SEQID {4}\n".format(chain_1, seqaln[0], chain_2, seqaln[1], calculate_seqid(seqaln))
			seqfasta_file.write(text)
			repo_info_seq_fasta[chain_1][chain_2] = text
		else:
			for line in repo_info_seq_fasta[chain_1][chain_2].split('\n'):
				seqfasta_file.write(line+"\n")
		seqfasta_file.write("END\n\n\n")

	write_on_repo(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat', repo_info_seq_fasta)

	# Writes on the main sequence file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	print_log(this_name, "Writing sequence alignments from structure")
	strfasta_file = open(strfasta_filename, 'w')
	for chain_2 in exelist:
		strfasta_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not chain_2 in repo_info_str_fasta[chain_1]:
			tmp_filename = fasta_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_fasta.tmp'
			if not os.path.exists(tmp_filename):
				continue
			tmp_file = open(tmp_filename)
			text = tmp_file.read().split('\n')
			tmp_file.close()
			repo_info_str_fasta[chain_1][chain_2] = ""
			for line in text:
				strfasta_file.write(line+'\n')
				repo_info_str_fasta[chain_1][chain_2] += line + '\n'
		else:
			for line in repo_info_str_fasta[chain_1][chain_2].split('\n'):
				strfasta_file.write(line + '\n')
		strfasta_file.write("END\n\n\n")
	write_on_repo(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_fasta.dat', repo_info_str_fasta)
	strfasta_file.close()
	time.sleep(1)
	
	# Writes on the main structure file. Each alignment begins with a "BEGIN" line, followed by two lines reporting the two chain names
	# (format: "CHAIN_X: chain_name", where X={1,2}), and ends with an "END" line.
	print_log(this_name, "Writing structure alignments")
	strpdb_file = open(strpdb_filename, 'w')
	for chain_2 in exelist:
		strpdb_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n")
		if not chain_2 in repo_info_str_pdb[chain_1]:
			tmp_filename = pdb_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_pdb.tmp'
			if not os.path.exists(tmp_filename):
				continue
			tmp_file = open(tmp_filename)
			text = tmp_file.read().split('\n')
			tmp_file.close()
			repo_info_str_pdb[chain_1][chain_2] = ""
			for line in text:
				strpdb_file.write(line+'\n')
				repo_info_str_pdb[chain_1][chain_2] += line + '\n'
		else:
			for line in repo_info_str_pdb[chain_1][chain_2].split('\n'):
				strpdb_file.write(line + '\n')
		strpdb_file.write("END\n\n\n")
	write_on_repo(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_pdb.dat', repo_info_str_pdb)
	strpdb_file.close()
	time.sleep(1)

	# Delete temporary files and folders
	for chain_2 in exelist:
		tmp_filename = fasta_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_fasta.tmp'
		if os.path.exists(tmp_filename):
			os.remove(tmp_filename)
		tmp_filename = pdb_tmpfolder_path + 'straln_' + chain_1 + '_' + chain_2 + '_pdb.tmp'
		if os.path.exists(tmp_filename):
			os.remove(tmp_filename)
	os.rmdir(fasta_tmpfolder_path)
	os.rmdir(pdb_tmpfolder_path)

	toc = time.time()

	if on_cluster and task == 'continue':
		if os.path.exists(fasta_tmpfolder_path_occ_filename):
			os.remove(fasta_tmpfolder_path_occ_filename)
		if os.path.exists(pdb_tmpfolder_path_occ_filename):
			os.remove(pdb_tmpfolder_path_occ_filename)

	print_log(this_name, "Stop {0}".format(chain_1))
	sys.stdout.flush()
	print_log(this_name, "Time spent: {0}".format(toc-tic))


def calculate_seqid(alignment):
	ntot = 0
	naln = 0
	for na in range(len(alignment[0])):
		if alignment[0][na] != '-' and alignment[1][na] != '-':
			ntot += 1
		if alignment[0][na] == alignment[1][na]:
			naln += 1
	if ntot > 0:
		return naln/ntot
	else:
		return 0


def make_new_table(locations, external_filename, instraln=False, equivalence={}, allow_incomplete=False):
	this_name = 'make_new_table'
	tic = time.time()
	names = {'alpha' : 'a', 'beta' : 'b'}

	if instraln:
		instructions_filename = locations['SYSFILES']['H_instraln']
	else:
		instructions_filename = locations['SYSFILES']['H_topologytype']
	instructions_file = open(instructions_filename, 'r')
	text = instructions_file.read().split('\n')
	instructions_file.close()

	instructions = {}
	for line in text:
		if not line:
			continue
		fields = line.split()
		if fields[0] not in instructions:
			instructions[fields[0]] = {}
		if int(fields[1]) not in instructions[fields[0]]:
			instructions[fields[0]][int(fields[1])] = []
		instructions[fields[0]][int(fields[1])].append(fields[2])

	table_filename = locations['FSYSPATH']['main'] + external_filename
	table_file = open(table_filename, 'w')
	table = {}
	not_completed = set()
	for toptype in sorted(list(instructions.keys())):
		table[toptype] = {}
		for top in sorted(list(instructions[toptype].keys())):
			table[toptype][top] = {}
			for chain_1 in sorted(instructions[toptype][top]):
				table[toptype][top][chain_1] = {}
				if chain_1 in equivalence:
					print_log(this_name, "Setting up equivalence restore for struct {0} from struct {1}".format(chain_1, equivalence[chain_1]))
					chain_equiv = equivalence[chain_1]
					seq_fasta_repo = repo_inspector(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_equiv + '_fasta.dat')[1]
					str_fasta_repo = repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_equiv + '_fasta.dat')[1]
				else:
					seq_fasta_repo = repo_inspector(locations['FSYSPATH']['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat')[1]
					str_fasta_repo = repo_inspector(locations['FSYSPATH']['repocstraln'] + 'str_' + chain_1 + '_fasta.dat')[1]
				for chain_2 in sorted(instructions[toptype][top]):
					if chain_1 == chain_2:
						continue
					if (chain_2 in list(equivalence.keys()) and equivalence[chain_2] == chain_1 or
					    chain_1 in list(equivalence.keys()) and equivalence[chain_1] == chain_2):
						seq_seqid = 1.0
						str_seqid = 1.0
						TMscore = 1.0
						RMSD = 0.0
					else:
						try:
							text = seq_fasta_repo[chain_2].split('\n')
						except:
							if allow_incomplete:
								not_completed.add(chain_1[:4])
								continue
							else:
								print(chain_1, chain_2)
								print(seq_fasta_repo)
								raise NameError("Incomplete data!!!")
						for line in text:
							if not line:
								continue
							fields = line.split()
							if fields[0] == 'seq_SEQID':
								seq_seqid = float(fields[1])
						try:
							text = str_fasta_repo[chain_2].split('\n')
						except:
							if allow_incomplete:
								not_completed.add(chain_1[:4])
								continue
							else:
								print(chain_1, chain_2)
								print(str_fasta_repo)
								raise NameError("Incomplete data!!!")
						for line in text:
							if not line:
								continue
							fields = line.split()
							if fields[0] == 'str_SEQID':
								str_seqid = float(fields[1])
							elif fields[0] == 'TM-score':
								TMscore = float(fields[1])
							elif fields[0] == 'RMSD':
								RMSD = float(fields[1])

					table[toptype][top][chain_1][chain_2] = (seq_seqid, str_seqid, TMscore, RMSD)
					if not allow_incomplete:
						table_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.8f}\t{5:10.8f}\t{6:10.8f}\t{7:10.6f}\n".format(names[toptype], str(top).zfill(3), chain_1, chain_2, seq_seqid, str_seqid, TMscore, RMSD))

	if allow_incomplete:
		for toptype in sorted(list(instructions.keys())):
			for top in sorted(list(instructions[toptype].keys())):
				for chain_1 in sorted(instructions[toptype][top]):
					if chain_1[:4] in not_completed:
						continue
					print("Writing", chain_1)
					for chain_2 in sorted(instructions[toptype][top]):
						if chain_1 == chain_2:
							continue
						table_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.8f}\t{5:10.8f}\t{6:10.8f}\t{7:10.6f}\n".format(names[toptype], str(top).zfill(3), chain_1, chain_2,  table[toptype][top][chain_1][chain_2][0],  table[toptype][top][chain_1][chain_2][1],  table[toptype][top][chain_1][chain_2][2],  table[toptype][top][chain_1][chain_2][3]))
						
	table_file.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return table
							

def make_table_from_repo(options, locations):
	### ADD EQUIVALENCES!!! FILE IS ALREADY STORED IN REPO, JUST COPY FROM OTHER ROUTINE
	this_name = 'make_table_from_repo'
	tic = time.time()
	names = {'alpha' : 'a', 'beta' : 'b'}

	instructions_filename = options['repository'] + '/.topology_classification.dat'
	instructions_file = open(instructions_filename, 'r')
	text = instructions_file.read().split('\n')
	instructions_file.close()

	instructions = {}
	for line in text:
		if not line:
			continue
		fields = line.split()
		if fields[0] not in instructions:
			instructions[fields[0]] = {}
		if fields[1] not in instructions[fields[0]]:
			instructions[fields[0]][fields[1]] = []
		instructions[fields[0]][fields[1]].append(fields[2])


	repolocations_filename = options['repository'] + os.path.basename(locations['SYSFILES']['H_repolocations'])
	if not os.path.exists(repolocations_filename):
		raise_error(this_name, "ERROR: hidden repository locations files not found in path {0}".format(options['repository']))
	repolocations_file = open(repolocations_filename, 'r')
	text = repolocations_file.read().split('\n')
	extrepolocations = {}
	extrepolocations['repository'] = options['repository']
	for line in text:
		if not line:
			continue
		fields = line.split()
		extrepolocations[fields[0]] = extrepolocations['repository'] + fields[1]


	fake_database = {}
	table = {}
	completed_chains_list = sorted([x[:-4] for x in os.listdir(extrepolocations['repochains']) if x[-4:] == '.pdb'])
	print_log(this_name, "Completed_chains_list length: {0}".format(len(completed_chains_list)))
	print_log(this_name, "Completed_chains_list extract {0}".format(completed_chains_list[:10]))
	for toptype in sorted(list(instructions.keys())):
		table[toptype] = {}
		for top in sorted(list(instructions[toptype].keys())):
			table[toptype][top] = {}
			for chain_1 in sorted(instructions[toptype][top]):
				if chain_1 not in completed_chains_list:
					continue
				table[toptype][top][chain_1] = {}
				print_log(this_name, "Chain {0}".format(chain_1))
				seq_fasta_repo = repo_inspector(extrepolocations['repocseqaln'] + 'seq_' + chain_1 + '_fasta.dat')[1]
				str_fasta_repo = repo_inspector(extrepolocations['repocstraln'] + 'str_' + chain_1 + '_fasta.dat')[1]
				for chain_2 in sorted(instructions[toptype][top]):
					if chain_2 not in completed_chains_list or chain_1 == chain_2:
						continue
					if chain_2 not in seq_fasta_repo or chain_2 not in str_fasta_repo:
						del table[toptype][top][chain_1]
						for revch in sorted(list(table[toptype][top].keys())):
							if chain_1 in table[toptype][top][revch]:
								del table[toptype][top][revch][chain_1]
						break
					text = seq_fasta_repo[chain_2].split('\n')
					for line in text:
						if not line:
							continue
						fields = line.split()
						if fields[0] == 'seq_SEQID':
							seq_seqid = float(fields[1])
					text = str_fasta_repo[chain_2].split('\n')
					for line in text:
						if not line:
							continue
						fields = line.split()
						if fields[0] == 'str_SEQID':
							str_seqid = float(fields[1])
						elif fields[0] == 'TM-score':
							TMscore = float(fields[1])
						elif fields[0] == 'RMSD':
							RMSD = float(fields[1])
					table[toptype][top][chain_1][chain_2] = (seq_seqid, str_seqid, TMscore, RMSD)
#					table_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.8f}\t{5:10.8f}\t{6:10.8f}\t{7:10.6f}\n".format(names[toptype], str(int(top)).zfill(3), chain_1, chain_2, seq_seqid, str_seqid, TMscore, RMSD))
				struct = chain_1[:4]
				chain = chain_1[5]
				if struct not in fake_database:
					fake_database[struct] = [{} ,{}]
					fake_database[struct][1]['FROM_PDB'] = {}
					PDB_dict = PDB_parser(locations, struct)
					fake_database[struct][1]['FROM_PDB'] = PDB_dict
					fake_database[struct][1]['CHAIN'] = {}
				if chain not in fake_database[struct][1]['CHAIN']:
					fake_database[struct][1]['CHAIN'][chain] = [{}]
					fake_database[struct][1]['CHAIN'][chain][0]['TYPE'] = toptype
					fake_database[struct][1]['CHAIN'][chain][0]['NUM_TM'] = top
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return table, fake_database


def structure_alignment(options, locations):
	this_name = 'structure_alignment'
	external_filename = options['output_tab']

	already_processed = []
	ex_list = {}
	hidden_repository_filename = locations['SYSFILES']['repocstraln']
	if os.path.exists(hidden_repository_filename):
		hidden_repository_file = open(hidden_repository_filename, 'r')
		text = hidden_repository_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			already_processed.append((fields[2], fields[3]))
		hidden_repository_file.close()
	repository_filename = locations['FSYSPATH']['main'] + external_filename
	if os.path.exists(repository_filename):
		repository_file = open(repository_filename, 'r')
		text = repository_file.read().split("\n")
		for line in text:
			if not line:
				continue
			fields = line.split()
			if not (fields[2], fields[3]) in already_processed:
				already_processed.append((fields[2], fields[3]))
		external_filename = 'new_' + external_filename
		repository_file.close()
	
	topologies = []
	for tt in 'alpha', 'beta':
		for i in os.listdir(locations['FSYSPATH'][tt]):
			if re.match('^\d*$', str(i)):
				topologies.append((tt, i))

	ex_check = {}
	exelist_filename = locations['SYSFILES']['H_scheduledalns']
	if os.path.exists(exelist_filename):
		exelist_file = open(exelist_filename, 'r')
		text = exelist_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			if fields[2] not in ex_list:
				ex_list[fields[2]] = []
			ex_list[fields[2]].append(fields[3])
		exelist_file.close()

	instraln_list = []
	if (options['instraln']):
		instructions_straln_filename = options['instraln']
	else:
		instructions_straln_filename = locations['SYSFILES']['H_instraln']
	if os.path.exists(instructions_straln_filename):
		instructions_straln_file = open(instructions_straln_filename, 'r')
		text = instructions_straln_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			if fields[2] not in instraln_list:
				instraln_list.append(fields[2])
		instructions_straln_file.close()
	elif (not os.path.exists(instructions_straln_filename)) and options['instraln']:
		raise NameError("Error: forced instruction_straln_file {0} but file not found".format(instructions_straln_filename))
	else:
		instraln_list.append('nolist')


	exelist = {}
	exelist_file = open(exelist_filename, 'a')
	for t in topologies:
		structs = [x[-10:-4] for x in os.listdir(locations['FSYSPATH'][t[0]] + t[1] + '/' + locations['TREE']['str']) if x[-4:] == '.pdb']
		if len(structs) > 1:
			for s1 in structs:
				if (not 'nolist' in instraln_list) and s1 not in instraln_list:
					continue
				exelist[s1] = {}
				for s2 in structs:
					if s1 == s2 or (s1, s2) in already_processed:
						continue
					exelist[s1][s2] = (t[0], t[1], s1, s2)
					if s1 not in ex_list or s2 not in ex_list[s1]:
						exelist_file.write("{0}\t{1}\t{2}\t{3}\n".format(t[0], t[1], s1, s2))
	exelist_file.close()

	if 'hostname' in options and options['hostname']:	
		host = 'cluster'
		hostname = options['hostname']
	else:
		host = 'local'

	data = []
	if host == 'cluster':
		commands_filename = locations['FSYSPATH']['main'] + 'FrTMjob_commands.txt'
		commands_file = open(commands_filename, 'w')
		is_done_dict = {}
	for i in sorted(list(exelist.keys())):
		exesublist = []
		if not list(exelist[i].keys()):
			continue
		for j in sorted(list(exelist[i].keys())):
			exesublist.append(exelist[i][j][3])
		jtmp = sorted(list(exelist[i].keys()))[0]
		# Creates, if needed, the alignment locations
		for name, val in locations['TREE'].items():
			topology_path = locations['FSYSPATH'][exelist[i][jtmp][0]] + exelist[i][jtmp][1] + '/'
			if 'aln' in name and not os.path.exists(topology_path + val):
				os.mkdir(topology_path + val)
		if host == 'cluster':
			# writes the input data in a log file to be given to the cluster routine
			data_filename = locations['FSYSPATH']['logs'] + exelist[i][jtmp][2] + '.dat'
			write_data(data_filename, (locations, (exelist[i][jtmp][0], exelist[i][jtmp][1], exelist[i][jtmp][2]), exesublist))
			# Compiles the command file here
			topology_path = locations['FSYSPATH'][exelist[i][jtmp][0]] + exelist[i][jtmp][1] + '/'
			fasta_tmpfolder_path = topology_path + locations['TREE']['straln'] + 'tmp_' + exelist[i][jtmp][2] + '_fasta/'
			pdb_tmpfolder_path = topology_path + locations['TREE']['straln'] + 'tmp_' + exelist[i][jtmp][2] + '_pdb/'
			if os.path.exists(fasta_tmpfolder_path) or os.path.exists(pdb_tmpfolder_path):
				commands_file.write("python FrTMjob_wrapper.py {0} continue\n".format(data_filename))
				if os.path.exists(fasta_tmpfolder_path + 'occupied'):
					os.remove(fasta_tmpfolder_path + 'occupied')
				if os.path.exists(pdb_tmpfolder_path + 'occupied'):
					os.remove(pdb_tmpfolder_path + 'occupied')
			else:
				commands_file.write("python FrTMjob_wrapper.py {0} from_scratch\n".format(data_filename))
			# Initializes the checklist
			is_done_dict[exelist[i][jtmp][2]] = [topology_path, False]
		elif host == 'local':
			data.append((locations, (exelist[i][jtmp][0], exelist[i][jtmp][1], exelist[i][jtmp][2]), exesublist))
	if host == 'cluster':
		commands_file.close()

	seq_repo_path = locations['FSYSPATH']['repocseqaln']
	str_repo_path = locations['FSYSPATH']['repocstraln']
	if not os.path.exists(seq_repo_path):
		os.mkdir(seq_repo_path)
	if not os.path.exists(str_repo_path):
		os.mkdir(str_repo_path)

	if host == 'cluster':
		num_jobs = int(options['number_of_jobs'])
		workdir = os.getcwd() + '/'
		print_log(this_name, "Create and submit first {0} pbs files".format(num_jobs))
		job_out_filename = {}
		for i in range(num_jobs):
			jobname = 'FrTMAlign_'+str(i)
			pbs_filename = locations['SYSFILES']['FrTMpbs'](jobname)
			jobname, walltime_hrs, job_out_filename[jobname] = create_pbs(hostname, locations, commands_filename, pbs_filename, jobname=jobname)
			subprocess.call(["ssh", "{0}".format(hostname), "qsub", "{0}".format(pbs_filename)], cwd=workdir)
			time.sleep(10)
		time.sleep(300)

#		previous_jobname = ''
#		current_jobname = jobname
#		current_instrfilename = commands_filename
		while False in [x[1] for x in list(is_done_dict.values())]:
#			print_log(this_name, "is_done_dict: {0}".format(is_done_dict.items()))
			for struct in sorted(list(is_done_dict.keys())):
				if not is_done_dict[struct][1]:
					topology_path = is_done_dict[struct][0]
					strfasta_filename = topology_path + locations['TREE']['straln'] + 'str_' + struct + '_fasta.dat'
					strpdb_filename = topology_path + locations['TREE']['straln'] + 'str_' + struct + '_pdb.dat'
					seqfasta_filename = topology_path + locations['TREE']['seqaln'] + 'seq_' + struct + '_fasta.dat'
					if os.path.exists(strfasta_filename) and os.path.exists(strpdb_filename) and os.path.exists(seqfasta_filename):
						print_log(this_name, "Alignments of struct {0} completed".format(struct))
						is_done_dict[struct][1] = True

#			print_log(this_name, "After check, is_done_dict: {0}".format(is_done_dict.items()))
			if not False in [x[1] for x in list(is_done_dict.values())]:
				print_log(this_name, "All alignments completed")
				break

			print_log(this_name, "Something is still missing: {0}".format("".join([x+" " for x in list(is_done_dict.keys()) if is_done_dict[x][1] == False])))
			job_out_filename = job_control(locations, walltime_hrs, hostname, job_out_filename, is_done_dict, num_jobs, 'FrTMAlign')
			time.sleep(300)
			
	elif host == 'local':
		pool = multiprocessing.Pool(processes=int(options['number_of_procs']))
		pool_outputs = pool.map(FrTMjob, data)
		pool.close()
		pool.join()

	equivalence_filename = locations['SYSFILES']['equivalences']
	equivalence = {}
	if os.path.exists(equivalence_filename):
		equivalence_file = open(equivalence_filename, 'r')
		text = equivalence_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			equivalence[fields[0]] = fields[1]
		equivalence_file.close()	

	table = {}
	table = make_new_table(locations, external_filename, instraln=True, equivalence=equivalence)
	return table
