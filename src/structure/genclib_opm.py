# Name: genclib.py
# Language: python3
# Libraries:
# Description: Generates HOMEP chain library
# Author: Edoardo Sarti
# Date: Aug 15 2016

import os
import shutil
import re
import multiprocessing
import subprocess
import copy
import spur
import spur.ssh
import numpy as np
from itertools import chain as iterchain
from Bio.PDB import *
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

from support_opm import *


# OPM checker
# Checks if the biological assembly of OPM is the same as the one reported in the PDB.
# If if is not, the correct biological assembly will be run in PPM.
def OPM_checker(PDB_dict_safe, PDBaddress, OPMpdb_dict, OPM_dict):
	this_name = 'OPM_checker'
	PDB_dict = copy.deepcopy(PDB_dict_safe)
	parser = PDBParser(PERMISSIVE=True)
	structure = parser.get_structure('self', PDBaddress)
	ppb = PPBuilder()
	PDB_NoC_withcoord = len([x for x in structure[0].get_chains() if ppb.build_peptides(x) and x.id in PDB_dict['BIOLOGICAL_UNIT']])
#	PDB_NoC = len(PDB_dict['BIOMATRIX'])*len(PDB_dict['BIOLOGICAL_UNIT'])
	PDB_NoC = len(PDB_dict['BIOMATRIX'])*PDB_NoC_withcoord
	OPM_NoC = len(OPMpdb_dict['CHAINS'])
	print_log(this_name, "Number of tm chains (PDB, OPM): {0} {1}".format(PDB_NoC, OPM_NoC))
	for chain in OPM_dict['tmchains']:
		if chain not in OPMpdb_dict['CHAINS']:
			print_log(this_name, "Tried to find OPM_dict chain in OPM pdb structure and failed".format(chain))
			return False
	if len(PDB_dict['BIOMATRIX']) > 1:
		if PDB_NoC == OPM_NoC:
			print_log(this_name, "OPM is fine, but there are equivalent chains. And the lesser evil is to PPM it.")
			return False
		else:
			print_log(this_name, "WARNING: This will lead to a mistake")
			return False
	error_found = False
	for chain in sorted(OPMpdb_dict['CHAINS']):
		print_log(this_name, "OPM_checker {0}".format(chain))
		try:
			PDB_dict['BIOLOGICAL_UNIT'].remove(chain)
		except:
			print_log(this_name, "OPM_checker tried to remove {0} and did not succeed".format(chain))
			error_found = True
			break
	if PDB_dict['BIOLOGICAL_UNIT']:
		print_log(this_name, "Checked the voidness of {0} and it is not".format(PDB_dict['BIOLOGICAL_UNIT']))
		error_found = True
	if not error_found:
		print_log(this_name, "OPM_checker OK")
		return True
	print_log(this_name, "OPM_checker NOPE")
	return False
	

# Structure checker
def check_whole_structures(locations, filters, OPM_data, run_in_PPM):
	this_name = 'check_whole_structures'
	tic = time.time()
	exclusions_filename = locations['SYSFILES']['excludedwhole'] 
	already_excluded = []
	if os.path.exists(exclusions_filename):
		exclusions_file = open(exclusions_filename, 'r')
		text = exclusions_file.read().split('\n')
		for line in text:
			if line and line[:6]:
				already_excluded.append(line[:6])
		exclusions_file.close()
		exclusions_file = open(exclusions_filename, 'a')
	else:
		exclusions_file = open(exclusions_filename, 'w')

	new_OPM_data = {}
	for struct in sorted(list(OPM_data.keys())):
		PDB_dict = PDB_parser(locations['FSYSPATH']['PDBpdbs'], struct)

		exclude_struct = []
		struct_filename = locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb'

		# NMR check
		if not filters['NMR'] and PDB_dict['TECHNIQUE'] == 'NMR':
			exclude_struct.append('NMR structure')

		# Theoretical model check
		if not filters['THM'] and PDB_dict['TECHNIQUE'] == 'THEORETICAL':
			print_log(this_name, "Struct {0} is theoretical".format(struct))
			exclude_struct.append('Theoretical model')

		# Resolution check
		if 'resolution' in filters and PDB_dict['TECHNIQUE'] != 'THEORETICAL' and PDB_dict['TECHNIQUE'] != 'NMR':
			if PDB_dict['RESOLUTION'] > float(filters['resolution']):
				exclude_struct.append('Resolution is higher than {0}'.format(filters['resolution']))
			elif PDB_dict['RESOLUTION'] == 0:
				exclude_struct.append('No resolution information found'.format(filters['resolution']))

		# WARNING: This filter excludes many NMR structures, so it is part of an NMR filter. Though, it is NOT an optional filter, because this information
		# is essential for the determination of the whole structure file. If you want to generalize to the NMRs, you have to find a way to fill this parameter
		# (and only then remove this filter)
		if 'BIOLOGICAL_UNIT' not in PDB_dict:
			exclude_struct.append('Chain {0}: the PDB does not have information regarding the biological unit, needed in this case'.format(chain))
			if struct not in already_excluded:
				exclusions_file.write(struct + '\t\t' + exclude_struct[0] + '\n')
			run_in_PPM.discard(struct)
			continue

		FASTA_dict = FASTA_parser(locations, struct, PDB_dict)

		# If there is the OPM pdb file, the checks end here. Otherwise, the checks on the PDB pdb file continue...
		if not os.path.exists(struct_filename):
#			if 'BIOLOGICAL_UNIT' not in PDB_dict:
#				exclude_struct.append('Chain {0}: the PDB does not have information regarding the biological unit, needed in this case'.format(chain))
#				if struct not in already_excluded:
#					exclusions_file.write(struct + '\t\t' + exclude_struct[0] + '\n')
#				run_in_PPM.discard(struct)
#				continue
			for chain in sorted(PDB_dict['BIOLOGICAL_UNIT']):
				if chain not in PDB_dict['CHAINS']:
					exclude_struct.append('Chain {0} is in biological unit but does not have coordinates'.format(chain))
					if struct not in already_excluded:
						exclusions_file.write(struct + '\t\t' + exclude_struct[0] + '\n')
					run_in_PPM.remove(struct)
					continue
				# Holes greater than hole threshold check
				if int(filters['hole_thr']) > 0:
					for n in range(1,len(PDB_dict[chain]['RESIDS'])):
						if PDB_dict[chain]['RESIDS'][n] - PDB_dict[chain]['RESIDS'][n-1] > int(filters['hole_thr']):
							exclude_struct.append('Chain {0} contains hole longer than {1} residues'.format(chain, filters['hole_thr']))
							break
						if PDB_dict[chain]['RESIDS'][n] - PDB_dict[chain]['RESIDS'][n-1] < 1:
							exclude_struct.append('Chain {0} contains disordered or repeated residue indexes'.format(chain))
							break

				# Strange residues
				for res in PDB_dict[chain]['RESNAMES']:
					if res == '0':
						exclude_struct.append('Chain {0}: one or more residues with unsupported names'.format(chain))
						break
		else:
			OPMpdb_dict = PDB_parser(locations['FSYSPATH']['OPMpdbs'], struct)
			for chain in sorted(OPMpdb_dict['CHAINS']):
				if int(filters['hole_thr']) > 0:
					for n in range(1,len(OPMpdb_dict[chain]['RESIDS'])):
						if OPMpdb_dict[chain]['RESIDS'][n] - OPMpdb_dict[chain]['RESIDS'][n-1] > int(filters['hole_thr']):
							exclude_struct.append('Chain {0} contains hole longer than {1} residues'.format(chain, filters['hole_thr']))
							break
				for res in OPMpdb_dict[chain]['RESNAMES']:
					if res == '0':
						exclude_struct.append('Chain {0}: one or more residues with unsupported names'.format(chain))
						break
		print_log(this_name, "Exclusion?")
		print_log(this_name, ''.join([x+'\n' for x in exclude_struct]))
		if not exclude_struct:
			print_log(this_name, "Struct {0} not excluded".format(struct))
			OPM_check = None
			if OPM_data[struct]['is_representative']:
				if os.path.exists(locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb'):
					OPMpdb_dict = PDB_parser(locations['FSYSPATH']['OPMpdbs'], struct)
					PDBaddress = locations['FSYSPATH']['PDBpdbs'] + struct + '.pdb'
					OPM_check =  OPM_checker(PDB_dict, PDBaddress, OPMpdb_dict, OPM_data[struct])
					if not OPM_check:
						print_log(this_name, "Add to run_in_PPM struct {0} because OPM_checker said so".format(struct))
						FASTA_dict = FASTA_parser(locations, struct, PDB_dict, force_PDB=True)
						run_in_PPM.add(struct)
			new_struct_filename = locations['FSYSPATH']['whole'] + struct + '_opm.pdb'
			if os.path.exists(struct_filename) and not os.path.exists(new_struct_filename):
				shutil.copy(struct_filename, new_struct_filename)
			elif not os.path.exists(struct_filename):
				print_log(this_name, "Add to run_in_PPM struct {0} because this file was not found: {1}".format(struct, struct_filename))
				run_in_PPM.add(struct)
			if struct not in new_OPM_data:
				new_OPM_data[struct] = OPM_data[struct]
			print_log(this_name,"Struct {0} has FROM_PDB".format(struct))
			new_OPM_data[struct]['OPMcheck'] = OPM_check
			new_OPM_data[struct]['FROM_PDB'] = PDB_dict
			new_OPM_data[struct]['FASTA'] = FASTA_dict
			print_log(this_name, "Struct {0} FASTA chains: {1}".format(struct, ''.join([x+' ' for x in sorted(list(new_OPM_data[struct]['FASTA'].keys())) if len(x)==1])))
		else:
			print_log(this_name, "Excluded")
			if struct not in already_excluded:
				print_log(this_name, "Write in excluded file")
				exclusions_file.write(struct + '\t\t' + exclude_struct[0] + '\n')
				for nl in range(1, len(exclude_struct)):
					exclusions_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_struct[nl] + '\n')
			if struct in run_in_PPM:
				print_log(this_name, "Remove from PPM set")
				run_in_PPM.remove(struct)

	exclusions_file.close()

	for struct in sorted(list(new_OPM_data.keys())):
		if not new_OPM_data[struct]['is_representative']:
			if new_OPM_data[struct]['related_to'] in new_OPM_data:
				new_OPM_data[struct]['OPMcheck'] = new_OPM_data[new_OPM_data[struct]['related_to']]['OPMcheck']

	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return new_OPM_data, run_in_PPM


def generate_full_structure(address_1, address_2, OPM_dict):
	this_name = 'generate_full_structure'
	tic = time.time()
	rotmatrices = []
	tarrays = []
	print_log(this_name, "{0}".format(address_1))
	sys.stdout.flush()	
	for nmat in range(len(OPM_dict['FROM_PDB']['BIOMATRIX'])):
		rotmatrices.append(np.array([OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][0][:3],
		                              OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][1][:3],
		                              OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][2][:3]]))
		tarrays.append(np.array([OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][0][3],
		                          OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][1][3],
		                          OPM_dict['FROM_PDB']['BIOMATRIX'][nmat][2][3]]))

	io = PDBIO()
	ppb = PPBuilder()
	alphabetical = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'
	ac = 0
	new_structure = Structure.Structure(0)
	new_model = Model.Model(0)
	used_chains = set()
	for nmat in range(len(OPM_dict['FROM_PDB']['BIOMATRIX'])):
		parser = PDBParser(PERMISSIVE=True)
		structure = parser.get_structure('self', address_1)
		for chain in [x for x in structure[0].get_chains() if x.id in OPM_dict['FROM_PDB']['BIOLOGICAL_UNIT']]:
			if nmat == 0:
				print_log(this_name, "Already used {0}".format(chain.id))
				used_chains.add(chain.id)
				if not ppb.build_peptides(chain):
					continue
			else:
				if not ppb.build_peptides(chain):
					continue
				OPM_dict['equivalent_chains'] = {}   # Contains all added chains which are identical to an already existing chain in the original structure
				for chid in alphabetical:
					if chid not in used_chains:
						print_log(this_name, "Use {0}".format(chid))
						print_log(this_name, "FASTA from {0} to {1}".format(chid, chain.id))
						sys.stdout.flush()
						OPM_dict['FASTA'][chid] = OPM_dict['FASTA'][chain.id]
						OPM_dict['FROM_PDB'][chid] = OPM_dict['FROM_PDB'][chain.id]
						OPM_dict['equivalent_chains'][chid] = chain.id
						chain.id = chid
						used_chains.add(chid)
						break
			for residue in chain:
				atomlist = []
				for atom in residue:
					atomlist.append((atom, atom.get_name(), atom.get_altloc()))
				atomlist = sorted(atomlist, key = lambda x : (x[2], x[1]))
#				print(residue.get_resname())
#				print(atomlist)
				atomnamelist = set()
				for at in atomlist:
					if at[1] in atomnamelist:
						print(residue.get_resname(), at, "discarded")
						continue
					coord = np.array(at[0].get_coord())
					new_coord = np.dot(rotmatrices[nmat], coord) + tarrays[nmat]
					at[0].set_coord(list(new_coord))
					atomnamelist.add(at[1])
			chain.detach_parent()
			new_model.add(chain)
	new_structure.add(new_model)
	io.set_structure(new_structure)
	io.save(address_2)
	print_log(this_name, "At the end there are these chains: {0}".format(''.join([x+' ' for x in sorted(list(OPM_dict['FASTA'].keys())) if len(x)==1])))
	toc = time.time()
	print(OPM_dict['FASTA'])
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return OPM_dict


def PPM_runner(data):
	def check_and_copy(struct, tmp_folder):
		this_name = 'PPM_runner:check_and_copy'
		out_filename = struct + '_results.dat'
		outsub_filename = struct + '_datasub1'
		outpar_filename = struct + '_datapar1'
		outpdb_filename = struct + 'out.pdb'

		# Check if it is peripheral
		is_peripheral = False
		if os.stat(tmp_folder + outsub_filename).st_size == 0:
			is_peripheral = True

		# Copy on repo
		print_log(this_name, "Copies")
		tic = time.time()
		shutil.copy(tmp_folder + out_filename, locations['FSYSPATH']['repoppmresults'])
		shutil.copy(tmp_folder + outsub_filename, locations['FSYSPATH']['repoppmdatasub1'] + outsub_filename + '.dat')
		shutil.copy(tmp_folder + outpar_filename, locations['FSYSPATH']['repoppmdatapar1'] + outpar_filename + '.dat')
		if not is_peripheral:
			shutil.copy(tmp_folder + outpdb_filename, locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb')
	
		# Move files
		shutil.move(tmp_folder + out_filename, locations['FSYSPATH']['PPM'])
		shutil.move(tmp_folder + outsub_filename, locations['FSYSPATH']['PPM'] + outsub_filename + '.dat')
		shutil.move(tmp_folder + outpar_filename, locations['FSYSPATH']['PPM'] + outpar_filename + '.dat')
		if not is_peripheral:
			shutil.move(tmp_folder + outpdb_filename, locations['FSYSPATH']['whole'] + struct + '_opm.pdb')

		shutil.rmtree(tmp_folder)
		toc = time.time()
		print_log(this_name, "Time spent {0}".format(toc-tic))

		if not is_peripheral:
			return True
		else:
			print_log(this_name, "Is peripheral")
			return False


	# Is it going to run in local parallel or in cluster?
	tic = time.time()
	if data[0] == 'local':
		host, locations, struct, OPM_dict = data
		this_name = 'PPM_runner_{0}'.format(struct)
		structs = [struct]
		OPM_data = {struct : OPM_dict}
	elif data[0] == 'cluster':
		host, options, locations, structs, OPM_data = data
		this_name = 'PPM_runner_cluster'
		num_jobs = int(options['number_of_jobs'])
	else:
		raise(NameError, "ERROR: host must be local or cluster")
	

	success = {}
	for struct in sorted(structs):
		success[struct] = False
		out_filename = struct + '_results.dat'
		outsub_filename = struct + '_datasub1'
		outpar_filename = struct + '_datapar1'

		if (os.path.exists(locations['FSYSPATH']['repoppmresults'] + out_filename) and
		os.path.exists(locations['FSYSPATH']['repoppmdatasub1'] + outsub_filename + '.dat') and
		os.path.exists(locations['FSYSPATH']['repoppmdatapar1'] + outpar_filename + '.dat') and
		os.path.exists(locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb') and
		OPM_data[struct]['OPMcheck'] == True):
			print_log(this_name, "Struct {0} already present, will be skipped".format(struct))
			shutil.copy(locations['FSYSPATH']['repoppmresults'] + out_filename, locations['FSYSPATH']['PPM'])
			shutil.copy(locations['FSYSPATH']['repoppmdatasub1'] + outsub_filename + '.dat', locations['FSYSPATH']['PPM'] + outsub_filename + '.dat')
			shutil.copy(locations['FSYSPATH']['repoppmdatapar1'] + outpar_filename + '.dat', locations['FSYSPATH']['PPM'] + outpar_filename + '.dat')
			shutil.copy(locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb', locations['FSYSPATH']['whole'] + struct + '_opm.pdb')
			success[struct] = True

	if host == 'cluster':
		instr_filename = locations['SYSFILES']['instr_filename'](0)
		instr_file = open(instr_filename, 'w')
		is_done_dict = {}
	for struct in sorted(structs):
		if not success[struct]:
			in_filename = struct + '_input.dat'
			out_filename = struct + '_results.dat'
			outsub_filename = struct + '_datasub1'
			outpar_filename = struct + '_datapar1'
			outpdb_filename = struct + 'out.pdb'

			tmp_folder = locations['FSYSPATH']['PPM'] + 'tmp_' + struct + '/'
			if os.path.exists(tmp_folder):
				shutil.rmtree(tmp_folder)

			# Make a working folder and copy the res.lib file
			os.mkdir(tmp_folder)
			shutil.copy(os.path.dirname(locations['OPT']['ppmexe']) + '/' + 'res.lib', tmp_folder)

			# Generate the full structure from the PDB pdb file and the BIOMATRIX instructions in it (if present)
			OPM_data[struct] = generate_full_structure(locations['FSYSPATH']['PDBpdbs'] + struct + '.pdb', tmp_folder + struct + '.pdb', OPM_data[struct])
			sys.stdout.flush()

			# The runs will always include only one structure
			in_file = open(tmp_folder + in_filename, 'w')
			in_file.write(" 0  in {0}.pdb\n".format(struct))
			in_file.close()


			# If the host is local, it means we are using multiprocess to parallelize and subprocess to run
			if host == 'local':

				# Preparing files for subprocess
				in_file = open(tmp_folder + in_filename, 'r')
				out_file = open(tmp_folder + out_filename, 'w')
				p = subprocess.Popen(locations['OPT']['ppmexe'], stdin=in_file, stdout=out_file, cwd=tmp_folder, stderr=subprocess.PIPE)
				p.wait()
				in_file.close()
				out_file.close()

				success[struct] = check_and_copy(struct, tmp_folder)

				# Cycle closes after one iteration: this is because in this case it's not really a cycle
				toc = time.time()
				print_log(this_name, "Time spent: {0}".format(toc-tic))
				return success[struct], struct, OPM_data[struct]

			# If the host is the cluster, then the routine we are in is not parallel. 
			# We are going to prepare the command to run in the cluster, launch the replicas and we wait for them to finish.
			elif host == 'cluster':
				instr_file.write("cd {0}; {1} < {2} > {3}\n".format(tmp_folder, locations['OPT']['ppmexe'], in_filename, out_filename))
				is_done_dict[struct] = ['', False]
		elif success[struct] and host == 'local':
			toc = time.time()
			print_log(this_name, "Time spent: {0}".format(toc-tic))
			return success[struct], struct, OPM_data[struct]

	if host == 'cluster':
		instr_file.close()
		hostname = options['hostname']
		workdir = os.getcwd() + '/'

		print_log(this_name, "Submitting to the cluster the first {0} jobs".format(num_jobs))
		job_out_filename = {}
		for i in range(num_jobs):
			jobname = 'PPM_' + str(i)
			pbs_filename = locations['SYSFILES']['PPMpbs'](jobname)
			jobname, walltime_hrs, job_out_filename[jobname] = create_pbs(hostname, locations, instr_filename, pbs_filename, jobname)
			subprocess.call(["ssh", "{0}".format(hostname), "qsub", "{0}".format(locations['SYSFILES']['PPMpbs'](jobname))], cwd=workdir)
			time.sleep(10)

		time.sleep(60)
		while False in [x[1] for x in list(is_done_dict.values())]:
			print_log(this_name, "is_done_dict: {0}".format(is_done_dict.items()))
			for struct in sorted(list(is_done_dict.keys())):
				if not is_done_dict[struct][1]:
					tmp_folder = locations['FSYSPATH']['PPM'] + 'tmp_' + struct + '/'
					if os.path.exists(tmp_folder + struct + 'out.pdb'):
						is_done_dict[struct][1] = True

			print_log(this_name, "After check, is_done_dict: {0}".format(is_done_dict.items()))
			if not False in [x[1] for x in list(is_done_dict.values())]:
				print_log(this_name, "All alignments completed")
				break

			print_log(this_name, "Something is still missing")
			job_out_filename = job_control(locations, walltime_hrs, hostname, job_out_filename, is_done_dict, num_jobs, 'PPM')
			time.sleep(60)

		PPM_success = set()
		for struct in sorted(structs):
			tmp_folder = locations['FSYSPATH']['PPM'] + 'tmp_' + struct + '/'
			# If it already succeeded it's because it was not processed (and thus it is already copied)
			if success[struct]:
				PPM_success.add(struct)
			# If it was not processed, then check_and_copy and see if it is a success
			else:
				success[struct] = check_and_copy(struct, tmp_folder)
				if success[struct]:
					PPM_success.add(struct)
				
		toc = time.time()
		print_log(this_name, "Time spent: {0}".format(toc-tic))
		sys.stdout.flush()	
		return PPM_success, structs, OPM_data


def PPM_extractor(locations, PPM_set, run_in_PPM, OPM_data):
	this_name = 'PPM_extractor'
	tic = time.time()
	print_log(this_name, "Set of PPM runs: {0}".format(''.join([x+' ' for x in sorted(list(PPM_set))])))
	exclusions_filename = locations['SYSFILES']['excludedchains']
	if os.path.exists(exclusions_filename):
		exclusions_file = open(exclusions_filename, 'a')
	else:
		exclusions_file = open(exclusions_filename, 'w')

	for struct in sorted(list(OPM_data.keys())):
		if struct in PPM_set:
			print_log(this_name, "Struct {0}".format(struct))
			OPM_data[struct]['tmchains'] = set()
			datasub1_filename = locations['FSYSPATH']['PPM'] + struct + '_datasub1.dat'
			datasub1_file = open(datasub1_filename, 'r')
			text = datasub1_file.read().split('\n')
			for line in text:
				if not line:
					continue
				sc_fields = line.split(';')
				chain = sc_fields[1]
				OPM_data[struct]['tmchains'].add(chain)
				OPM_data[struct][chain] = find_segments(sc_fields[3])
		elif struct in run_in_PPM and struct not in PPM_set:
			exclusions_file.write(struct + '\t\t' + 'Inconsistencies in the PPM run' + '\n')
			del OPM_data[struct]
		elif struct not in run_in_PPM:
			pass
	exclusions_file.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return OPM_data

	
def check_tmdoms_and_transfer(struct_filename, chain_filename, chain, segments, resolution):
	this_name = 'check_tmdoms_and_transfer'
	struct_file = open(struct_filename, 'r')
	text = struct_file.read().split('\n')
	struct_file.close()

	chain_file = open(chain_filename, 'w')
	limit_inf, limit_sup = None, None
	lzR = {}
	hzR = {}
	residmax = 'X'
	residmin = 'X'
	for line in text:
		if line[0:4] == 'ATOM' and line[21] == chain:
			chain_file.write(line + '\n')
			if line[13:16].strip() == 'CA':
				resid = int(line[22:26].strip())
				residmax = resid
				if residmin == 'X':
					residmin = resid
				if not resid in lzR or lzR[resid] > float(line[46:54]):
					lzR[resid] = float(line[46:54])
				if not resid in hzR or hzR[resid] < float(line[46:54]):
					hzR[resid] = float(line[46:54])
		if (line[0:6]=='HETATM' or line[0:4]=='ATOM') and line[17:20]=='DUM':
			if not (limit_inf and limit_sup):
				if line[12:16].strip()=='N':
					limit_inf = float(line[46:54]) + resolution/2.0
				elif line[12:16].strip()=='O':
					limit_sup = float(line[46:54]) - resolution/2.0
	chain_file.close()

	l = list(iterchain.from_iterable(segments))
	prev = -1
	loops = []
	for n in range(1,len(l)):
		if prev == -1:
			prev = int(l[n])
		else:
			loops.append((prev, int(l[n])))
			prev = -1


	print_log(this_name, "{0}\t{1}\n".format(struct_filename, chain_filename))
	# Calculates start of initial segment
	rev_goes_outside_at = 'X'
	for i in reversed(range(residmin, int(segments[0][0].strip())+1)):
		if i in lzR:
			print_log(this_name, "{0}".format(i))
			if lzR[i] < limit_inf or hzR[i] > limit_sup:
				rev_goes_outside_at = i
				print_log(this_name, "goes outside")
				break
	if rev_goes_outside_at == 'X':
		rev_goes_outside_at = residmin
	tmdoms = [(rev_goes_outside_at, 'X')]
	# Calculates end of nth segment and start of (n+1)th segment
	for n, loop in enumerate(loops):
		goes_outside_at = 'X'
		goes_back_inside_at = 'X'
		init, end = loop
		print_log(this_name, "loop {0} {1}".format(init, end))
		for i in range(init, end+1):
			if i in lzR:
				print_log(this_name, "{0}".format(i))
				if goes_outside_at == 'X' and (lzR[i] < limit_inf or hzR[i] > limit_sup):
					print_log(this_name, "goes outside")
					goes_outside_at = i
				if goes_outside_at != 'X' and goes_back_inside_at == 'X' and lzR[i] > limit_inf and hzR[i] < limit_sup:
					print_log(this_name, "goes inside")
					goes_back_inside_at = i
				if goes_back_inside_at != 'X' and (lzR[i] < limit_inf or hzR[i] > limit_sup):
					print_log(this_name, "last 'goes inside' was fake, here it goes outside again")
					goes_back_inside_at = 'X'
		if goes_outside_at != 'X':
			tmdoms[-1] = (tmdoms[-1][0], goes_outside_at)
			if goes_back_inside_at == 'X':
				goes_back_inside_at = end
			tmdoms.append((goes_back_inside_at, 'X'))
		else:
			print_log(this_name, "fake loop")
	goes_outside_at = 'X'
	# Calculates end of final segment
	for i in range(int(segments[-1][1].strip()), residmax+1):
		if i in lzR:
			print_log(this_name, "{0}".format(i))
			if lzR[i] < limit_inf or hzR[i] > limit_sup:
				goes_outside_at = i
				print_log(this_name, "goes outside")
				break
	if goes_outside_at == 'X':
		goes_outside_at = residmax+1
	tmdoms[-1] = (tmdoms[-1][0], goes_outside_at)

	print_log(this_name, "Chain {0}".format(chain))
	print_log(this_name, "TM domains: {0}".format(tmdoms))
	return tmdoms					


def check_chainwise(locations, filters, OPM_data):
	this_name = 'check_chainwise'
	tic = time.time()
	instructions = {}
	instructions_filename = locations['SYSFILES']['H_topologytype']
	instructions_file = open(instructions_filename, 'w')
	instructions_straln = {}
	instructions_straln_filename = locations['SYSFILES']['H_instraln']
	instructions_straln_file = open(instructions_straln_filename, 'w')
	exclusions_filename = locations['SYSFILES']['excludedchains']
	exclusions_straln_filename = locations['SYSFILES']['excludedstralnchains']
	equivalence_filename = locations['SYSFILES']['equivalences']
	already_excluded = []
	if os.path.exists(exclusions_filename):
		exclusions_file = open(exclusions_filename, 'r')
		text = exclusions_file.read().split('\n')
		for line in text:
			if line and line[:6]:
				already_excluded.append(line[:6])
		exclusions_file.close()
		exclusions_file = open(exclusions_filename, 'a')
	else:
		exclusions_file = open(exclusions_filename, 'w')
	already_excluded_straln = []
	if os.path.exists(exclusions_straln_filename):
		exclusions_straln_file = open(exclusions_straln_filename, 'r')
		text = exclusions_straln_file.read().split('\n')
		for line in text:
			if line and line[:6]:
				already_excluded_straln.append(line[:6])
		exclusions_straln_file.close()
		exclusions_straln_file = open(exclusions_straln_filename, 'a')
	else:
		exclusions_straln_file = open(exclusions_straln_filename, 'w')
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
		equivalence_file = open(equivalence_filename, 'a')
	else:
		equivalence_file = open(equivalence_filename, 'w')
		

	for struct in sorted(list(OPM_data.keys())):
		for chain in sorted(OPM_data[struct]['tmchains']):
			chain_filename = locations['FSYSPATH']['chains'] + struct + '_' + chain + '.pdb'
			if not os.path.exists(chain_filename):
				struct_filename = locations['FSYSPATH']['whole'] + struct + '_opm.pdb'
				segments = OPM_data[struct][chain]['segments']
				resolution = float(OPM_data[struct]['FROM_PDB']['RESOLUTION'])
				OPM_data[struct][chain]['tmdoms'] = check_tmdoms_and_transfer(struct_filename, chain_filename, chain, segments, resolution)
				OPM_data[struct][chain]['ntm'] = str(len(OPM_data[struct][chain]['tmdoms']))
			s_type = OPM_data[struct]['class']
			n_pdbtm = OPM_data[struct][chain]['ntm']


	for struct in sorted(list(OPM_data.keys())):
		print_log(this_name, "PDB code: {0}".format(struct))
		for chain in sorted(OPM_data[struct]['tmchains']):
			print_log(this_name,"Struct {0}_{1}".format(struct, chain))
			exclude_chain = []
			exclude_chain_straln = []
			if not exclude_chain:
				s_type = OPM_data[struct]['class']
				n_pdbtm = OPM_data[struct][chain]['ntm']
				if 'equivalent_chains' in OPM_data[struct] and chain in OPM_data[struct]['equivalent_chains'] and OPM_data[struct]['equivalent_chains'][chain] in OPM_data[struct]['tmchains'] and OPM_data[struct][OPM_data[struct]['equivalent_chains'][chain]]['ntm'] == OPM_data[struct][chain]['ntm']:
					if struct + '_' + chain not in equivalence:
						equivalence_file.write('{0}\t{1}\n'.format(struct + '_' + chain, struct + '_' + OPM_data[struct]['equivalent_chains'][chain]))
					exclude_chain_straln.append("Chain {0} is equivalent to chain {1}\n".format(chain, OPM_data[struct]['equivalent_chains'][chain]))
				if not exclude_chain_straln:
					instructions_straln_file.write("{0}\t{1}\t{2}\n".format(s_type, n_pdbtm, struct+'_'+chain))
					if struct not in instructions_straln:
						instructions_straln[struct] = {}
					instructions_straln[struct][chain] = (s_type, n_pdbtm, OPM_data[struct]['FASTA'][chain])
				else:
					if struct + '_' + chain not in already_excluded_straln:
						exclusions_straln_file.write(struct + '_' + chain + '\t\t' + exclude_chain_straln[0] + '\n')
						for nl in range(1, len(exclude_chain_straln)):
							exclusions_straln_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_chain_straln[nl] + '\n')
				instructions_file.write("{0}\t{1}\t{2}\n".format(s_type, n_pdbtm, struct+'_'+chain))
				if struct not in instructions:
					instructions[struct] = {}
				instructions[struct][chain] = (s_type, n_pdbtm, OPM_data[struct]['FASTA'][chain])
			else:
				if struct + '_' + chain not in already_excluded:
					exclusions_file.write(struct + '_' + chain + '\t\t' + exclude_chain[0] + '\n')
					for nl in range(1, len(exclude_chain)):
						exclusions_file.write(' '*len(struct) + ' ' + ' '*len(chain) + '\t\t' + exclude_chain[nl] + '\n')
	equivalence_file.close()
	exclusions_file.close()
	exclusions_straln_file.close()
	instructions_file.close()
	instructions_straln_file.close()
	shutil.copy(instructions_filename, locations['FSYSPATH']['repository'])
	shutil.copy(instructions_straln_filename, locations['FSYSPATH']['repository'])
	shutil.copy(equivalence_filename, locations['FSYSPATH']['repository'])
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))

	pymolpath = options['pymol_path'] 
	for pdbname in sorted(list(OPM_data.keys())):
		os.system("cp create_whole_png.cmd {0}; cd {0}; cp {1} .; sed -i 's/STRUCT/{3}/g' create_whole_png.cmd; {5} -ucq create_whole_png.cmd; mv {2}.png {4}.png".format(locations['FSYSPATH']['wholegif'], locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb', 'f'+str(1).zfill(4), pdbname + '_opm.pdb', pdbname, pymolpath))
		for chain in sorted(OPM_data[pdbname]['tmchains']):
			os.system("cp create_chain_png.cmd {0}; cd {0}; cp {1} .; sed -i 's/STRUCT/{3}/g' create_chain_png.cmd; sed -i 's/XXX/{5}/g' create_chain_png.cmd; {6} -ucq create_chain_png.cmd; mv {2}.png {4}.png".format(locations['FSYSPATH']['chainsgif'], locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb', 'f'+str(1).zfill(4), pdbname + '_opm.pdb', pdbname+'_'+chain, chain, pymolpath))
	return OPM_data, instructions


def structure_sorter(locations, instructions):
	this_name = 'structure_sorter'
	tic = time.time()
	ttd = {'alpha' : 'a', 'beta' : 'b'}

	for tt in 'alpha', 'beta':
		for i in sorted(os.listdir(locations['FSYSPATH'][tt])):
			if os.path.exists(locations['FSYSPATH'][tt] + i + '/' + locations['TREE']['seq']):
				shutil.rmtree(locations['FSYSPATH'][tt] + i + '/' + locations['TREE']['seq'])
			if os.path.exists(locations['FSYSPATH'][tt] + i + '/' + locations['TREE']['str']):
				shutil.rmtree(locations['FSYSPATH'][tt] + i + '/' + locations['TREE']['str'])
	for struct in sorted(instructions):
		for chain in sorted(instructions[struct]):
			tt = instructions[struct][chain][0]
			ntm = instructions[struct][chain][1]
			destination_dir = locations['FSYSPATH'][tt] + str(ntm) + '/'
			if not os.path.exists(destination_dir):
				os.mkdir(destination_dir)
			if not os.path.exists(destination_dir + locations['TREE']['str']):
				os.mkdir(destination_dir + locations['TREE']['str'])
			if not os.path.exists(destination_dir + locations['TREE']['seq']):
				os.mkdir(destination_dir + locations['TREE']['seq'])
			pdb_filename = struct + '_' + chain + '.pdb'
			shutil.copy(locations['FSYSPATH']['chains'] + pdb_filename, destination_dir + locations['TREE']['str'])
			fasta_filename = locations['FSYSPATH'][tt] + str(ntm) + '/' + locations['TREE']['seq'] + struct + '_' + chain + '.fa'
			fasta_file = open(fasta_filename, 'w')
			fasta_file.write(">{0}_{1}\n{2}\n".format(struct, chain, instructions[struct][chain][2]))
			fasta_file.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))


def gif_creator(locations, instructions):
	this_name = 'gif_creator'
	tic = time.time()
	script_folder = os.getcwd() + '/'
	fnull = open(os.devnull, 'w')
	
	for struct in instructions:
		if os.path.exists(locations['FSYSPATH']['repochainsgif'] + '{0}_gif'.format(struct)):
			shutil.copy(locations['FSYSPATH']['repochainsgif'] + '{0}_gif'.format(struct), locations['FSYSPATH']['wholegif'])
			continue
		os.mkdir(locations['FSYSPATH']['whole'] + '{0}_gif'.format(struct))
		os.chdir(locations['FSYSPATH']['whole'] + '{0}_gif'.format(struct))
		shutil.copy(locations['FSYSPATH']['whole'] + struct + '_opm.pdb', './')
		shutil.copy(script_folder + 'create_whole_gif.cmd', locations['FSYSPATH']['whole'] + '{0}_gif'.format(struct))
		bashcommand = "sed -i 's/STRUCT/{0}/g' create_whole_gif.cmd".format(struct + '_opm.pdb')
		os.system(bashcommand)
		p = subprocess.Popen(['pymol', '-ucq', 'create_whole_gif.cmd'], stdout=fnull, stderr=fnull, cwd=(locations['FSYSPATH']['whole']+'{0}_gif'.format(struct)))
		p.wait()
		nframe = 1
		bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
		os.system(bashcommand)
		for nframe in range(2,51):
			bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
			os.system(bashcommand)
			bashcommand = "convert {0}.png {1}.gif".format('f'+str(nframe).zfill(4), 'f'+str(102-nframe).zfill(4))
			os.system(bashcommand)
		nframe = 51
		bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
		os.system(bashcommand)
		bashcommand = "convert -loop 0 -delay 3 -type palette f*.gif {0}.gif".format(struct)
		os.system(bashcommand)
		shutil.copy('{0}.gif'.format(struct), locations['FSYSPATH']['wholegif'])
		shutil.copy('{0}.gif'.format(struct), locations['FSYSPATH']['repowholegif'])
		os.chdir(script_folder)
		shutil.rmtree(locations['FSYSPATH']['whole'] + '{0}_gif'.format(struct))

	for struct in instructions:
		for chain in instructions[struct]:
			if os.path.exists(locations['FSYSPATH']['repochainsgif'] + '{0}_{1}.gif'.format(struct, chain)):
				shutil.copy(locations['FSYSPATH']['repochainsgif'] + '{0}_{1}.gif'.format(struct, chain), locations['FSYSPATH']['chainsgif'])
				continue
			os.mkdir(locations['FSYSPATH']['chains'] + '{0}_{1}_gif'.format(struct, chain))
			os.chdir(locations['FSYSPATH']['chains'] + '{0}_{1}_gif'.format(struct, chain))
			shutil.copy(locations['FSYSPATH']['whole'] + struct + '_opm.pdb', './')  # It's ok, we need the whole structure
			shutil.copy(script_folder + 'create_chain_gif.cmd', locations['FSYSPATH']['chains'] + '{0}_{1}_gif'.format(struct, chain))
			bashcommand = "sed -i 's/STRUCT/{0}/g' create_chain_gif.cmd".format(struct + '_opm.pdb')
			os.system(bashcommand)
			bashcommand = "sed -i 's/XXX/{0}/g' create_chain_gif.cmd".format(chain)
			os.system(bashcommand)
			p = subprocess.Popen(['pymol', '-ucq', 'create_chain_gif.cmd'], stdout=fnull, stderr=fnull, cwd=locations['FSYSPATH']['chains']+'{0}_{1}_gif'.format(struct, chain))
			p.wait()
			nframe = 1
			bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
			os.system(bashcommand)
			for nframe in range(2,51):
				bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
				os.system(bashcommand)
				bashcommand = "convert {0}.png {1}.gif".format('f'+str(nframe).zfill(4), 'f'+str(102-nframe).zfill(4))
				os.system(bashcommand)
			nframe = 51
			bashcommand = "convert {0}.png {0}.gif".format('f'+str(nframe).zfill(4))
			os.system(bashcommand)
			bashcommand = "convert -loop 0 -delay 3 -type palette f*.gif {0}_{1}.gif".format(struct, chain)
			os.system(bashcommand)
			shutil.copy('{0}_{1}.gif'.format(struct, chain), locations['FSYSPATH']['chainsgif'])
			shutil.copy('{0}_{1}.gif'.format(struct, chain), locations['FSYSPATH']['repochainsgif'])
			os.chdir(script_folder)
			shutil.rmtree(locations['FSYSPATH']['chains'] + '{0}_{1}_gif'.format(struct, chain))
	fnull.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))


# Library function
def generate_chain_pdb_files(options, filters, locations, OPM_data, run_in_PPM):
	# Hardcoded variables
	this_name = 'genclib'
	indent = " "*len(header(this_name))
	version = 3.1

	# Checks
	for path_name in [x[1] for n, x in enumerate(locations['FSYSPATH'].items()) if n > 0]:
		if not os.path.exists(path_name):
			logmsg = header(this_name) + "ERROR: The directory path {0} does not exist. Please generate the file system first.".format(path_name)
			write_log(this_name, logmsg)	
			raise NameError(logmsg)

	# Whole structures checker
	OPM_data, run_in_PPM = check_whole_structures(locations, filters, OPM_data, run_in_PPM)
	write_data(locations['FSYSPATH']['main'] + '.run_in_PPM.dat', run_in_PPM)

	# PPM runner
	pool_outputs = []
	if run_in_PPM:
		if 'hostname' in options and options['hostname']:
			host = 'cluster'
			data = host, options, locations, run_in_PPM, OPM_data
			PPM_success, placeholder, OPM_data = PPM_runner(data)
		else:
			host = 'local'
			data = []
			np = int(options['number_of_procs'])
			for struct in sorted(run_in_PPM):
				data.append((host, locations, struct, OPM_data[struct]))

			pool = multiprocessing.Pool(processes=np)
			pool_outputs = pool.map(PPM_runner, data)
			pool.close()
			pool.join()

			PPM_success = set()
			for s in pool_outputs:
				if not s:
					continue
				success, struct, OPM_dict = s
				print(success, struct)
				print(list(OPM_dict.keys()))
				if success:
					PPM_success.add(struct)
					OPM_data[struct] = OPM_dict
#				else:
#					if struct in OPM_data:
#						del OPM_data[struct]
			sys.stdout.flush()

		OPM_data_success_PPM = PPM_extractor(locations, PPM_success, run_in_PPM, OPM_data)

	OPM_data, instructions = check_chainwise(locations, filters, OPM_data)

	OPM_archive(locations, OPM_data, archive_filename=locations['FSYSPATH']['chains']+'OPM_archive.dat')

	# Filesystem sorting
	structure_sorter(locations, instructions)

#	gif_creator(locations, instructions)

	write_data(locations['SYSFILES']['opmdata'], OPM_data)

	return OPM_data
