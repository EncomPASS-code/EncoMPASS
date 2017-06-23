# Name: support.py
# Language: python3
# Libraries: datetime
# Description: Support functions used by modules
# Author: Edoardo Sarti
# Date: Aug 10 2016

import datetime
import shutil
import os
import sys
import urllib.request
import gzip
import time
import codecs
import copy
import re
import time
import subprocess
import operator as op
from functools import reduce
from bs4 import BeautifulSoup
from collections import OrderedDict
from colors import *
from Bio.PDB import *
import warnings
from Bio.PDB.PDBParser import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

# Support functions
def write_log(name, text):
	log_filename = name + '.log'
	log_file = open(log_filename, 'w')
	log_file.write(text)
	log_file.close()


def header(name):
	return "[" + name + " " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "] "


def raise_error(name, text):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0] + '\n'
	for line in lines[1:]:
		logmsg += indent + line + '\n'
	write_log(name, logmsg)
	raise NameError(logmsg)


def print_log(name, text, quiet=False):
	indent = " "*len(header(name))
	lines = text.split('\n')
	logmsg = header(name) + lines[0]
	for line in lines[1:]:
		logmsg += '\n' + indent + line
	if not quiet:
		print(logmsg)
	return logmsg + '\n'


def archive_old_file(locations, filenames):
	if type(filenames) == str:
		filenames = [filenames]
	for filename in filenames:
		filename_noext = os.path.basename(os.path.splitext(filename)[0])
		extension = os.path.splitext(filename)[1]
		one_hundred_files = False
		for i in range (0,100):
			new_filename = locations['FSYSPATH']['old'] + filename_noext + '_' + str(i) + extension
			if not os.path.exists(new_filename):
				shutil.move(filename, new_filename)
				break
			if i == 99:
				one_hundred_files = True
		if one_hundred_files:
			for i in range(1, 100):
				old_filename = locations['FSYSPATH']['old'] + filename_noext + '_' + str(i) + extension
				new_filename = locations['FSYSPATH']['old'] + filename_noext + '_' + str(i-1) + extension
				shutil.move(old_filename, new_filename)
			new_filename = locations['FSYSPATH']['old'] + filename_noext + '_' + str(99) + extension
			shutil.move(filename, new_filename)

def repo_inspector(repo_filename):
	this_name = 'repo_inspector'
	repo_info = {}
	chain_1 = ""
	if os.path.exists(repo_filename):
		repo_file = open(repo_filename, 'r')
		text = repo_file.read().split('\n')
		repo_file.close()
		record = False
		for line in text:
			fields = line.split()
			if line and fields[0] == 'BEGIN':
				chain_1 = fields[2]
				chain_2 = fields[4]
				record = True
				recorded_text = ""
				continue
			if line and fields[0] == 'END':
				if chain_2 in repo_info:
					print_log(this_name, "WARNING: multiple occurrences of couple {0} {1} in file {2}".format(chain_1, chain_2, repo_filename))
				repo_info[chain_2] = recorded_text
				record = False
				continue
			if record:
				recorded_text += line + '\n'
	if repo_info and not chain_1:
		raise_error('repo_inspector', "ERROR: chain_1 is empty but the repository dictionary is not")
	return chain_1, repo_info


def write_on_repo(repo_filename, textdict, append=False):
	if append:
		repo_file = open(repo_filename, 'a')
	else:
		repo_file = open(repo_filename, 'w')
	for chain_1 in sorted(list(textdict.keys())):
		for chain_2 in sorted(list(textdict[chain_1].keys())):
			repo_file.write("BEGIN\t\tCHAIN_1: " + chain_1  + "\tCHAIN_2: " + chain_2 + "\n" + textdict[chain_1][chain_2] + "\nEND\n\n\n")
	repo_file.close()


def complete_info_from_repo(repo_2dict, repo_location, suffix):
	this_name = 'complete_info_from_repo'
	lensuf = len(suffix)
	for filename in os.listdir(repo_location):
		if filename[-lensuf:] == suffix:
			print_log(this_name, "{0}".format(repo_location + filename))
			chain_1, temp_repo_dict = repo_inspector(repo_location + filename)
			if chain_1 not in repo_2dict:
				repo_2dict[chain_1] = temp_repo_dict
	return repo_2dict


def simple_download(url, local_filename):
	try:
		response = urllib.request.urlopen(url)
	except:
		return False

	local_file = open(local_filename, 'wb')
	shutil.copyfileobj(response, local_file)
	local_file.close()
	return True


def reiterated_simple_download(url, local_filename):
	this_name = 'reiterated_simple_download'
	max_r = 3
	error_found = False
	downloaded = False
	reiterations = 0
	while (error_found or not downloaded):
		if reiterations >= max_r:
			return False
		if reiterations > 0:
			time.sleep(reiterations**2)
			print_log(this_name, "Try number: {0}".format(reiterations))
		downloaded = simple_download(url, local_filename)
		if downloaded:
			with codecs.open(local_filename, "r",encoding='utf-8', errors='ignore') as local_file:
				for line in local_file:
					if "Error: Could not connect to database. Please try again later." in line:
						error_found = True
		reiterations += 1
	return True


def gzip_download(url, local_filename):
	try:
		response = urllib.request.urlopen(url)
	except:
		return False

	with gzip.GzipFile(fileobj=response) as uncompressed, open(local_filename, 'wb') as local_file:
		shutil.copyfileobj(uncompressed, local_file)
		local_file.close()
	return True


def reiterated_gzip_download(url, local_filename):
	this_name = 'reiterated_gzip_download'
	max_r = 3
	error_found = False
	downloaded = False
	reiterations = 0
	while (error_found or not downloaded):
		if reiterations >= max_r:
			return False
		if reiterations > 0:
			time.sleep(reiterations**2)
			print_log(this_name, "try number: {0}".format(reiterations))
		downloaded = gzip_download(url, local_filename)
		if downloaded:
			with codecs.open(local_filename, "r",encoding='utf-8', errors='ignore') as local_file:
				for line in local_file:
					if "Error: Could not connect to database. Please try again later." in line:
						error_found = True
		reiterations += 1
	return True


def find_segments(line):
	result_dict = {}
	in_par = False
	segments = []
	for c in line:
		if c == '(':
			in_par = True
			tc = []
		if c == ')' and in_par:
			in_par = False
			segments.append(''.join(tc))
		if in_par:
			tc.append(c)
	result_dict['nseg'] = str(len(segments))
	result_dict['segments'] = []
	for s in segments:
		initial = s.split('-')[0].replace('(','').strip()
		final = s.split('-')[1]
		result_dict['segments'].append((initial, final))
	return result_dict

def OPM_archive(locations, OPM_data, archive_filename=None):
	this_name = 'OPM_archive'
	tic = time.time()
	attributes = ['title', 'class', 'family', 'superfamily', 'TCDBfamily', 'TCDBsuperfamily', 'tmchains']
	if not archive_filename:
		archive_filename = locations['SYSFILES']['OPMarchive']
	archive_file = open(archive_filename, 'w')
	for struct in sorted(list(OPM_data.keys())):
		print_log(this_name, "Struct {0}".format(struct))
		if not OPM_data[struct]['is_representative']:
			continue
		line = "{0:10} {1:20} ".format(struct, "representative")
		if 'tmchains' in OPM_data[struct] and OPM_data[struct]['tmchains']:
			nch = len(OPM_data[struct]['tmchains'])
		else:
			nch = 0
		line += "{0:<10} ".format(nch)
		for key in attributes:
			if key == 'tmchains':
				continue
			datum = ""
			if key in OPM_data[struct] and OPM_data[struct][key]:
				if type(OPM_data[struct][key]) == list or type(OPM_data[struct][key]) == tuple:
					datum = "".join([x+' ' for x in OPM_data[struct][key]])
					datum = '"'+datum.strip()+'"'
				else:
					datum = str(OPM_data[struct][key])
			line += "{0:50} ".format(datum)
		archive_file.write(line.rstrip()+'\n')
		if 'tmchains' in OPM_data[struct] and OPM_data[struct]['tmchains']:
			for chain in sorted(OPM_data[struct]['tmchains']):
				line = "{0:10} {1:20} ".format(struct, "repr_chain")
				line += "{0:<10} ".format(chain)
				if 'ntm' in OPM_data[struct][chain]:
					line += "{0:<10} ".format(int(OPM_data[struct][chain]['ntm']))
					subline = "{0} ".format("".join(['({0},{1}) '.format(x[0], x[1]) for x in OPM_data[struct][chain]['tmdoms']]))
				else:
					line += "{0:<10} ".format(int(OPM_data[struct][chain]['nseg']))
					subline = "{0} ".format("".join(['({0},{1}) '.format(x[0], x[1]) for x in OPM_data[struct][chain]['segments']]))
				line += '['+subline.strip()+']'
				archive_file.write(line.rstrip()+'\n')
		if 'related' in OPM_data[struct] and OPM_data[struct]['related']:
			for relstruct in sorted(OPM_data[struct]['related']):
				line = "{0:10} {1:20} {2:10} ".format(relstruct, "related", struct)
				archive_file.write(line.rstrip()+'\n')
	archive_file.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))


def write_data(data_filename, opm_data):
	this_name = 'write_data'
	tic = time.time()
	def go_through(data_file, spath, s):
		if type(s) == dict or type(s) == OrderedDict:
			for key in list(s.keys()):
				if type(s[key]) == list:
					data_file.write('list\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				elif type(s[key]) == tuple:
					data_file.write('tuple\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				elif type(s[key]) == set:
					data_file.write('set\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				elif type(s[key]) == frozenset:
					data_file.write('fset\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				elif type(s[key]) == dict:
					data_file.write('dict\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				elif type(s[key]) == OrderedDict:
					data_file.write('odict\t{0}\n'.format(key))
					spath = go_through(data_file, spath, s[key])
				else:
					data_file.write('entr\t{0}\t"{1}"\n'.format(key, s[key]))
			if type(s) == dict:
				data_file.write('\\dict\n')
			elif type(s) == OrderedDict:
				data_file.write('\\odict\n')
		elif type(s) == list or type(s) == tuple or type(s) == set or type(s) == frozenset:
			for item in s:
				if type(item) == list:
					data_file.write('list\n'.format())
					spath = go_through(data_file, spath, item)
				elif type(item) == tuple:
					data_file.write('tuple\n'.format())
					spath = go_through(data_file, spath, item)
				elif type(item) == set:
					data_file.write('set\n'.format())
					spath = go_through(data_file, spath, item)
				elif type(item) == frozenset:
					data_file.write('fset\n'.format())
					spath = go_through(data_file, spath, item)
				elif type(item) == dict:
					data_file.write('dict\n'.format())
					spath = go_through(data_file, spath, item)
				elif type(item) == OrderedDict:
					data_file.write('odict\n'.format())
					spath = go_through(data_file, spath, item)
				else:
					data_file.write('entr\t"{0}"\n'.format(item))
			if type(s) == list:
				data_file.write('\\list\n')
			elif type(s) == tuple:
				data_file.write('\\tuple\n')
			elif type(s) == set:
				data_file.write('\\set\n')
			elif type(s) == frozenset:
				data_file.write('\\fset\n')
		return spath

	print_log(this_name, "Writing on {0}".format(data_filename))
	data_file = open(data_filename, 'w')

	spath = []
	s = copy.deepcopy(opm_data)
	spath.append(s)

	# First step
	if type(s) == list:
		data_file.write('list\n')
	elif type(s) == tuple:
		data_file.write('tuple\n')
	elif type(s) == set:
		data_file.write('set\n')
	elif type(s) == frozenset:
		data_file.write('fset\n')
	elif type(s) == dict:
		data_file.write('dict\n')
	elif type(s) == OrderedDict:
		data_file.write('odict\n')

	ends = go_through(data_file, spath, s)
	data_file.close()
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))

def splice_line(line):
	openframe = False
	fields = []
	field = ""
	for char in line:
		if char == "\"" or char == "'":
			openframe = not(openframe)
		elif char == " " or char == "\t" or char == "\n":
			if openframe:
				field += char
			else:
				if field:
					fields.append(field)
				field = ""
		else:
			field += char
	if field:
		fields.append(field)
	return fields
			

def read_data(data_filename):
	this_name = 'read_data'
	tic = time.time()
	print_log(this_name, "Reading from {0}".format(data_filename))

	data_file = open(data_filename, 'r')
	text = data_file.read().split('\n')
	data_file.close()

	buildlist = []
	for line in text:
		fields = splice_line(line)
		if not fields:
			continue
#		print(line)
		if fields[0] == 'entr':
			if type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict:
				if fields[2] == 'True':
					entr = True
				elif fields[2] == 'False':
					entr = False
				else:
					entr = fields[2]
				buildlist[-1][fields[1]] = entr
			else:
				if fields[1] == 'True':
					entr = True
				elif fields[1] == 'False':
					entr = False
				else:
					entr = fields[1]
				if type(buildlist[-1]) == list:
					buildlist[-1].append(fields[1])
				elif type(buildlist[-1]) == tuple:
					tmptpl = buildlist[-1]
					buildlist[-1] = (tuple([x for x in tmptpl] + [fields[1]]))
				elif type(buildlist[-1]) == set:
					buildlist[-1].add(fields[1])
				elif type(buildlist[-1]) == frozenset:
					buildlist[-1] = frozenset([x for x in buildlist[-1]] + [fields[1]])
		elif fields[0] == 'list':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append([])
		elif fields[0] == 'tuple':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append(())
		elif fields[0] == 'set':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append(set())
		elif fields[0] == 'fset':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append(frozenset())
		elif fields[0] == 'dict':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append({})
		elif fields[0] == 'odict':
			if len(buildlist) > 0 and (type(buildlist[-1]) == dict or type(buildlist[-1]) == OrderedDict):
				buildlist[-1][fields[1]] = None
			buildlist.append(OrderedDict())
		elif fields[0][0] == '\\':
			if len(buildlist) > 1:
				if type(buildlist[-2]) == dict or type(buildlist[-2]) == OrderedDict:
					for k, v in buildlist[-2].items():
						if v == None:
							buildlist[-2][k] = buildlist[-1]
							del buildlist[-1]
				elif type(buildlist[-2]) == list:
					buildlist[-2].append(buildlist[-1])
					del buildlist[-1]
				elif type(buildlist[-2]) == tuple:
					tmptpl = buildlist[-2]
					buildlist[-2] = tuple([x for x in tmptpl] + [buildlist[-1]])
					del buildlist[-1]
				elif type(buildlist[-2]) == set:
					buildlist[-2].add(buildlist[-1])
					del buildlist[-1]
				elif type(buildlist[-2]) == frozenset:
					tmptpl = buildlist[-2]
					buildlist[-2] = frozenset([x for x in tmptpl] + [buildlist[-1]])
	toc = time.time()
	print_log(this_name, "Time spent: {0}".format(toc-tic))
	return buildlist[0]


def write_roman(num):

	roman = OrderedDict()
	roman[1000] = "M"
	roman[900] = "CM"
	roman[500] = "D"
	roman[400] = "CD"
	roman[100] = "C"
	roman[90] = "XC"
	roman[50] = "L"
	roman[40] = "XL"
	roman[10] = "X"
	roman[9] = "IX"
	roman[5] = "V"
	roman[4] = "IV"
	roman[1] = "I"

	def roman_num(num):
		for r in roman.keys():
			x, y = divmod(num, r)
			yield roman[r] * x
			num -= (r * x)
			if num > 0:
				roman_num(num)
			else:
				break

	return "".join([a for a in roman_num(num)])

def from3to1(resname):
	f3t1 = {'ALA' : 'A',
	        'ARG' : 'R',
	        'ASN' : 'N',
	        'ASP' : 'D',
	        'CYS' : 'C',
	        'GLN' : 'Q',
	        'GLU' : 'E',
	        'GLY' : 'G',
	        'HIS' : 'H',
	        'ILE' : 'I',
	        'LEU' : 'L',
	        'LYS' : 'K',
	        'MET' : 'M',
	        'PHE' : 'F',
	        'PRO' : 'P',
	        'SER' : 'S',
	        'THR' : 'T',
	        'TRP' : 'W',
	        'TYR' : 'Y',
	        'VAL' : 'V'}

	if resname in list(f3t1.keys()):
		return f3t1[resname]
	else:
		return '0'

# PDB parser
def PDB_parser(location, struct):
	this_name = 'PDB_parser'
	struct_filename = location + struct + '.pdb'
	struct_file = open(struct_filename, 'r')
	text = struct_file.read().split("\n")
	struct_file.close()

	PDB_dict = {}
	res_ids = {}
	b_factor = {}
	b_norm = {}
	biomolecule_occurrences = 0
	tech_list = ['NMR', 'X-RAY', 'THEORETICAL', 'ELECTRON']
	print_log(this_name, "Struct {0}".format(struct))
	for line in text:
		if line[0:6] == 'EXPDTA':
			for tech in tech_list:
				if tech in line:
					PDB_dict['TECHNIQUE'] = tech
			if 'TECHNIQUE' not in PDB_dict:
				PDB_dict['TECHNIQUE'] = 'OTHER'
			print_log(this_name, "Struct {0}    Technique {1}".format(struct, PDB_dict['TECHNIQUE']))
		elif line[0:10] == 'REMARK   2' and 'RESOLUTION' in line:
			fields = line.split()
			for nf in range(len(fields)):
				if 'ANGSTROM' in fields[nf]:
					if fields[nf-1] != "NULL":
						PDB_dict['RESOLUTION'] = float(fields[nf-1])
					else:
						PDB_dict['RESOLUTION'] = 0
			if 'RESOLUTION' not in PDB_dict and PDB_dict['TECHNIQUE'] != 'THEORETICAL' and PDB_dict['TECHNIQUE'] != 'NMR':
				raise NameError("ERROR: Resolution annotation of pdb {0} is badly formatted: {1}".format(struct_filename, line))
		elif line[0:10] == 'REMARK   3' and 'FREE R VALUE' in line and 'ERROR' not in line and 'SET' not in line and (line.split()[3] == 'FREE' or line.split()[3] == 'BIN'):
			try:
				PDB_dict['RFACTOR'] = float(line.split()[-1])
			except ValueError:
				PDB_dict['RFACTOR'] = 'NULL'
		elif line[0:4] == 'ATOM':
			if not line[21]:
				raise NameError("ERROR: There is an ATOM without chain name: {0}".format(line))
			ch_name = line[21]
			if ch_name not in res_ids:
				res_ids[ch_name] = []
				b_factor[ch_name] = 0
				b_norm[ch_name] = 0
			if ch_name in res_ids and (not res_ids[ch_name] or  res_ids[ch_name][-1][0] != int(line[22:26])):
				res_ids[ch_name].append((int(line[22:26]), line[17:20]))
			if line[60:66]:
				b_factor[ch_name] += float(line[60:66])
				b_norm[ch_name] += 1
		elif line[0:5] == 'TITLE':
			if 'TITLE' not in PDB_dict:
				PDB_dict['TITLE'] = line[10:].rstrip()
			else:
				PDB_dict['TITLE'] += line[10:].rstrip()
		elif 'REMARK 350' in line:
			print_log(this_name, "{0}   {1}".format(struct, line))
			if 'BIOMOLECULE:' in line:
				biomolecule_occurrences += 1
			if biomolecule_occurrences > 1:
				continue	
			if 'APPLY THE FOLLOWING TO CHAINS:' in line or 'AND CHAINS:' in line:
				if 'BIOLOGICAL_UNIT' not in PDB_dict:
					PDB_dict['BIOLOGICAL_UNIT'] = set()
				fields = line.split()
				are_chains = False
				for field in fields:
					if are_chains:
						PDB_dict['BIOLOGICAL_UNIT'].add(field[0])
					if field == 'CHAINS:':
						are_chains = True
			elif 'BIOMT1' in line or 'BIOMT2' in line or 'BIOMT3' in line:
				if not 'BIOMATRIX' in PDB_dict:
					PDB_dict['BIOMATRIX'] = {}
				fields = line.split()
				if not int(fields[3])-1 in PDB_dict['BIOMATRIX']:
					PDB_dict['BIOMATRIX'][int(fields[3])-1] = {}
				Mline = int(fields[2][-1])-1
				PDB_dict['BIOMATRIX'][int(fields[3])-1][Mline] = [float(fields[4]), float(fields[5]), float(fields[6]), float(fields[7])]
		elif line[0:6] == 'ENDMDL':
			break

	if 'RFACTOR' not in PDB_dict:
		PDB_dict['RFACTOR'] = 'NULL'

	if 'TECHNIQUE' in PDB_dict:
		if PDB_dict['TECHNIQUE'] == 'NMR' and not 'BIOLOGICAL_UNIT' in PDB_dict:
			PDB_dict['BIOLOGICAL_UNIT'] = set()
			for chain in sorted(list(res_ids.keys())):
				PDB_dict['BIOLOGICAL_UNIT'].add(chain)
	else:
		PDB_dict['TECHNIQUE'] = 'UNKNOWN'

	PDB_dict['CHAINS'] = set()
	for chain in sorted(list(res_ids.keys())):
		PDB_dict['CHAINS'].add(chain)
		if b_factor[chain] == 0:
			b_factor[chain] = 'NULL'
		else:
			b_factor[chain] = b_factor[chain] / b_norm[chain]
		PDB_dict[chain] = {'AVG_BFACTOR' : b_factor[chain],
		                   'NRES' : len(res_ids[chain]),
		                   'RESIDS' : [x[0] for x in res_ids[chain]],
		                   'RESNAMES' : [from3to1(x[1]) for x in res_ids[chain]]}
	return PDB_dict


def FASTA_parser(locations, struct, PDB_dict, force_PDB=False):
	this_name = "FASTA_parser"

	fasta_dict = {}
	if (os.path.exists(locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb') or os.path.exists(locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb')) and not force_PDB:
		fasta_dict['extracted_from_structure'] = True
		parser = PDBParser(PERMISSIVE=True)
		if os.path.exists(locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb'):
			address = locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb'
		elif os.path.exists(locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb'):
			address = locations['FSYSPATH']['repoppmpdb'] + struct + '_opm.pdb'
		structure = parser.get_structure('self', address)
		ppb = PPBuilder()
		for model in structure:
			for chain in model:
				for pp in ppb.build_peptides(chain):
					fasta_dict[chain.id] = '{0}'.format(pp.get_sequence())
			break
	else:
		fasta_dict['extracted_from_structure'] = False
		fasta_filename = locations['FSYSPATH']['PDBfasta'] + struct + '.fa'
		fasta_file = open(fasta_filename, 'r')
		text = fasta_file.read().split("\n")
		fasta_file.close()

		title_line = False
		for line in text:
			if not (line and re.match('^\S', line)):
				title_line = False
				continue
			if line[0] == '>':
				structname = line[1:5]
				if struct.lower() != structname.lower():
					raise_error(this_name, "ERROR: filename and structure name do not match in FASTA file {0}".format(fasta_filename))
				chainname = line[6]
				fasta_dict[chainname] = ""
				title_line = True
			elif line[0] != '>' and title_line:
				fasta_dict[chainname] += line.strip()

	print_log(this_name, "FASTA chains of struct {0}".format(struct))
	count_idc = {}
	tobeconfirmed = set()
	for k in sorted(list(fasta_dict.keys())):
		if len(k) == 1:
			if k not in PDB_dict['BIOLOGICAL_UNIT']:
				print_log(this_name, "Chain {0} not found in biological unit, check for biomatrix compatibility".format(k))
				tobeconfirmed.add(k)
			else:
				count_idc[k] = 1
				print_log(this_name, "{0}    {1}".format(k, fasta_dict[k]))
	for k in sorted(tobeconfirmed):
		found_match = False
		for kk in sorted(list(count_idc.keys())):
			if fasta_dict[k] == fasta_dict[kk] and count_idc[kk] < len(PDB_dict['BIOMATRIX']):
				print_log(this_name, "Chain {0} explained by biomatrix, saved".format(k))
				found_match = True
				count_idc[kk] += 1
				break
		if not found_match:
			print_log(this_name, "Chain {0} not explained by biomatrix, deleted".format(k))
			del fasta_dict[k]
			continue
		print_log(this_name, "{0}    {1}".format(k, fasta_dict[k]))

	print(struct)
	print(fasta_dict)
	return fasta_dict


def retrieve_OPM_family_classification(force_redo=False, data_folder='OPM_families_html/'):
	if not os.path.exists(data_folder):
		os.mkdir(data_folder)

	home_url = "http://opm.phar.umich.edu/"
	proteinlist_url = home_url + "proteins.php"
	proteinlist_filename = data_folder + "proteinlist.txt"
	downloaded = reiterated_simple_download(proteinlist_url, proteinlist_filename)

	if not downloaded:
		raise NameError("Protein list not downloaded", proteinlist_url)

	proteinlist_file = open(proteinlist_filename, 'r')
	soup = BeautifulSoup(proteinlist_file)

	classification = {}
	family_dict = {}

	data = []
	table = soup.find('table', attrs={'class':'data'})
	
	table_header = table.find('thead')
	rows = table_header.find_all('tr')
	for row in rows[:1]:
		cols = row.find_all('th')
		chl = []
		for col in cols:
			chl.append(col.text)

	table_body = table.find('tbody')
	rows = table_body.find_all('tr')
	for row in rows:
		cols = row.find_all('td')
		for col in cols:
			if col == cols[2]:
				struct = col.text
				classification[struct] = {}
		for ncol in range(len(cols)):
			if ncol != 2:
				classification[struct][chl[ncol]] = cols[ncol].text
			if chl[ncol] == 'Family':
				hrefs = cols[ncol].find_all('a', href=True)
				classification[struct]['Family_href'] = home_url + hrefs[0]['href']
				family_dict[classification[struct]['Family']] = classification[struct]['Family_href']
	proteinlist_file.close()

	return classification

def create_pbs(hostname, locations, command_filename, pbs_filename, jobname, npw=16, hwtw=48, mwtw=0, maxppn=16):
	# npw: number of Processors in Worker
	# hwtw: Hour WallTime for Worker

	job_out_filename = locations['FSYSPATH']['logs'] + jobname + '_out.txt'

	log_path = error_path = locations['FSYSPATH']['logs']
	workdir = os.getcwd() + '/'

	command_file = open(command_filename, 'r')
	text = [l for l in command_file.read().split('\n') if l]
	command_file.close()
	if len(text) < npw:
		npw = len(text)
	# This conversion depends on how much memory you want to reserve for each process, and how much RAM do each processor have
	# In LoBoS, the node in Haswell machine have 16 GB and 16 cores, but we want at least 2 GB per process, so 8 processes per node.
	nn = (npw-1)//16 + 1

	if npw >= maxppn:
		ppn = maxppn
	else:
		ppn = npw

	# This text is specific for the Haswell machine in LoBos. Modify when needed:
	# key hwell
	# module paths
	# bash_profile path
	# Note that mpiq.exe must be in the scripts folder.
	if mwtw>0 and mwtw<60:
		pbstext = """#!/bin/bash
#PBS -u esarti
#PBS -N {7}
#PBS -l nodes={0}:hwell:ppn={8}
#PBS -l walltime=00:{1}:00
#PBS -l pmem=2000mb
#PBS -j oe
#PBS -o {2}
#PBS -e {3}
	
module load intel/17.0.0
module load openmpi/2.0.1
source /u/esarti/.bash_profile
source /u/esarti/.bashrc
	
cd {4}

mpirun -np {5} mpiq.exe {6} &> {9}
		""".format(nn, str(mwtw).zfill(2), log_path, error_path, workdir, npw, command_filename, jobname, ppn, job_out_filename)
	else:
		pbstext = """#!/bin/bash
#PBS -u esarti
#PBS -N {7}
#PBS -l nodes={0}:hwell:ppn={8}
#PBS -l walltime={1}:00:00
#PBS -l pmem=2000mb
#PBS -j oe
#PBS -o {2}
#PBS -e {3}
	
module load intel/17.0.0
module load openmpi/2.0.1
source /u/esarti/.bash_profile
source /u/esarti/.bashrc
	
cd {4}

mpirun -np {5} mpiq.exe {6} &> {9}
		""".format(nn, hwtw, log_path, error_path, workdir, npw, command_filename, jobname, ppn, job_out_filename)

	pbs_file = open(pbs_filename, 'w')
	pbs_file.write(pbstext)
	pbs_file.close()
	return jobname, hwtw, job_out_filename

def job_control(locations, walltime_hrs, hostname, job_out_filename, is_done_dict, num_of_jobs, mode):
	workdir = os.getcwd() + '/'

	# qstat -a on a log file
	qstat_logfilename = locations['FSYSPATH']['logs'] + "last_qstat.txt"
	qstat_logfile = open(qstat_logfilename, 'w')
	subprocess.call(["ssh", "{0}".format(hostname), "qstat", "-a"], stdout=qstat_logfile)
	qstat_logfile.close()

	# Read qstat -a
	qstat_logfile = open(qstat_logfilename, 'r')
	text = qstat_logfile.read().split('\n')
	qstat_logfile.close()

	unfinished = set()
	stillon = set()
	for jobname in sorted(list(job_out_filename.keys())):
		status = None
		for line in text:
			fields = line.split()
			if jobname in line:
				status = fields[-2]
		if status == None:
			if mode == 'FrTMAlign':
				job_out_file = open(job_out_filename[jobname], 'r')
				otext = job_out_file.read().split('\n')
				job_out_file.close()

				at_least_one_Start = False
				for oline in otext:
					if 'Start' not in oline:
						continue
					at_least_one_Start = True
				if not at_least_one_Start:
					raise NameError("ERROR: Job {0} crashed before starting. Check {1}".format(jobname, job_out_filename[jobname]))

				for oline in otext:
					if not oline:
						continue
					ofields = oline.split()
					if len(ofields)>=2 and ofields[-2] == 'Start':
						unfinished.add(ofields[-1])
					elif len(ofields)>=2 and ofields[-2] == 'Stop':
						unfinished.remove(ofields[-1])
		else:
			stillon.add(jobname)

	if (not unfinished) and (False not in [is_done_dict[x][1] for x in list(is_done_dict.keys())]):
		return

	if len(stillon) < num_of_jobs:
		i = 0
		current_instrfilename = ''
		while not current_instrfilename:
			tentative_instrfilename = locations['SYSFILES']['instr_filename'](i)
			if not os.path.exists(tentative_instrfilename):
				current_instrfilename = tentative_instrfilename
			i += 1

		current_instrfile = open(current_instrfilename, 'w')
		if mode == 'FrTMAlign':
			for s in sorted(unfinished):
				data_filename = locations['FSYSPATH']['logs'] + s + '.dat'
				fasta_tmpfolder_path_occ = is_done_dict[s][0] + locations['TREE']['straln'] + 'tmp_' + s + '_fasta/' + 'occupied'
				if os.path.exists(fasta_tmpfolder_path_occ):
					os.remove(fasta_tmpfolder_path_occ)
				current_instrfile.write("python FrTMjob_wrapper.py {0} complete\n".format(data_filename))
			for s in sorted([x for x in list(is_done_dict.keys()) if not is_done_dict[x][1]]):
				data_filename = locations['FSYSPATH']['logs'] + s + '.dat'
				current_instrfile.write("python FrTMjob_wrapper.py {0} from_scratch\n".format(data_filename))
			pbs_key = 'FrTMpbs'
		elif mode == 'PPM':
			for s in sorted([x for x in list(is_done_dict.keys()) if not is_done_dict[x][1]]):
				in_filename = s + '_input.dat'
				out_filename = s + '_results.dat'
				tmp_folder = locations['FSYSPATH']['PPM'] + 'tmp_' + s + '/'
				current_instrfile.write("cd {0}; {1} < {2} > {3}\n".format(tmp_folder, locations['OPT']['ppmexe'], in_filename, out_filename))
			pbs_key = 'PPMpbs'
		current_instrfile.close()

		num_new_jobs = num_of_jobs - len(stillon)
		jn = 0
		for jobname in sorted(list(job_out_filename.keys())):
			new_jn = int(jobname.split('_')[-1])
			if jn < new_jn:
				jn = new_jn
			jobname_base = jobname[:-(len(str(new_jn)))]

		for i in range(jn+1, jn+num_new_jobs+1):
			jobname = jobname_base + str(i)
			pbs_filename = locations['SYSFILES'][pbs_key](jobname)
			jobname, walltime_hrs, job_out_filename[jobname] = create_pbs(hostname, locations, current_instrfilename, pbs_filename, jobname=jobname)
			subprocess.call(["ssh", "{0}".format(hostname), "qsub", "{0}".format(pbs_filename)], cwd=workdir)

	return job_out_filename

def ncr(n, r):
#	print("ncr", n, r)
	r = min(r, n-r)
	if r == 0:
		return 1
#	print("reduce", range(n, n-r, -1))
	numer = reduce(op.mul, range(n, n-r, -1))
	denom = reduce(op.mul, range(1, r+1))
	return numer//denom
