# Name: genopmlib.py
# Language: python3
# Libraries:  
# Description: Generates HOMEP raw pdb library with OPM data
# Author: Edoardo Sarti
# Date: Aug 26 2016

from support_opm import *
import re, os, multiprocessing
from bs4 import BeautifulSoup

def retrieve_opm_list(locations):
	opmlib_url = locations['FIXED']['opmliburl']
	local_filename = locations['FSYSPATH']['OPMhtmls'] + 'structure_list.html'
	if not os.path.exists(local_filename):
		print(opmlib_url, local_filename)
		downloaded = reiterated_simple_download(opmlib_url, local_filename)
	
	local_file = open(local_filename, 'r')
	text = local_file.read().split('\n')
	local_file.close()
	struct_list = []
	for line in text:
		if not line:
			continue
		if '<br />' in line and line[0] != '<':
			pdbcode = [x.strip() for x in re.findall(r'\s*(.*)<br />', line, re.DOTALL) if x.strip()]
			if pdbcode:
				struct_list.append(pdbcode[0])
	return struct_list


def crawler(data):
	options, locations, struct_list, this_name, rel_list = data

	refstruct_semiurl = locations['FIXED']['opmrefsemiurl']
	OPM_data = {}
	run_in_PPM = set()
	for struct in struct_list:
		logmsg = print_log(this_name, "Examining structure {0}".format(struct))
		OPM_data[struct] = {}

		# Download HTML page
		refstruct_url = refstruct_semiurl + struct
		info_filename = locations['FSYSPATH']['OPMhtmls'] + struct + '_info.html'
		# This cannot come from repo!
		if not os.path.exists(info_filename):
			downloaded = reiterated_simple_download(refstruct_url, info_filename)
			if downloaded:
				logmsg += print_log(this_name, "OPM HTML page stored in {0}".format(info_filename), quiet=True)
			else:
				logmsg += print_log(this_name, "WARNING: OPM HTML page not found!")
				write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
				del OPM_data[struct]
				continue

		# Download PDB file from PDB
		pdb_url = locations['FIXED']['pdbsemiurl'] + struct + '.pdb.gz'
		pdb_filename = locations['FSYSPATH']['PDBpdbs'] + struct + '.pdb'
		r_pdb_filename = locations['FSYSPATH']['repoPDBpdbs'] + struct + '.pdb'
		if os.path.exists(r_pdb_filename):
			shutil.copy(r_pdb_filename, pdb_filename)
		if not os.path.exists(pdb_filename):
			downloaded = reiterated_gzip_download(pdb_url, pdb_filename)
			if downloaded:
				logmsg += print_log(this_name, "PDB pdb file stored in {0}".format(pdb_filename), quiet=True)
			else:
				logmsg += print_log(this_name, "WARNING: PDB pdb file not found!")
				write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
				del OPM_data[struct]
				continue

		# Download FASTA file from PDB
		fasta_url = locations['FIXED']['pdbfastasemiurl'] + struct
		fasta_filename = locations['FSYSPATH']['PDBfasta'] + struct + '.fa'
		r_fasta_filename = locations['FSYSPATH']['repoPDBfasta'] + struct + '.fa'
		if os.path.exists(r_fasta_filename):
			shutil.copy(r_fasta_filename, fasta_filename)
		if not os.path.exists(fasta_filename):
			downloaded = reiterated_simple_download(fasta_url, fasta_filename)
			if downloaded:
				logmsg += print_log(this_name, "PDB fasta file stored in {0}".format(fasta_filename), quiet=True)
			else:
				logmsg += print_log(this_name, "WARNING: PDB fasta file not found!")
				write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
				del OPM_data[struct]
				continue

		# Download PDB file from OPM
		opmpdb_url = locations['FIXED']['opmpdbsemiurl'] + struct + '.pdb'
		opmpdb_filename = locations['FSYSPATH']['OPMpdbs'] + struct + '.pdb'
		r_opmpdb_filename = locations['FSYSPATH']['repoOPMpdbs'] + struct + '.pdb'
		if os.path.exists(r_opmpdb_filename):
			shutil.copy(r_opmpdb_filename, opmpdb_filename)
		if not os.path.exists(opmpdb_filename):
			downloaded = reiterated_simple_download(opmpdb_url, opmpdb_filename)
			if downloaded:
				logmsg += print_log(this_name, "OPM pdb file stored in {0}".format(opmpdb_filename), quiet=True)
			else:
				logmsg += print_log(this_name, "WARNING: OPM pdb page not found!")
				print("crawler: run in PPM struct", struct)
				print("crawler: because OPM structure couldn't be downloaded")
				run_in_PPM.add(struct)

		# Read HTML file
		info_file = open(info_filename, 'r')
		try:
			text = info_file.read()
			text = text.split('\n')
		except:
			write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
			del OPM_data[struct]
			continue
#		text = info_file.read().split('\n')
		info_file.close()
		activate_namesearch = False
		activate_relatedsearch = False
		activate_clssearch = False
		activate_chainsearch = False
		activate_familysearch = False
		activate_supfamsearch = False
		abort = False
		for line in text:
			if not line:
				continue
			if '<!-- body section -->' in line:
				activate_namesearch = True
				activate_clssearch = True
				activate_familysearch = True
				activate_supfamsearch = True
			if struct in line and activate_namesearch:
				OPM_data[struct]['title'] = [x.strip() for x in re.findall(r'&raquo;(.*?)</h1>$', line, re.DOTALL)]
				if OPM_data[struct]['title']:
					for n in range(3):
						if type(OPM_data[struct]['title']) == list:
							OPM_data[struct]['title'] = OPM_data[struct]['title'][0]
					if OPM_data[struct]['title'].rfind(',') > -1:
						OPM_data[struct]['title'] = OPM_data[struct]['title'][:OPM_data[struct]['title'].rfind(',')]
				activate_namesearch = False
			if 'Class' in line and activate_clssearch:
				soup = BeautifulSoup(line, 'lxml')
				cls = [link.text for link in soup.find_all('a')][0]
				if 'Alpha' in cls:
					OPM_data[struct]['class'] = 'alpha'
				elif 'Beta' in cls:
					OPM_data[struct]['class'] = 'beta'
				else:
					OPM_data[struct]['class'] = 'other'
			if 'Superfamily' in line and activate_supfamsearch:
				soup = BeautifulSoup(line, 'lxml')
				superfamilies = [link.text for link in soup.find_all('a')]
				for s in superfamilies:
					fields = s.split()
					if len(fields) < 2:
						continue
					if 'TCDB' in s:
						OPM_data[struct]['TCDB_superfamily'] = fields[0]
					else:
						OPM_data[struct]['superfamily'] = (fields[0], ''.join(fields[1:]))
				activate_supfamsearch = False
			if 'Family' in line and activate_familysearch:
				soup = BeautifulSoup(line, 'lxml')
				families = [link.text for link in soup.find_all('a')]
				for s in families:
					fields = s.split()
					if len(fields) < 2:
						continue
					if 'TCDB' in s:
						OPM_data[struct]['TCDB_family'] = fields[0]
					else:
						OPM_data[struct]['family'] = (fields[0], ''.join(fields[1:]))
				activate_familysearch = False
			if 'Other PDB entries representing this structure' in line:
				activate_relatedsearch = True
			if '?extrapdb=' in line and activate_relatedsearch:
				soup = BeautifulSoup(line, 'lxml')
				OPM_data[struct]['related'] = [link.text for link in soup.find_all('a') if link.text != struct]
				log = "{0} has {1} related structures: ".format(struct, len(OPM_data[struct]['related']))
				for rel in OPM_data[struct]['related']:
					log += rel + " "
				logmsg += print_log(this_name, log)
				activate_relatedsearch = False
			if 'transmembrane subunit' in line:
				activate_chainsearch = True
			if '<td align=' in line and activate_chainsearch:
				chain = [x.strip() for x in re.findall(r'<b>(\s*.\s*?)</b>', line, re.DOTALL)][0]
				if 'tmchains' not in OPM_data[struct]:
					OPM_data[struct]['tmchains'] = set()
					logmsg += print_log(this_name, "{0} chain list:".format(struct))
				if chain in OPM_data[struct]['tmchains']:
					if chain.lower() == chain:
						logmsg += print_log(this_name, "Warning: two chains with the same name ({0}) found in the OPM html description of struct {1}. Structure aborted.".format(chain, struct))
						write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
						del OPM_data[struct]
						abort = True
						break
					else:
						logmsg += print_log(this_name, "Warning: two chains with the same name ({0}) found in the OPM html description of struct {1}. Second occurrence changed in {2}".format(chain, struct, chain.lower()))
						chain = chain.lower()
				OPM_data[struct]['tmchains'].add(chain)
				print_log(this_name, "tm_chain {0} {1}".format(struct, chain))
				OPM_data[struct][chain] = find_segments(line)
				logmsg += print_log(this_name, "{0} chain {1} has {2} segments".format(struct, chain, OPM_data[struct][chain]['nseg']))
			if '</table>' in line and activate_chainsearch:
				activate_chainsearch = False
		if abort:
			continue
		if not 'tmchains' in OPM_data[struct]:
			logmsg += print_log(this_name, "No TM domains found in {0}".format(struct))
			write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
			del OPM_data[struct]
			continue
		if not 'class' in OPM_data[struct]:
			logmsg += print_log(this_name, "Class specification not found in {0}".format(struct))
			write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
			del OPM_data[struct]
			continue
		if OPM_data[struct]['class'] == 'other':
			logmsg += print_log(this_name, "Structure {0} is neither alpha nor beta".format(struct))
			write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)
			del OPM_data[struct]
			continue
		OPM_data[struct]['is_representative'] = True
		if 'related' in OPM_data[struct] and not options['no_related']:
			if not OPM_data[struct]['related']:
				raise NameError("ERROR: internal")
			activate_relnamesearch = False
			for relstruct in OPM_data[struct]['related']:
				if rel_list and relstruct not in rel_list:
					continue
				relstruct_url = locations['FIXED']['opmrelsemiurl'] + relstruct
				logmsg += print_log(this_name, "------------------------------", quiet=True)
				logmsg += print_log(this_name, "{0} related structure: {1}".format(struct, relstruct))

				if relstruct in OPM_data:
					logmsg += print_log(this_name, "Struct {0} already found as representative".format(relstruct))
					OPM_data[struct]['related'].remove(relstruct)
					continue

				# Download HTML file from OPM
				relinfo_filename = locations['FSYSPATH']['relOPMhtmls'] + relstruct + '_info.html'
				r_relinfo_filename = locations['FSYSPATH']['reporelOPMhtmls'] + relstruct + '_info.html'
				if os.path.exists(r_relinfo_filename):
					shutil.copy(r_relinfo_filename, relinfo_filename)
				if not os.path.exists(relinfo_filename):
					downloaded = reiterated_simple_download(relstruct_url, relinfo_filename)
					if downloaded:
						logmsg += print_log(this_name, "OPM HTML page stored in {0}".format(relinfo_filename), quiet=True)
					else:
						logmsg += print_log(this_name, "WARNING: OPM HTML page not found!")

				# Download PDB file from PDB
				pdb_url = locations['FIXED']['pdbsemiurl'] + relstruct + '.pdb.gz'
				pdb_filename = locations['FSYSPATH']['PDBpdbs'] + relstruct + '.pdb'
				r_pdb_filename = locations['FSYSPATH']['repoPDBpdbs'] + relstruct + '.pdb'
				if os.path.exists(r_pdb_filename):
					shutil.copy(r_pdb_filename, pdb_filename)
				if not os.path.exists(pdb_filename):
					downloaded = reiterated_gzip_download(pdb_url, pdb_filename)
					if downloaded:
						logmsg += print_log(this_name, "PDB pdb file stored in {0}".format(pdb_filename), quiet=True)
					else:
						logmsg += print_log(this_name, "WARNING: PDB pdb file not found!")
						continue

				# Download FASTA file from PDB
				fasta_url = locations['FIXED']['pdbfastasemiurl'] + relstruct
				fasta_filename = locations['FSYSPATH']['PDBfasta'] + relstruct + '.fa'
				r_fasta_filename = locations['FSYSPATH']['repoPDBfasta'] + relstruct + '.fa'
				if os.path.exists(r_fasta_filename):
					shutil.copy(r_fasta_filename, fasta_filename)
				if not os.path.exists(fasta_filename):
					downloaded = reiterated_simple_download(fasta_url, fasta_filename)
					if downloaded:
						logmsg += print_log(this_name, "PDB fasta file stored in {0}".format(fasta_filename), quiet=True)
					else:
						logmsg += print_log(this_name, "WARNING: PDB fasta file not found!")
						continue

				OPM_data[relstruct] = {}
				for key in 'family', 'TCDB_family', 'superfamily', 'TCDB_superfamily', 'class':
					if key in OPM_data[struct]:
						OPM_data[relstruct][key] = OPM_data[struct][key]

				# The HTML page is only used to read the title, which is otherwise imported from the representative structure
				try:
					relinfo_file = open(relinfo_filename, 'r')
					rtext = relinfo_file.read()
					rtext = rtext.split('\n')
					for rline in rtext:
						if '<!-- body section -->' in rline:
							activate_relnamesearch = True
						if relstruct in rline and activate_relnamesearch:
							OPM_data[relstruct]['title'] = [x.strip() for x in re.findall('&raquo;(.*?)</h1>', rline, re.DOTALL)]
							if OPM_data[relstruct]['title']:
								for n in range(3):
									if type(OPM_data[relstruct]['title']) == list:
										OPM_data[relstruct]['title'] = OPM_data[relstruct]['title'][0]
							activate_relnamesearch = False
					relinfo_file.close()
				except:
					OPM_data[relstruct]['title'] = OPM_data[struct]['title']

				OPM_data[relstruct]['is_representative'] = False
				OPM_data[relstruct]['related_to'] = struct


				print_log(this_name, "Run in PPM struct {0}".format(relstruct))
				run_in_PPM.add(relstruct)

#				# OPM pdb files for related structures are not downloaded, since they lack anyways of important information
#				if not os.path.exists(locations['FSYSPATH']['repoppmpdb'] + relstruct + '_opm.pdb'):
#					print("crawler: run in PPM struct", relstruct)
#					print("crawler: because it is a related structure without backup")
#					run_in_PPM.add(relstruct)

		logmsg += print_log(this_name, "==============================\n", quiet=True)
		write_log(locations['FSYSPATH']['logs'] + this_name, logmsg)

		if len(struct_list) == 1:
			write_data(locations['FSYSPATH']['logs'] + struct_list[0] + '_OPMdata.log', OPM_data)
	return OPM_data, run_in_PPM


def update_raw_pdb_library(options, locations):
	struct_list = retrieve_opm_list(locations)

	data = []
	i = 0
	if os.path.exists(locations['FSYSPATH']['logs'] + 'crawler.log'):
		shutil.move(locations['FSYSPATH']['logs'] + 'crawler.log', locations['FSYSPATH']['logs'] + 'crawler_old.log')
	for s in struct_list:
		if not os.path.exists(locations['FSYS']['whole'] + s + '.pdb'):
			log_filename = 'crawler_'+str(i)
			data.append((locations, [s]))
			if os.path.exists(locations['FSYSPATH']['logs'] + log_filename):
				os.remove(locations['FSYSPATH']['logs'] + log_filename)
			i += 1
#			data.append


def generate_raw_pdb_library(options, locations):
	struct_list = retrieve_opm_list(locations)

#	toni_list = ["1okc", "1af6", "1pv6", "2nwl", "3tij", "3v8g"] ## TEMPORARY!!! REMOVE!!!
#	toni_list = ["5c78", "4z7f", "1lda", "3kly", "4pj0", "3tdp", "2ycw"] #CONFORMATIONAL
#	toni_list = ["5c78", "4z7f", "1lda", "3kly", "3tdp", "2ycw"] #CONFORMATIONAL2
#	toni_list = ["4fz1", "1qj8", "2fgq", "1ezv", "1ijd", "1l0l", "1nek", "1q16", "2a65", "2b6o", "2ns1", "2x72", "3f7v", "3f7y", "3fb6", "3fb7", "3m9i", "3oe0", "3t56", "3wme", "3x3b", "4fxz", "4fz1", "4h1d", "4pxf", "4j05", "4kfm", "4ntx", "4e7c"] #ALLTEST
#	toni_list = ["1rwt", "1yew", "2yxq", "2z73", "4c9h", "4fxz", "4kyt", "4o9p", "5azd", "5fvn", "5jae"]
#	if options['hostname']:
#		struct_list = struct_list[1:50]

	data_struct_list = []
	log_filename = 'crawler_0'
	for s in struct_list:
#		if s not in toni_list:   ## TEMPORARY!!! REMOVE!!!
#			continue	## TEMPORARY!!! REMOVE!!!
		data_struct_list.append(s)
		if os.path.exists(locations['FSYSPATH']['logs'] + log_filename):
			os.remove(locations['FSYSPATH']['logs'] + log_filename)

	tic = time.time()
#	rel_list = ["3usi"]  ### MUST BE VOID FOR THE WHOLE RUN!!!!
	rel_list = []
	output, run_in_PPM = crawler((options, locations, data_struct_list, log_filename, rel_list))

	totlog_filename = locations['FSYSPATH']['logs'] + 'crawler.log'
	with open(totlog_filename, 'w') as totlog_file:
		for struct in data_struct_list:
			in_filename = locations['FSYSPATH']['logs'] + struct + '.log'
			if os.path.exists(in_filename):
				with open(in_filename, 'r') as in_file:
					for line in in_file:
						totlog_file.write(line)
				os.remove(in_filename)

	# OPM_data contains as entries only reference structures! A list of related structures is appended.	
	OPM_data = {}
	for key in list(output.keys()):
		OPM_data[key] = output[key]
	toc = time.time()
	print("crawler: time spent", toc-tic)

	# Generate an archive file here
	tic = time.time()
	OPM_archive(locations, OPM_data)
	toc = time.time()
	print("OPM_archive: time spent", toc-tic)

	# Write general archive file
	write_data(locations['SYSFILES']['opmdata'], OPM_data)

	return OPM_data, run_in_PPM
