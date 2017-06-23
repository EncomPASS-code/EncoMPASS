# Name: genfsys_opm_symmetry.py
# Language: python 3.5/2.7
# Libraries: sys, os, datetime, argparse; support_opm
# Description: Generates the excerpt of the HoTMeSS file system that affects symmetry
# Author: Edoardo Sarti; Antoniya A. Aleksandrova
# Date: Dec 14 2016

import sys
import os
import datetime
import argparse
import collections
from support_opm import *


# Command line parser
# Return two dictionaries:
#  options - contains all options specified in the command line and/or in 
#            defaults
#  filters - contains all optional or adjustable criteria of selection that
#            will be used to filter the database in genclib.
def main_parser():
	this_name = 'main_parser'
	parser = argparse.ArgumentParser()

	# Compulsory arguments
	parser.add_argument('-pdbtm', '--pdbtm_file_path', nargs=1)

	# Either of these must be set
	parser.add_argument('-d', '--install_dir', nargs=1)  # Either of these
	parser.add_argument('-m', '--main_dir', nargs='?')   # must be set

	# These are compulsory if -d is set, will not be considered if -m is set
	parser.add_argument('-s', '--straln_path', nargs='?')
	parser.add_argument('-ppm', '--ppm_path', nargs='?')
	parser.add_argument('-seq', '--seqaln_path', nargs='?')
	parser.add_argument('-np', '--number_of_procs', nargs='?')
	parser.add_argument('-st', '--subunit_thr', nargs='?')
	parser.add_argument('-ct', '--cluster_thr', nargs='?')
	parser.add_argument('-rf', '--resolution_filter', nargs='?')

	# Optional
	parser.add_argument('-ht', '--hole_thr', nargs='?')
	parser.add_argument('-oh', '--output_homep', nargs='?')
	parser.add_argument('-otab', '--output_tab', nargs='?')
	parser.add_argument('-dpcl', '--DPClustering', nargs='?')
	parser.add_argument('-repo', '--repository', nargs='?')
	parser.add_argument('--oldtable', nargs='?')
	parser.add_argument('-host', '--hostname', nargs='?')
	parser.add_argument('-is', '--instraln', nargs='?')
	parser.add_argument('-with_nmr', action='store_true')
	parser.add_argument('-with_theoretical', action='store_true')
	parser.add_argument('-no_overwrite', action='store_true')
	parser.add_argument('-no_related', action='store_true')


	# Default values for optional arguments
	parser.set_defaults(install_dir = '')
	parser.set_defaults(main_dir = '')
	parser.set_defaults(hole_thr = '100')
	parser.set_defaults(output_tab = 'structure_alignments.dat')
	parser.set_defaults(output_homep = 'HOMEP4.dat')
	parser.set_defaults(repository = None)
	parser.set_defaults(oldtable = '')

	parsed = parser.parse_args()


	# Create 'options' dictionary containing all selected options as
	# strings
	options = {}
	for x in sorted(parsed.__dict__):
		if type(parsed.__dict__[x]) == list:
#			print(x, parsed.__dict__[x][0])
			options[x] = parsed.__dict__[x][0]
		else:
#			print(x, parsed.__dict__[x])
			options[x] = parsed.__dict__[x]

	for x in list(options.keys()):
		if options[x] == 'False':
			options[x] = False
		elif options[x] == 'True':
			options[x] = True

	# Either flag -d or -m must be set
	if not (options['install_dir'] or options['main_dir']):
		raise_error(this_name, "ERROR: either -d (--install_dir) or -m (--main_dir) must be specified")

	if options['install_dir'] and options['no_overwrite']:
		raise_error(this_name, "ERROR: option -no_overwrite cannot be combined with option -d (--install_dir)")

	# If -d flag is set, these flags must be set too
	if options['install_dir'] and not (options['seqaln_path'] or options['straln_path'] or options['number_of_procs'] or 
	 options['object_thr'] or options['cluster_thr'] or options['resolution_filter'] or options['ppm_path']):
		raise_error(this_name, "ERROR: if flag -d (--install_dir) is set, these flags must be set too:\n" +
		                       "       -s   (--straln_path) <structure_aln_program_path>\n" +
		                       "       -seq (--seqaln_path) <sequence_aln_program_path>\n" +
		                       "       -ppm (--ppm_path) <PPM_program_path>\n" +
		                       "       -np  (--number_of_procs) <number of processors>\n" +
		                       "       -ot  (--object_thr) <seq id threshold to determine Object>\n" +
		                       "       -ct  (--cluster_thr) <TM-score threshold to clusterize Objects>\n" +
		                       "       -rf  (--resolution_filter) <resolution filter (in A) to select chains>\n")

	# Create 'filters' dictionary containing all optional and tunable
	# criteria of selection for genclib
	filters = {}
	if options['resolution_filter'] != None:
		filters['resolution'] = float(options['resolution_filter'])
	if options['with_nmr'] != None:
		filters['NMR'] = options['with_nmr']
	if options['with_theoretical'] != None:
		filters['THM'] = options['with_theoretical']
	if int(options['hole_thr']) != None:
		filters['hole_thr'] = int(options['hole_thr'])

	print(options, filters)

	return options, filters


def write_hidden_files(locations):



	# Write the .locations.dat hidden file in the main directory of the file system.
	# It contains the same information contained in the 'locations' dictionary.
	locations_filename = locations['SYSFILES']['H_locations']
	locations_file = open(locations_filename, 'w')
	for key, value in list(locations.items()):
		for kkey, vvalue in list(value.items()):
			locations_file.write("{0}\t\t{1}\t\t{2}\n".format(key, kkey, vvalue))
	locations_file.close()

	return

def read_and_merge(filenames, dictionaries):
	locations_filename = filenames
	locations = dictionaries

	# Read the .locations.dat file and returns the 'locations' dictionary
	locations = {}
	locations_file = open(locations_filename, 'r')
	text = locations_file.read().split('\n')
	for line in text:
		if line:
			fields = line.split()
			if fields[0] not in locations:
				locations[fields[0]] = {}
			locations[fields[0]][fields[1]] = fields[2]

	return locations


# --- Library functions ---------------------------------------------------- #

# Generate the file system for the HOMEP database 
def generate_filesystem():
	# Hardcoded variables
	this_name = 'generate_filesystem'
	indent = " "*len(header(this_name))
	version = 4
	other_versions_allowed = True

	# Run command line parser
#	options, filters = main_parser()
	home=os.getcwd()
	install_path = os.path.abspath(home) + '/'
	main_dir = 'HoTMeSS_' + str(version) + '_' + datetime.datetime.now().strftime("%Y_%m_%d") + '/'
	main_path = install_path + main_dir

	# Run checks over names and addresses
	if not os.path.exists(install_path):
		raise_error(this_name, "ERROR: The installation directory path {0} does not exist. Please specify an existing path.".format(install_path))
	if os.path.exists(main_path):
		raise_error(this_name, "ERROR: In the installation directory path {0} there already is a folder named {1}\n".format(install_path, main_dir))
	for filename in os.listdir(install_path):
		if filename[0:5] == 'HoTMeSS' and not other_versions_allowed:
			raise_error(this_name, "ERROR: In the installation directory path {0} there are other versions of HoTMeSS.\n".format(install_path) +
			              indent + "       If you want to continue, you have to set the internal variable other_versions_allowed as True.")

	# Create 'locations' nested dictionary of ordered dictionaries.
	# Under the keyword 'FSYS' should go all names relative to the file system;
	# Under the keyword 'FSYSPATH' should go all paths relative to the file system (they are one more than the FSYS entries, since there is also the install path;
	# Under the keyword 'TREE' should go all denominations of the recurrent sequence/structure/alignment tree
	# Under the keyword 'SYSFILES' should go all system files (hidden system files have keys starting with "H_")
	# Under the keyword 'OPT' should go all other locations and paths it's convenient to save.
	locations = {'TREE' : collections.OrderedDict(), 'FSYS' : collections.OrderedDict(), 'FSYSPATH' : collections.OrderedDict(), 'SYSFILES': collections.OrderedDict(), 'OPT' : collections.OrderedDict(), 'FIXED' : collections.OrderedDict()}
	# TREE

	# FSYS
	locations['FSYS']['main'] = main_dir
	locations['FSYS']['database'] = 'database/'                                                         # database/
	locations['FSYS']['layers'] = locations['FSYS']['database'] + 'layers/'                             # database/layers/
#	locations['FSYS']['tree'] = locations['FSYS']['layers'] + 'tree/'                                   # database/layers/tree/
#	locations['FSYS']['alpha'] = locations['FSYS']['tree'] + 'alpha/'                                   # database/layers/tree/alpha/
#	locations['FSYS']['beta'] = locations['FSYS']['tree'] + 'beta/'                                     # database/layers/tree/beta/
	locations['FSYS']['symmetries'] = locations['FSYS']['layers'] + 'symmetries/'                       # database/layers/symmetries/
#	locations['FSYS']['sequences'] = locations['FSYS']['layers'] + 'sequences/'                         # database/layers/sequences/
	locations['FSYS']['selection'] = locations['FSYS']['database'] + 'selection/'                       # database/selection/
	locations['FSYS']['whole'] = locations['FSYS']['selection'] + 'whole_structs/'                      # database/selection/whole_structs/
	locations['FSYS']['chains'] = locations['FSYS']['selection'] + 'chains/'                            # database/selection/chains/
	locations['FSYS']['old'] = locations['FSYS']['selection'] + '.old/'                                 # database/selection/.old/
	locations['FSYS']['repository'] = 'repository/'                                                     # repository/
	locations['FSYS']['repoPDBpdbs'] = locations['FSYS']['repository'] + 'PDBpdbs/'
	locations['FSYS']['repoPDBfasta'] = locations['FSYS']['repository'] + 'PDBfasta/'
	locations['FSYS']['repoOPMpdbs'] = locations['FSYS']['repository'] + 'OPMpdbs/'
	locations['FSYS']['repoOPMhtmls'] = locations['FSYS']['repository'] + 'OPMhtmls/'
	locations['FSYS']['repowhole'] = locations['FSYS']['repository'] + 'whole_structs/'                 # repository/whole_structs/
	locations['FSYS']['repowholegif'] = locations['FSYS']['repowhole'] + 'gifs/'                        # repository/whole_structs/gifs/
	locations['FSYS']['repochains'] = locations['FSYS']['repository'] + 'chains/'                       # repository/chains/
	locations['FSYS']['repochainsgif'] = locations['FSYS']['repochains'] + 'gifs/'                      # repository/chains/gifs/
#	locations['FSYS']['repocaln'] = locations['FSYS']['repochains'] + locations['TREE']['aln']          # repository/chains/alignments/
#	locations['FSYS']['repocseqaln'] = locations['FSYS']['repochains'] + locations['TREE']['seqaln']    # repository/chains/alignments/seq_alns/
#	locations['FSYS']['repocstraln'] = locations['FSYS']['repochains'] + locations['TREE']['straln']    # repository/chains/alignments/str_alns/
	locations['FSYS']['repoppm'] = locations['FSYS']['repository'] + 'PPM/'                             # repository/PPM/
	locations['FSYS']['repoppmresults'] = locations['FSYS']['repoppm'] + 'results/'                     # repository/PPM/results/
	locations['FSYS']['repoppmdatasub1'] = locations['FSYS']['repoppm'] + 'datasub1/'                   # repository/PPM/datasub1/
	locations['FSYS']['repoppmdatapar1'] = locations['FSYS']['repoppm'] + 'datapar1/'                   # repository/PPM/datapar1/
	locations['FSYS']['repoppmpdb'] = locations['FSYS']['repoppm'] + 'pdbs/'                            # repository/PPM/pdbs/
	locations['FSYS']['PDB'] = 'PDB/'                                                                   # PDB/
	locations['FSYS']['PDBpdbs'] = locations['FSYS']['PDB'] + 'pdbs/'                                   # PDB/pdbs/
	locations['FSYS']['PDBfasta'] = locations['FSYS']['PDB'] + 'fasta/'                                 # PDB/fasta/
	locations['FSYS']['OPM'] = 'OPM/'                                                                   # OPM/
	locations['FSYS']['OPMpdbs'] = locations['FSYS']['OPM'] + 'pdbs/'                                   # OPM/pdbs/
	locations['FSYS']['OPMhtmls'] = locations['FSYS']['OPM'] + 'htmls/'                                 # OPM/htmls/
	locations['FSYS']['PPM'] = locations['FSYS']['OPM'] + 'PPM/'                                        # OPM/PPM/
	#logs/
	locations['FSYS']['logs'] = 'logs/'                                                                 # logs/
	locations['FSYS']['symmetry_log']=locations['FSYS']['logs']+'symmetry/'                             # logs/symmetry
	locations['FSYS']['cesymm_log']=locations['FSYS']['symmetry_log']+'cesymm/'                         # logs/symmetry/cesymm
	locations['FSYS']['cesymm_low_thr_log']=locations['FSYS']['symmetry_log']+'cesymm_lower_threshold/'      # logs/symmetry/cesymm_lower_threshold
	locations['FSYS']['cesymm_from_symd_log']=locations['FSYS']['symmetry_log']+'cesymm_from_symd/'     # logs/symmetry/cesymm_from_symd
	locations['FSYS']['symd_log']=locations['FSYS']['symmetry_log']+'symd/'                             # logs/symmetry/symd
	locations['FSYS']['webarchive'] = 'webarchive/'							    						# webarchive/
	locations['FSYS']['wholegif'] = locations['FSYS']['webarchive'] + 'whole_structs_gifs/'             # webarchive/whole_structs_gifs/
	locations['FSYS']['chainsgif'] = locations['FSYS']['webarchive'] + 'chains_gifs/'                   # webarchive/chains_gifs/
	locations['FSYS']['symm_wholepngs']=locations['FSYS']['webarchive'] + 'symm_whole_structs_pngs/'     # webarchive/symm_whole_structs_pngs/
	locations['FSYS']['symm_chainspngs'] = locations['FSYS']['webarchive'] + 'symm_chains_pngs/'         # webarchive/symm_chains_pngs/
	locations['FSYS']['symmetry'] = 'symmetry/'							    							# symmetry/
	locations['FSYS']['cesymm'] = locations['FSYS']['symmetry'] + 'cesymm/'                   			# symmetry/cesymm/
	locations['FSYS']['symd'] = locations['FSYS']['symmetry'] + 'symd/'                   				# symmetry/symd/
	locations['FSYS']['cesymm_low_thr'] = locations['FSYS']['symmetry'] + 'cesymm_lower_threshold/'   	# symmetry/cesymm_lower_threshold/
	locations['FSYS']['cesymm_from_symd'] = locations['FSYS']['symmetry'] + 'cesymm_from_symd/'   		# symmetry/cesymm_from_symd/
	locations['FSYS']['symm_tables'] = locations['FSYS']['symmetry'] + 'tables/'   						# symmetry/tables/	
	

	# FSYSPATH
	locations['FSYSPATH']['install'] = install_path
	for name, val in locations['FSYS'].items():
		if name == 'main':
			locations['FSYSPATH'][name] = locations['FSYSPATH']['install'] + val
		else:
			locations['FSYSPATH'][name] = locations['FSYSPATH']['install'] + locations['FSYS']['main'] + val
	# SYSFILES
#	locations['SYSFILES']['H_options'] = locations['FSYSPATH']['main'] + '.options.dat'
#	locations['SYSFILES']['H_filters'] = locations['FSYSPATH']['main'] + '.filters.dat'
	locations['SYSFILES']['H_locations'] = locations['FSYSPATH']['main'] + '.locations.dat'
	locations['SYSFILES']['H_topologytype'] = locations['FSYSPATH']['database'] + '.topology_classification.dat'
#	locations['SYSFILES']['H_instraln'] = locations['FSYSPATH']['database'] + '.straln_instructions.dat'
#	locations['SYSFILES']['H_scheduledalns'] = locations['FSYSPATH']['database'] + '.scheduled_alignments.dat'
#	locations['SYSFILES']['H_repolocations'] = locations['FSYSPATH']['repository'] + '.repo_locations.dat'
	locations['SYSFILES']['OPMarchive'] = locations['FSYSPATH']['OPM'] + 'OPM_archive.dat'
#	locations['SYSFILES']['excludedwhole'] = locations['FSYSPATH']['whole'] + 'exclusions.txt'
#	locations['SYSFILES']['excludedchains'] = locations['FSYSPATH']['chains'] + 'exclusions.txt'
#	locations['SYSFILES']['excludedstralnchains'] = locations['FSYSPATH']['chains'] + 'exclusions_from_instructions.txt'
#	locations['SYSFILES']['equivalences'] = locations['FSYSPATH']['chains'] + 'equivalences.txt'
#	locations['SYSFILES']['chaindata'] = locations['FSYSPATH']['chains'] + 'chain_database.txt'
#	locations['SYSFILES']['missingpdbfiles'] = locations['FSYSPATH']['PDBpdbs'] + 'missing_files.txt'
#	locations['SYSFILES']['missingfastafiles'] = locations['FSYSPATH']['PDBfasta'] + 'missing_files.txt'
#	locations['SYSFILES']['repocstraln'] = locations['FSYSPATH']['repocstraln'] + 'structure_alignments.dat'
#	locations['SYSFILES']['repocseqaln'] = locations['FSYSPATH']['repocseqaln'] + 'sequence_alignments.dat'
#	locations['SYSFILES']['webwhole'] = locations['FSYSPATH']['webarchive'] + 'web_archive_whole.xml'
#	locations['SYSFILES']['webchain'] = locations['FSYSPATH']['webarchive'] + 'web_archive_chain.xml'
	locations['SYSFILES']['symmwebwhole'] = locations['FSYSPATH']['webarchive'] + 'symmetry_web_archive_whole.xml'
	locations['SYSFILES']['symmwebchain'] = locations['FSYSPATH']['webarchive'] + 'symmetry_web_archive_chain.xml'
	locations['SYSFILES']['opmdata'] = locations['FSYSPATH']['main'] + '.opm_archive.txt'
	locations['SYSFILES']['homepdata'] = locations['FSYSPATH']['main'] + '.homep_archive.txt'
	locations['SYSFILES']['homeptable'] = locations['FSYSPATH']['main'] + '.homep_table_archive.txt'
#	locations['SYSFILES']['pbsfile'] = locations['FSYSPATH']['main'] + 'allin1.pbs'
#	locations['SYSFILES']['FrTMpbs'] = locations['FSYSPATH']['main'] + 'FrTM_allin1.pbs'
#	locations['SYSFILES']['clusterinstr'] = locations['FSYSPATH']['main'] + 'cluster_commands.txt'
#	locations['SYSFILES']['blosum62'] = locations['FSYSPATH']['main'] + 'blosum62.mat'
	locations['SYSFILES']['pdb_tmchains'] = locations['FSYSPATH']['symmetry'] + '.pdb_tmchains.dat'
	locations['SYSFILES']['cesymm_completed'] = locations['FSYSPATH']['cesymm'] + 'cesymm_completed.list'
	locations['SYSFILES']['symd_completed'] = locations['FSYSPATH']['symd'] + 'symd_completed.list'
	locations['SYSFILES']['cesymm_low_thr_completed'] = locations['FSYSPATH']['cesymm_low_thr'] + 'cesymm_low_thr_completed.list'
	locations['SYSFILES']['cesymm_from_symd_completed'] = locations['FSYSPATH']['cesymm_from_symd'] + 'cesymm_from_symd_completed.list'

	# OPT
#	if options['straln_path']:
#		locations['OPT']['stralnexe'] = options['straln_path']
#	if options['ppm_path']:
#		locations['OPT']['ppmexe'] = options['ppm_path']
#	if options['seqaln_path']:
#		locations['OPT']['seqalnexe'] = options['seqaln_path']
	locations['OPT']['cesymmpath']=install_path+'NIH_work/cesymm-2.0.0-SNAPSHOT/runCESymm.sh'
	locations['OPT']['symdpath']=install_path+'NIH_work/symd1.61/src/symd1.61-linux'
	locations['OPT']['temppath']=install_path+'hotmess_scripts/temp/'
	# FIXED
#	locations['FIXED']['opmliburl'] = 'http://opm.phar.umich.edu/classes_dl.php?type=1'
#	locations['FIXED']['opmrefsemiurl'] = 'http://opm.phar.umich.edu/protein.php?pdbid='
#	locations['FIXED']['opmrelsemiurl'] = 'http://opm.phar.umich.edu/protein.php?extrapdb='
#	locations['FIXED']['opmpdbsemiurl'] = 'http://opm.phar.umich.edu/pdb/'
#	locations['FIXED']['pdbmainsemiurl'] = 'http://www.rcsb.org/pdb/explore/explore.do?structureId='
#	locations['FIXED']['pdbsemiurl'] = 'http://www.rcsb.org/pdb/files/'
#	locations['FIXED']['pdbfastasemiurl'] = 'http://www.rcsb.org/pdb/files/fasta.txt?structureIdList='


	# Generate filesystem
	log = ""
	os.mkdir(locations['FSYSPATH']['main'])
	log += print_log(this_name, "Main directory created: {0}\n".format(locations['FSYSPATH']['main']))

	c = 0
	for index, duple in enumerate(locations['FSYSPATH'].items()):
		if index > 1:
			os.mkdir(duple[1])
			log += print_log(this_name, "Directory {0} has been created.".format(duple[0]))
	write_log(this_name, log)
	write_hidden_files(locations)
	return locations
'''
	# Import repository
	if 'repository' in options and options['repository']:
		if not os.path.exists(options['repository']):
			raise_error(this_name, "ERROR: external repository {0} not found".format(options['repository']))
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
		for loc in list(extrepolocations.keys()):
			if loc in locations['FSYSPATH']:
				loc_files = os.listdir(extrepolocations[loc])
				for filename in loc_files:
					full_filename = os.path.join(extrepolocations[loc], filename)
					if (os.path.isfile(full_filename)):
						if os.path.exists(locations['FSYSPATH'][loc] + filename):
							repo_dict = repo_inspector(locations['FSYSPATH'][loc] + filename)
							ext_repo_dict = repo_inspector(full_filename)
							repo_dict.update(ext_repo_dict)
							for line in list(repo_dict.keys())[0].split('\n'):
								if line and line.split()[0] == 'BEGIN':
									chain_1 = line.split()[2]
									break
							repo_2dict = {chain_1 : repo_dict}
							write_on_repo(locations['FSYSPATH'][loc] + filename, repo_2dict)
						else:
							shutil.copy(full_filename, locations['FSYSPATH'][loc])				
'''
	# Write the three hidden files	



# Retrieves all infromation about the file system

def filesystem_info(main_dir):
	# Hardcoded variables
	this_name = 'filesystem_info'
	indent = " "*len(header(this_name))
	version = 4

	# Run command line parser

	main_path = main_dir
	
	# Perform checks on paths
	if not os.path.exists(main_path):
		raise_error(this_name, "ERROR: Main directory {0} not found.".format(main_path))

	# If the .filters.dat, .locations.dat or the .options.dat files are not found, it returns error
	locations_filename = main_path + '/.locations.dat'
	if not os.path.exists(locations_filename):
		raise_error(this_name, "ERROR: File {0} not found. Filesystem corrupted.".format(locations_filename))


	locations = {}
	locations = read_and_merge(locations_filename, locations)
#	if not options['no_overwrite']:
#		archive_old_file(locations, [options_filename, filters_filename, locations_filename])
#		write_hidden_files(options, filters, locations)

	return locations

if __name__=='__main__':
	generate_filesystem()
