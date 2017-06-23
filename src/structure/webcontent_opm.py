# Name: webcontent_external.py
# Language: python3
# Libraries: 
# Description: Creates the web archive 
# Author: Edoardo Sarti
# Date: Aug 17 2016

import genfsys_opm
import math
import copy
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from support_opm import *
import pylab as pl
from PIL import Image

def neighborlists_creator(options, locations, opm_data, table):
	si_thr = float(options['subunit_thr'])
	tm_thr = float(options['cluster_thr'])

	if not opm_data:
		opm_data = read_data(locations['SYSFILES']['opmdata'])

	print(si_thr, tm_thr)

	# Compile table for whole struct page

	neighbors = {}
	convab = {'a' : 'alpha', 'b' : 'beta'}
	se_n = {}
	st_n = {}
	tot_n = {}
	if not table:
		table = {}
		sumtablefilename = locations['FSYSPATH']['main'] + options['output_tab']
		sumtablefile = open(sumtablefilename, 'r')
		text = sumtablefile.read().split('\n')
		sumtablefile.close()
		for line in text:
			if not line:
				continue
			tmpcl, tt, n1, n2, sesi, stsi, tm, rmsd = line.split()
			if n1 not in neighbors:
				se_nn = 0
				st_nn = 0
				tot_nn = 0
				cl = convab[tmpcl]
				se_n[n1] = []
				st_n[n1] = []
				tot_n[n1] = []
				tt = int(tt)
			else:
				cl, tt, se_nn, st_nn, tot_nn = neighbors[n1]
			sesi = float(sesi)
			tm = float(tm)
			if sesi > si_thr:
				se_nn += 1
				se_n[n1].append(n2)
			if tm > tm_thr:
				st_nn += 1
				st_n[n1].append(n2)
			if sesi > si_thr or tm > tm_thr:
				tot_nn += 1
				tot_n[n1].append(n2)

			### TEST
#			if not (cl == 'alpha' and tt >= 14):
#				continue
			###

			if cl not in list(table.keys()):
				table[cl] = {}
			if tt not in list(table[cl].keys()):
				table[cl][tt] = {}
			if n1 not in list(table[cl][tt].keys()):
				table[cl][tt][n1] = {}

			table[cl][tt][n1][n2] = (sesi, float(stsi), tm, float(rmsd))
			neighbors[n1] = (cl, tt, se_nn, st_nn, tot_nn)
	else:
		for cl in sorted(list(table.keys())):
			for tt in sorted(list(table[cl].keys())):
				for s1 in sorted(list(table[cl][tt].keys())):
					se_nn = 0
					st_nn = 0
					tot_nn = 0
					se_n[s1] = []
					st_n[s1] = []
					tot_n[s1] = []
					for s2 in sorted(list(table[cl][tt][s1].keys())):
						if table[cl][tt][s1][s2][0] > si_thr:
							se_nn += 1
							se_n[s1].append(s2)
						if table[cl][tt][s1][s2][2] > tm_thr:
							st_nn += 1
							st_n[s1].append(s2)
						if table[cl][tt][s1][s2][0] > si_thr or table[cl][tt][s1][s2][2] > tm_thr:
							tot_nn += 1
							tot_n[s1].append(s2)
						neighbors[s1] = (cl, int(tt), se_nn, st_nn, tot_nn)				


	equivalence_filename = locations['SYSFILES']['equivalences']
	equivalences = {}
	if os.path.exists(equivalence_filename):
		equivalence_file = open(equivalence_filename, 'r')
		text = equivalence_file.read().split('\n')
		for line in text:
			if not line:
				continue
			fields = line.split()
			equivalences[fields[0]] = fields[1]
		equivalence_file.close()

	neighborstablefilename = locations['SYSFILES']['neighborstable']
	neighborstablefile = open(neighborstablefilename, 'w')
	neighborstablefile.write("{0:10} {1:10} {2:4}    {3:18} {4:18} {5:18}\n".format("structure", "class", "NTM", "SequenceNeighbors", "StructureNeighbors", "TotalNeighbors"))
	seneighborslistfilename = locations['SYSFILES']['neighborsselist']
	seneighborslistfile = open(seneighborslistfilename, 'w')
	stneighborslistfilename = locations['SYSFILES']['neighborsstlist']
	stneighborslistfile = open(stneighborslistfilename, 'w')
	totneighborslistfilename = locations['SYSFILES']['neighborstlist']
	totneighborslistfile = open(totneighborslistfilename, 'w')
	names = {}
	for struct in sorted(list(neighbors.keys())):
		pdbname = struct[:4]
		if pdbname not in list(names.keys()):
			names[pdbname] = set()
		names[pdbname].add(struct)
		cl, tt, se_nn, st_nn, tot_nn = neighbors[struct]
		seneighborslistfile.write("{0:10} {1}\n".format(struct, "".join([x+';' for x in sorted(se_n[struct])])))
		stneighborslistfile.write("{0:10} {1}\n".format(struct, "".join([x+';' for x in sorted(st_n[struct])])))
		totneighborslistfile.write("{0:10} {1}\n".format(struct, "".join([x+';' for x in sorted(tot_n[struct])])))
		neighborstablefile.write("{0:10} {1:10} {2:4d}    {3:18d} {4:18d} {5:18d}\n".format(struct, cl, tt, se_nn, st_nn, tot_nn))
	neighborstablefile.close()
	seneighborslistfile.close()
	stneighborslistfile.close()
	totneighborslistfile.close()

#	if True:
#		return neighbors, names, opm_data, table

	# Histogram of number of TM chains per complex
	NTMCfilename = locations['FSYSPATH']['webarchive'] + 'tmchains_hist.dat'
	NTMC_hist = []
	NTMC_dict = {}
	FTMC_hist = []
	ncmax = 0
	for pdbname in sorted(list(names.keys())):
		if pdbname not in opm_data or 'tmchains' not in opm_data[pdbname]:
			continue
		nc = len(opm_data[pdbname]['FROM_PDB']['BIOMATRIX'])*len(opm_data[pdbname]['FROM_PDB']['BIOLOGICAL_UNIT'])
		if nc > ncmax:
			ncmax = nc
	stackhist = - np.ones((len(list(names.keys())), ncmax))
	stackcounter = 0
	for pdbname in sorted(list(names.keys())):
		if pdbname not in opm_data:
			print("{0} not in opm_data".format(pdbname))
			continue
		if 'tmchains' not in opm_data[pdbname]:
			print("{0} does not have tmchains entry in opm_data".format(pdbname))
			continue	
		ntmc = len(list(opm_data[pdbname]['tmchains']))
		nc = len(opm_data[pdbname]['FROM_PDB']['BIOMATRIX'])*len(opm_data[pdbname]['FROM_PDB']['BIOLOGICAL_UNIT'])
		if ntmc > nc:
			print("Struct {0} has more TM chains ({1})than chains ({2})! (...?)".format(pdbname, ntmc, nc))
		NTMC_hist.append(ntmc)
		stackhist[stackcounter][nc-1] = ntmc
		FTMC_hist.append(ntmc/nc)
		stackcounter += 1
		if ntmc not in NTMC_dict:
			NTMC_dict[ntmc] = 0
		NTMC_dict[ntmc] += 1
	NTMCfile = open(NTMCfilename, 'w')
	for n in range(1, max(NTMC_hist)+1):
		if n not in list(NTMC_dict.keys()):
			NTMCfile.write("{0}\t{1}\n".format(n, 0))
		else:
			NTMCfile.write("{0}\t{1}\n".format(n, NTMC_dict[n]))

	fig, ax = plt.subplots()
	bins = np.arange(0.5, max(NTMC_hist) + 1.5, 1)
	weights = np.ones_like(NTMC_hist)/len(NTMC_hist)
	ax.hist(NTMC_hist, bins=bins, weights=weights)
	ax.set_xlabel('Number of transmembrane subunits')
	ax.set_ylabel('Fraction of structures')
	NTMCfig = locations['FSYSPATH']['webarchive'] + 'tmchains_hist.png'
	plt.savefig(NTMCfig)

	fig, ax = plt.subplots()
	bins = np.arange(0.5, max(NTMC_hist) + 1.5, 1)
	weights = np.ones_like(stackhist)
	legend = []
	sumw = 0
	for ic in range(ncmax):
		for j in range(len(stackhist.T[ic])):
			if stackhist.T[ic][j] > 0:
				sumw += 1
		legend.append(str(ic+1))
	weights = np.ones_like(stackhist)/sumw
	colhues = 4
	colors = []
	for i in range(1,colhues+1):
		for ii in range(1,colhues+1):
			for iii in range(1,colhues+1):
				colors.append((i/colhues, ii/colhues, iii/colhues))
	colors = sorted(colors, key = lambda x : (-(x[0]+x[1]+x[2]), x[0], x[1], x[2]))
	print(colors)
	print(stackhist.size)
	x, xx, patches = ax.hist(stackhist, bins=bins, weights=weights, stacked=True, color=colors[ncmax:0:-1], label=legend)
#	jet = pl.get_cmap('jet', len(patches))
#	for i in range(len(patches)):
#		patches[i].set_facecolor(jet(i))
	print("stackhist", stackhist)
	ax.set_xlabel('Number of transmembrane subunits')
	ax.set_ylabel('Fraction of structures')
	NTMCfig = locations['FSYSPATH']['webarchive'] + 'stacked_tmchains_hist.png'
	art = []
	lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=6)
	art.append(lgd)
	plt.savefig(NTMCfig, additional_artists=art, bbox_inches="tight")


	fig, ax = plt.subplots()
	binsize = 0.1
	bins = np.arange(0 - binsize/2, 1 + binsize*1.5, binsize)
	weights = np.ones_like(FTMC_hist)/len(FTMC_hist)
	ax.hist(FTMC_hist, bins=bins, weights=weights)
	ax.set_xlabel('Fraction of transmembrane subunits in complex')
	ax.set_ylabel('Fraction of structures')
	FTMCfig = locations['FSYSPATH']['webarchive'] + 'fraction_tmchains_hist.png'
	plt.savefig(FTMCfig)

	# Number of TM regions

	if 'alpha' in list(table.keys()):
		TMR_hist_alpha = []
		for tt in range(1, max([int(x) for x in list(table['alpha'].keys())])+1):
			if tt in list(table['alpha'].keys()):
				n = len(list(table['alpha'][tt].keys()))
			else:
				n = 0
			for i in range(n):
				TMR_hist_alpha.append(tt)
		if TMR_hist_alpha:
			fig, ax = plt.subplots()
			bins = np.arange(min(TMR_hist_alpha)-0.5, max(TMR_hist_alpha)+1.5)
			weights = np.ones_like(TMR_hist_alpha)/len(TMR_hist_alpha)
			ax.hist(TMR_hist_alpha, bins=bins, weights=weights)
			ax.set_xlabel('Number of transmembrane regions in a subunit')
			ax.set_ylabel('Fraction of structures')
			aTMRfig = locations['FSYSPATH']['webarchive'] + 'alpha_tmregions_hist.png'
			plt.savefig(aTMRfig)

	if 'beta' in list(table.keys()):
		TMR_hist_beta = []
		for tt in range(1, max([int(x) for x in list(table['beta'].keys())])+1):
			if tt in list(table['beta'].keys()):
				n = len(list(table['beta'][tt].keys()))
			else:
				n = 0
			for i in range(n):
				TMR_hist_beta.append(tt)
		if TMR_hist_beta:
			fig, ax = plt.subplots()
			bins = np.arange(min(TMR_hist_beta)-0.5, max(TMR_hist_beta)+1.5)
			weights = np.ones_like(TMR_hist_beta)/len(TMR_hist_beta)
			ax.hist(TMR_hist_beta, bins=bins, weights=weights)
			ax.set_xlabel('Number of transmembrane regions in a subunit')
			ax.set_ylabel('Fraction of structures')
			bTMRfig = locations['FSYSPATH']['webarchive'] + 'beta_tmregions_hist.png'
			plt.savefig(bTMRfig)

	# Coverage of each chain (divided by cl and tt)

	coverage = {}
	for pdbname in list(names.keys()):
		for struct in names[pdbname]:
			ntmr = 0
			chain = struct[5]
			if chain not in opm_data[pdbname]['FROM_PDB']:
				print("Chain {0} not in FROM_PDB of entry {1}".format(chain, pdbname))
				continue
			for resid in opm_data[pdbname]['FROM_PDB'][chain]['RESIDS']:
				for tmr in opm_data[pdbname][chain]['tmdoms']:
					if int(resid) >= int(tmr[0]) and int(resid) <= int(tmr[1]):
						ntmr += 1
						break
			key = neighbors[struct][0]+str(neighbors[struct][1]).zfill(2)
			if key not in list(coverage.keys()):
				coverage[key] = []
			coverage[key].append(ntmr/len(opm_data[pdbname]['FROM_PDB'][chain]['RESIDS']))

	ybins = np.mgrid[0:1:50j]
	nybins = len(ybins)
	H = - np.ones((len(list(coverage.keys())), nybins))
	i = 0
	for cltt in sorted(list(coverage.keys())):
		print(cltt)
		weights = np.ones_like(coverage[cltt])/len(coverage[cltt])
		tmp_h, e = np.histogram(coverage[cltt], bins=ybins, weights=weights)
		H[i] = np.array(list(tmp_h) + [0]*(50-len(tmp_h)))
		i += 1
	xbins = np.arange(i)
	nxbins = i
	print(H)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	covfig = locations['FSYSPATH']['webarchive'] + 'coverage.png'

	X,Y = np.meshgrid(ybins, xbins)
	alphalimit = len([sorted(list(coverage.keys()))[i] for i in xbins if 'alpha' in list(coverage.keys())[i]])
	print(X[:][:alphalimit].shape, Y[:][:alphalimit].shape, H[:][:alphalimit].shape)
	ax.plot_surface(X[:][:alphalimit], Y[:][:alphalimit], H[:][:alphalimit], rstride=1, cstride=1000, color='r', shade=False, lw=.5)
	ax.plot_surface(X[:][alphalimit:], Y[:][alphalimit:], H[:][alphalimit:], rstride=1, cstride=1000, color='cyan', shade=False, lw=.5)
#	ax.set_xlabel("Class and topology")
	if nxbins < 6:
		xt = xbins
		xtl = [sorted(list(coverage.keys()))[i] for i in xbins]
	else:
		xt = np.arange(0, nxbins, 5)
		xtl = [sorted(list(coverage.keys()))[i] for i in xt]
	ax.set_yticks(xt)
	ax.set_yticklabels(xtl, rotation='vertical')
	ax.set_xticks(np.mgrid[0:1:6j])
	ax.set_xlabel("Fraction covered by TM regions")
	ax.set_zlabel("Fraction of structures")

#	ax.set_zlim(0, 5)
#	ax.set_xlim(-51, 51)
#	ax.set_zlabel("Intensity")
	ax.view_init(20,-120)
	plt.savefig(covfig)

	# Length of each TM region

	tmreglen = []
	structured_x = []
	structured_y = []
	for pdbname in list(names.keys()):
		for struct in names[pdbname]:
			chain = struct[5]
			for tmr in opm_data[pdbname][chain]['tmdoms']:
				if abs(int(tmr[0])-int(tmr[1])) > 50:
					print(struct, tmr, abs(int(tmr[0])-int(tmr[1])))
				tmreglen.append(abs(int(tmr[0])-int(tmr[1])) + 1)
				inside = 0
				for nr in range(int(tmr[0]), int(tmr[1])+1):
					for seg in opm_data[pdbname][chain]['segments']:
						if nr >= int(seg[0].strip()) and nr <= int(seg[1].strip()):
							inside += 1
				structured_y.append(inside / (abs(int(tmr[0])-int(tmr[1])) + 1))
	print(structured_y)
	structured_x[:] = tmreglen[:]
			
	fig, ax = plt.subplots()
	bins = np.arange(min(tmreglen)-0.5, max(tmreglen) + 1.5)
	weights = np.ones_like(tmreglen)/len(tmreglen)
	h, bins, patches = ax.hist(tmreglen, bins=bins)
	print(h)
	h, bins, patches = ax.hist(tmreglen, bins=bins, weights=weights)
	ax.set_xlabel('Length of the TM region')
	ax.set_ylabel('Fraction of TM regions')
	TMRLfig = locations['FSYSPATH']['webarchive'] + 'tmrlength_hist.png'
	plt.savefig(TMRLfig)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	xbins = np.arange(min(tmreglen)-0.5, max(tmreglen) + 1.5)
	ybins =  np.arange(0, 1.1, 0.1)

	hist, xedges, yedges = np.histogram2d(structured_x, structured_y, bins=(xbins,ybins), normed=False)
	print(hist)
	for x in range(len(hist)):
		s = 0
		for y in range(len(hist[x])):
			ss = hist[x][y]
			hist[x][y] = ss + s
			s += ss
	print(hist)
	xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])
	xpos = xpos.flatten()/2.
	ypos = ypos.flatten()/2.
	zpos = np.zeros_like (xpos)
	colors = (['b']*(len(hist)))*(len(hist[0])-1) + (['r']*(len(hist)))
	
	dx = xedges [1] - xedges [0]
	dy = (yedges [1] - yedges [0])/3
	dz = hist.T.flatten()

	ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors)
#	bins = np.arange(min(tmreglen)-0.5, max(tmreglen) + 1.5)
#	weights = np.ones_like(tmreglen)/len(tmreglen)
#	ax.hist(tmreglen, bins=bins, weights=weights, zdir='y')
	TMRLfig = locations['FSYSPATH']['webarchive'] + 'structured_tmrlength_hist.png'
	plt.savefig(TMRLfig)

	
	return neighbors, names, opm_data, table, equivalences

def cdist(c1, c2):
	return ((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)**0.5

def find_disposition(tm, pcoord, placed, j):
#	print(tm, j)
	n = len(placed)

	# Transform polar coordinates in cartesian coordinates for calculations
	ccoord = np.copy(pcoord)
	for a in range(pcoord.shape[0]):
		ccoord[a][0] = pcoord[a][0]*math.cos(pcoord[a][1])
		ccoord[a][1] = pcoord[a][0]*math.sin(pcoord[a][1])

	sampled_pcoord = np.zeros((360, 2))
	sampled_ccoord = np.zeros((360, 2))
	for a in range(sampled_pcoord.shape[0]):
		theta = 2*math.pi*a/sampled_pcoord.shape[0]
#		print(theta)
		sampled_pcoord[a] = (tm[0][j], theta)
		sampled_ccoord[a] = (tm[0][j] * math.cos(theta), tm[0][j] * math.sin(theta))
#		print(sampled_ccoord[a])

	# Sum the potential
	V = np.zeros(360)
	for a in range(V.shape[0]):
		for i in placed:
#			print(a, i, ccoord[i], sampled_ccoord[a])
			d = cdist(ccoord[i], sampled_ccoord[a])
			d0 = tm[i][j]
			V[a] += (d - d0)**2
#		print(sampled_ccoord[a][0], sampled_ccoord[a][1], V[a])

	am = V.argmin()
	return sampled_pcoord[am], sampled_ccoord[am]

def polar(options, locations, table):
	this_name = 'polar'

	for cl in sorted(list(table.keys())):
		for tt in sorted(list(table[cl].keys())):
			allstruct= set()
			for s1 in sorted(list(table[cl][tt].keys())):
				allstruct.add(s1)
				for s2 in sorted(list(table[cl][tt][s1].keys())):
					allstruct.add(s2)
			# N is the total number of structures
			N = len(allstruct)

			print_log(this_name, "Found {0} structures in set {1} {2}".format(N, cl, tt))

			# all_list is a sorted list of structure names
			all_list = sorted(list(allstruct))

			# s_to_i is the table converting names to indices
			s_to_i = {}
			i_to_s = {}
			for i, s in enumerate(all_list):
				s_to_i[s] = i
				i_to_s[i] = s

			# Extract TMscore and sequence identity AS METRICS
			tm = - np.ones((N, N))
			SI = - np.ones((N, N))
			for s1 in sorted(list(table[cl][tt].keys())):
				i1 = s_to_i[s1]
				for s2 in sorted(list(table[cl][tt][s1].keys())):
					i2 = s_to_i[s2]
					tm[i1][i2] = 1.0 - table[cl][tt][s1][s2][2]
					SI[i1][i2] = 1.0 - table[cl][tt][s1][s2][0]
			
			tm = tm - np.diag(np.diag(tm))
			SI = SI - np.diag(np.diag(SI))

			# Thresholds for TMscore and sequence identity AS METRICS
			tm_thr = 1.0 - float(options['cluster_thr'])
			SI_thr = 1.0 - float(options['subunit_thr'])

			# Create main list for each struct
			# Contains ref numbers
			I_list = {}
			for I in range(N):
				I_list[I] = set()
				for i in range(N):
					if I != i and tm[I][i] < tm_thr:
						I_list[I].add(i)

			# Main loop
			for I in range(N):
				print("Structure {0} has {1} neighbors:".format(i_to_s[I], len(I_list[I])))
				sys.stdout.flush()

				# This is the sorted list for this subproblem
				sub_list = [I]
				for i in sorted(I_list[I]):
					sub_list.append(i)

				# This is the conversion table for this subproblem
				i_to_j = {I : 0}
				j_to_i = {0 : I}
				for j, i in enumerate(sub_list):
					i_to_j[i] = j
					j_to_i[j] = i

				# TMscore and SI matrices for the subproblem
				n = len(I_list[I])
				subtm = np.zeros((n+1,n+1))
				subSI = np.zeros((n+1,n+1))

				# Transforms global reference numbers to subproblem reference numbers
				# (without holes, so that a matrix can be compiled)
				for i in sub_list:
					j = i_to_j[i]
					for ii in sub_list:
						jj = i_to_j[ii]
						subtm[j][jj] = tm[i][ii]
						subSI[j][jj] = SI[i][ii]

				# Polar coordinates
				pcoords = - np.ones((n+1,2))
				ccoords = - np.ones((n+1,2))

				# Coords of I
				pcoords[0][0] = 0.0
				pcoords[0][1] = 0.0

				ccoords[0][0] = 0.0
				ccoords[0][1] = 0.0

				build_matrix = 10*np.ones((n+1,n+1))
				build_matrix[0] = np.copy(subtm[0])
				placed = set()
				placed0 = set([0])
				while len(placed) + len(placed0) < n+1:
					linear_form = []
					for nn in range(build_matrix.shape[0]):
						for mm in range(build_matrix.shape[1]):
							linear_form.append((nn, mm, build_matrix[nn][mm]))
					linear_form = sorted(linear_form, key = lambda x : (x[2], x[0]))
					for tr in linear_form:
						if j_to_i[tr[1]] not in I_list[I]:
							continue
						if tr[0] == 0 and tr[2] == 0 and tr[1] not in placed0:
							placed0.add(tr[1])
							pcoords[tr[1]][0] = 0.0
							pcoords[tr[1]][1] = 0.0
							ccoords[tr[1]][0] = 0.0
							ccoords[tr[1]][1] = 0.0
							print("{0} placed!".format(j_to_i[tr[1]]))
							sys.stdout.flush()
							continue
						elif tr[1] not in placed and tr[1] not in placed0:
							pcoords[tr[1]], ccoords[tr[1]] = find_disposition(build_matrix, pcoords, placed, tr[1])
							build_matrix[tr[1]] = np.copy(subtm[tr[1]])
							placed.add(tr[1])
							print("{0} placed!".format(j_to_i[tr[1]]))
							sys.stdout.flush()
							break

				# Writing coordinates file and 
				coordfilename = locations['FSYSPATH']['polar'] + '{0}_neigh.pcoord'.format(i_to_s[I])
				coordfile = open(coordfilename, 'w')
				coordfile.write("{0} 0.0 0.0\n".format(i_to_s[I]))
				theta = []
				rho = []
				for i in I_list[I]:
					s = i_to_s[i]
					r = (random.random()-0.5)*2*math.pi/360
					theta.append(pcoords[i_to_j[i]][1]+r)
					rho.append(pcoords[i_to_j[i]][0])
					coordfile.write("{0} {1} {2}\n".format(s, pcoords[i_to_j[i]][0], pcoords[i_to_j[i]][1]+r))
				coordfile.close()

				if not theta:
					continue

				# Creating polar plot
				fig = plt.figure(1)
				ax = plt.subplot(111, projection='polar')
				ax.scatter(theta, rho)
				ax.set_rmax(0.4)
				ax.axes.get_xaxis().set_ticks([])
				plt.yticks([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4], ['', '0.9', '', '0.8', '', '0.7', '', '0.6'])
				figname = locations['FSYSPATH']['polar'] + 'p_' + i_to_s[I] + ".png"
				plt.savefig(figname)
				plt.close(fig)


def scatterplots(options, locations, neighbors, table):
	# Read table and compile histogram (200x200 bins)
	tablefilename = locations['FSYSPATH']['main'] + options['output_tab']
	f = open(tablefilename, 'r')
	text = f.read().split('\n')
	f.close()
	Dx = []
	Dy = []
	for cl in sorted(list(table.keys())):
		for tt in sorted(list(table[cl].keys())):
			allstruct= set()
			for s1 in sorted(list(table[cl][tt].keys())):
				if s1 not in list(neighbors.keys()):
					continue
				for s2 in sorted(list(table[cl][tt][s1].keys())):
					if s2 not in list(neighbors.keys()):
						continue
					Dx.append(float(table[cl][tt][s1][s2][0]))
					Dy.append(float(table[cl][tt][s1][s2][2]))
	Dx = np.array(Dx)
	Dy = np.array(Dy)
	xmin = Dx.min()
	xmax = Dx.max()
	ymin = Dy.min()
	ymax = Dy.max()
	xedges = np.mgrid[0:1:200j]
	yedges = np.mgrid[0:1:200j]
	H, xedges, yedges = np.histogram2d(Dy, Dx, bins=(xedges, yedges))

	# Take log of histogram, adjust -inf to 0
	H = np.log(H)
	H[H == -np.inf] = 0

	# Build KDE on log-histogram, by creating a fake set of data
	DDx = []
	DDy = []
	for x in range(len(H)):
		for y in range(len(H[x])):
			for i in range(math.ceil(H[x][y])):
				DDx.append(x/len(H))
				DDy.append(y/len(H[x]))
				DDx.append(2-(x/len(H)))
				DDy.append(y/len(H[x]))
				DDx.append(x/len(H))
				DDy.append(2-(y/len(H[x])))
				DDx.append(2-(x/len(H)))
				DDy.append(2-(y/len(H[x])))
	DDx = np.array(DDx)
	DDy = np.array(DDy)
	xmin = DDx.min()
	xmax = DDx.max()
	ymin = DDy.min()
	ymax = DDy.max()
	X, Y = np.mgrid[0:2:400j, 0:2:400j]
#	X, Y = np.mgrid[0.05:0.95:200j, 0.05:0.95:200j]
	positions = np.vstack([X.ravel(), Y.ravel()])
	values = np.vstack([DDx, DDy])
	kernel = stats.gaussian_kde(values)
	Z = np.reshape(kernel(positions), X.shape)
#	ar = (xmax-xmin)/(ymax-ymin)
	Z = Z.T

	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02
	binwidth = 0.05

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	for struct in sorted(list(neighbors.keys())):
		figfilename = locations['FSYSPATH']['densityscatter'] + 'ds_' + struct + '.png'
#		if os.path.exists(figfilename):
#			continue
		print(struct)
		ZZ = np.array(Z)
		fig = plt.figure(1, figsize=(8, 8))
		ax = plt.axes(rect_scatter)
		axHistx = plt.axes(rect_histx)
		axHisty = plt.axes(rect_histy)
#		ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,extent=[0, 1, 0, 1], aspect=1)
		ax.set_xlim([0, 1])
		ax.set_ylim([0, 1])
		levels = np.array([0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8])
		cfset = ax.contourf(X, Y, ZZ, levels=levels, cmap="Blues")
		cset = ax.contour(X, Y, ZZ, levels=levels, colors='k')
		plt.locator_params(nbins=6)
		ax.imshow(np.rot90(ZZ), cmap="Blues",extent=[0, 2, 0, 2], aspect=1)
		XX = []
		YY = []
		cl = neighbors[struct][0]
		tt = neighbors[struct][1]	
		for s2 in sorted(list(table[cl][tt][struct].keys())):
			XX.append(float(table[cl][tt][struct][s2][0]))
			YY.append(float(table[cl][tt][struct][s2][2]))
		if not XX:
			plt.close(fig)
			continue
		ax.scatter(XX, YY, color='pink')
		ax.set_xlabel("Sequence Identity")
		ax.set_ylabel("TM-Score")
		
		bins = np.arange(0, 1 + binwidth, binwidth)
		plt.locator_params(nbins=4)
		axHistx.hist(XX, bins=bins, color='pink')
		axHistx.set_xlim(ax.get_xlim())
		axHistx.set_xticks([])

		plt.locator_params(nbins=4)
		axHisty.hist(YY, bins=bins, orientation='horizontal', color='pink')
		axHisty.set_ylim(ax.get_ylim())
		axHisty.set_yticks([])

		fig.savefig(figfilename)
		plt.close(fig)
#		break

def gif_instr_creator(options, locations, neighbors, names):
	this_name = 'gif_creator'
	tic = time.time()
	script_folder = os.getcwd() + '/'
	fnull = open(os.devnull, 'w')	

	pymolpath = options['pymol_path']

	shutil.copy(script_folder + 'create_whole_gif.cmd', locations['FSYSPATH']['webarchive'])
	shutil.copy(script_folder + 'create_chain_gif.cmd', locations['FSYSPATH']['webarchive'])
	shutil.copy(script_folder + 'create_whole_png.cmd', locations['FSYSPATH']['webarchive'])
	shutil.copy(script_folder + 'create_chain_png.cmd', locations['FSYSPATH']['webarchive'])

	for pdbname in sorted(list(names.keys())):
		wholecreatorfilename = locations['FSYSPATH']['wholegif'] + pdbname + '_gif_instr.sh'
		wholecreatorfile = open(wholecreatorfilename, 'w')
		wholecreatorfile.write("mkdir {0}_pngs/\n".format(pdbname))
		wholecreatorfile.write("cd {0}_pngs/\n".format(pdbname))
		wholecreatorfile.write("cp ../../../{0} .\n".format(locations['FSYS']['whole'] + pdbname + '_opm.pdb'))
		wholecreatorfile.write("sed 's/STRUCT/{0}/g' ../../create_whole_gif.cmd > create_whole_gif.cmd\n".format(pdbname + '_opm.pdb'))
		wholecreatorfile.write("pymol -ucq create_whole_gif.cmd\n")
		nframe = 1
		wholecreatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
		for nframe in range(2,51):
			wholecreatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
			wholecreatorfile.write("convert {0}.png {1}.gif\n".format('f'+str(nframe).zfill(4), 'f'+str(102-nframe).zfill(4)))
		nframe = 51
		wholecreatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
		wholecreatorfile.write("convert -loop 0 -delay 3 -type palette f*.gif {0}.gif\n".format(pdbname))
		wholecreatorfile.write("mv {0}.gif ../\n".format(pdbname))
		wholecreatorfile.write("cd ../\n")
		wholecreatorfile.write("rm -rf {0}_pngs/\n".format(pdbname))
		wholecreatorfile.close()
			
		os.system("cp create_whole_png.cmd {0}; cd {0}; cp {1} .; sed -i 's/STRUCT/{3}/g' create_whole_png.cmd; {5} -ucq create_whole_png.cmd; mv {2}.png {4}.png".format(locations['FSYSPATH']['wholegif'], locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb', 'f'+str(1).zfill(4), pdbname + '_opm.pdb', pdbname, pymolpath))


	for struct in sorted(list(neighbors.keys())):
		pdbname = struct[:4]
		creatorfilename = locations['FSYSPATH']['chainsgif'] + struct + '_gif_instr.sh'
		creatorfile = open(creatorfilename, 'w')
		creatorfile.write("mkdir {0}_pngs/\n".format(struct))
		creatorfile.write("cd {0}_pngs/\n".format(struct))
		creatorfile.write("cp ../../../{0} .\n".format(locations['FSYS']['whole'] + pdbname + '_opm.pdb'))
		creatorfile.write("sed 's/STRUCT/{0}/g' ../../create_chain_gif.cmd > create_chain_gif.cmd\n".format(pdbname + '_opm.pdb'))
		creatorfile.write("sed -i 's/XXX/{0}/g' create_chain_gif.cmd\n".format(struct[5]))
		creatorfile.write("pymol -ucq create_chain_gif.cmd\n")
		nframe = 1
		creatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
		for nframe in range(2,51):
			creatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
			creatorfile.write("convert {0}.png {1}.gif\n".format('f'+str(nframe).zfill(4), 'f'+str(102-nframe).zfill(4)))
		nframe = 51
		creatorfile.write("convert {0}.png {0}.gif\n".format('f'+str(nframe).zfill(4)))
		creatorfile.write("convert -loop 0 -delay 3 -type palette f*.gif {0}.gif\n".format(struct))
		creatorfile.write("mv {0}.gif ../\n".format(struct))
		creatorfile.write("cd ../\n")
		creatorfile.write("rm -rf {0}_pngs/\n".format(struct))

		os.system("cp create_chain_png.cmd {0}; cd {0}; cp {1} .; sed -i 's/STRUCT/{3}/g' create_chain_png.cmd; sed -i 's/XXX/{5}/g' create_chain_png.cmd; {6} -ucq create_chain_png.cmd; mv {2}.png {4}.png".format(locations['FSYSPATH']['chainsgif'], locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb', 'f'+str(1).zfill(4), pdbname + '_opm.pdb', struct, struct[5], pymolpath))
	
	toc = time.time()	
	print_log(this_name, "Time spent: {0}".format(toc-tic))

def transfer_estimator(locations, neighbors, opm_data, equivalences):
	this_name = 'transfer_estimator'
	tic = time.time()
	def k_func(d):
		return (1 - d**2)/(1 - d**4)

	if not opm_data:
		opm_data = read_data(locations['SYSFILES']['opmdata'])

	# Create topology groups for first loop
	groups = {'alpha' : {}, 'beta' : {}}
	for struct in list(neighbors.keys()):
		cl, tt, se_nn, st_nn, tot_nn = neighbors[struct]
		tt = int(tt)
		if tt not in groups[cl]:
			groups[cl][tt] = set()
		groups[cl][tt].add(struct)

	print_log(this_name, "Creating distance dictionary")
	for cl in sorted(list(groups.keys())):
		for tt in sorted(list(groups[cl].keys())):
			if len(groups[cl][tt]) < 2:
				continue
			for struct in sorted(groups[cl][tt]):
				is_equivalent = False
				if struct in list(equivalences.keys()):
					is_equivalent = True
#				if struct != '1bcc_N':
#					continue
				figfilename = locations['FSYSPATH']['distributions'] + 'distr_' + struct + '.png'
#				if os.path.exists(figfilename):
#					continue
				print_log(this_name, "Creating distance dictionary for {0}".format(struct))
				sys.stdout.flush()
				anc = {}
				eqstruct = struct
				if is_equivalent:
					eqstruct = equivalences[struct]
				strpdbfilename = locations['FSYSPATH'][cl] + str(tt) + '/' + locations['TREE']['straln'] + 'str_' + eqstruct + '_pdb.dat'
				print(strpdbfilename)
				strpdbfile = open(strpdbfilename, 'r')
				text = strpdbfile.read().split('\n')
				strpdbfile.close()
				jump_this_entry = False
				for line in text:
					if not line:
						continue
					fields = line.split()
					if 'BEGIN' in line:
						s1 = struct
						s2 = fields[4]
						if s1 == s2:
							jump_this_entry = True
						else:
							jump_this_entry = False
						if s1 not in anc:
							anc[s1] = {}
						anc[s1][s2] = []
						ter = False
					elif 'ATOM' in line and not jump_this_entry:
						if not ter:
							nres = int(fields[4])
							cx, cy, cz = float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:52].strip())
							anc[s1][s2].append((nres, cx, cy, cz))
							#print('FIRST', nres)
							#print(anc[s1][s2][-1])
						else:
							nres = anc[s1][s2][cr][0]
							dist = ((anc[s1][s2][cr][1] - float(line[30:38].strip()))**2 + (anc[s1][s2][cr][2] - float(line[38:46].strip()))**2 + (anc[s1][s2][cr][3] - float(line[46:52].strip()))**2)**0.5
							anc[s1][s2][cr] = (nres, dist)
							#print(cr, 'MATCH IT', nres)
							#print(anc[s1][s2][cr])
							cr += 1
					elif 'TER' in line and not jump_this_entry:
						ter = True
						cr = 0

				if is_equivalent and struct in list(anc.keys()):
					if eqstruct not in list(anc[struct].keys()):
						anc[struct][eqstruct] = []
					for nres in opm_data[struct[:4]]['FROM_PDB'][struct[5]]['RESIDS']:
						anc[struct][eqstruct].append((int(nres), 0.0))

				print_log(this_name, "Creating distributions and figures")
				sys.stdout.flush()
				S = struct
				figfilename = locations['FSYSPATH']['distributions'] + 'distr_' + S + '.png'
				distrfilename = locations['FSYSPATH']['distributions'] + 'distr_' + S + '.dat'
				distrfile = open(distrfilename, 'w')
				DYY = {}
				ymax = 0
				for s1 in sorted(groups[cl][tt]):
					if s1 == S:
						continue
					distrfile.write('BEGIN CHAIN1: {0} CHAIN2: {1}\n'.format(S, s1))
					for nres, dist in sorted(anc[S][s1], key = lambda x : x[0]):
						delta = dist
						if nres not in DYY:
							DYY[nres] = []
						DYY[nres].append(delta)
						if delta > ymax:
							ymax = delta
						distrfile.write("{0:10d} {1:20.6f}\n".format(nres, delta))
					distrfile.write('\nEND\n\n')
				distrfile.close()
				if ymax == 0:
					print("WARNING: {0} not generated because is compared only with identical structures!".format(struct))
					continue
				if not DYY:
					print("WARNING: {0} not generated because DYY is void!".format(struct))
					continue
				xmin = sorted(list(DYY.keys()))[0]
				xmax = sorted(list(DYY.keys()))[-1]
				M = np.zeros((1+xmax, 100))*(-math.inf)
				for nres in range(xmin, xmax+1):
					if not nres in DYY:
						DYY[nres] = [0]*len(groups[cl][tt])
					elif len(DYY[nres]) < len(groups[cl][tt]):
						for nn in range(len(groups[cl][tt]) - len(DYY[nres])):
							DYY[nres].append(0)
					DYY[nres] = np.array(sorted(DYY[nres]))
					try:
						kernel = stats.gaussian_kde(DYY[nres])
						Y = np.mgrid[0:ymax:100j]
						M[nres] = kernel(Y)
						maxline = M[nres].max() # THIS IS NOT FORMALLY CORRECT, BUT ENHANCES THE VIEW WITH THE CHOSEN COLOR SCALE
						M[nres] = M[nres]/maxline
					except:
						M[nres] = np.zeros(100)
				fig, ax = plt.subplots()
				ar = (xmax-xmin)/(ymax-0)
				ax.set_xlim([xmin, xmax])
				ax.set_ylim([0, ymax])
#				ax.imshow(np.rot90(M), cmap=plt.cm.gist_earth_r, extent=[0, xmax, 0, ymax], aspect=ar/2)
				ax.imshow(np.rot90(M), cmap="Blues", extent=[0, xmax, 0, ymax], aspect=ar/2)
				ax.set_xlabel('Residue Number')
				ax.set_ylabel('Distance')
				ax.set_title('Density of distances of each residue from the \ncorresponding ones in the set of alignments\n ')
				if struct[:4] in opm_data and struct[5] in opm_data[struct[:4]] and 'tmdoms' in opm_data[struct[:4]][struct[5]]: # and struct[5] in opm_data[struct[:4]]['FROM_PDB']:
					for tmdom in opm_data[struct[:4]][struct[5]]['tmdoms']:
						a = int(tmdom[0].strip())
						b = int(tmdom[1].strip())
#						a = opm_data[struct[:4]]['FROM_PDB'][struct[5]]['RESIDS'].index(tmdom[0].strip())
#						b = opm_data[struct[:4]]['FROM_PDB'][struct[5]]['RESIDS'].index(tmdom[1].strip())
						plt.axvspan(a, b, facecolor='0.1', alpha=0.2)
				else:
					print("shading failed for struct {0}".format(struct))
					if struct[:4] not in opm_data:
						print("{0} not in opm_data".format(struct[:4]))
					else:
						if struct[5] not in opm_data[struct[:4]]:
							print("Chain {0} not in opm_data[{1}]".format(struct[5],struct[:4]))
						else:
							if 'tmdoms' not in opm_data[struct[:4]][struct[5]]:
								print("tmdoms not found in opm_data[{1}][{0}]".format(struct[5],struct[:4]))
					#	if struct[5] not in opm_data[struct[:4]]['FROM_PDB']:
					#		print("Chain {0} not in opm_data[{1}]['FROM_PDB']".format(struct[5],struct[:4]))
				plt.savefig(figfilename)
				plt.close(fig)
			

def webcontent(locations, neighbors, names, opm_data):
#	names = {}
#	for struct in sorted(list(neighbors.keys())):
#		pdbname = struct[:4]
#		if pdbname not in list(names.keys()):
#			names[pdbname] = set()
#		names[pdbname].add(struct[5])

	if not opm_data:
		opm_data = read_data(locations['SYSFILES']['opmdata'])
			
	webarchive_filename = locations['SYSFILES']['webwhole'] 
	webarchive_file = open(webarchive_filename, 'w')
	webarchive_file.write('<Database>\n')
	for pdbname in sorted(names.keys()):
		print(pdbname)
		wholetmp_filename = locations['FSYSPATH']['webarchive'] + 'temp/' + pdbname + '_whole.xml'
		cesymm_text = ''
		symd_text = ''
		analysis_text = ''
		if os.path.exists(wholetmp_filename):
			wholetmp_file = open(wholetmp_filename, 'r')
			text = wholetmp_file.read().split('\n')
			wholetmp_file.close()
			cesymm = False
			symd = False
			for line in text:
				if not line:
					continue
				if cesymm:
					cesymm_text += '\t'+line+'\n'
				elif symd:
					symd_text += '\t'+line+'\n'
				elif '<Size' in line:
					sizevalue = line.split()[1]
				if '<CE-Symm>' in line:
					cesymm = True
				elif '</CE-Symm>' in line:
					cesymm = False
				elif '<SymD>' in line:
					symd = True
				elif '</SymD>' in line:
					symd = False
		else:
			print("Not available in Toni's files: {0}".format(pdbname))
			continue
		wholeantmp_filename = locations['FSYSPATH']['webarchive'] + 'temp/analysis/' + pdbname + '_analysis_whole.xml'
		if os.path.exists(wholeantmp_filename):
			wholeantmp_file = open(wholeantmp_filename, 'r')
			text = wholeantmp_file.read().split('\n')
			wholeantmp_file.close()
			analysis = False
			analysis_text = ''
			for line in text:
				if not line:
					continue
				if '</Analysis>' in line:
					analysis = False
				if analysis:
					analysis_text += line+'\n'
				if '<Analysis>' in line:
					analysis = True
			
		webarchive_file.write('\t<Structure ID="{0}">\n'.format(pdbname))
		# HEADER
		webarchive_file.write('\t\t<Header>\n')
		if opm_data[pdbname]['is_representative']:
			webarchive_file.write('\t\t\t<Title> {0} </Title>\n'.format(opm_data[pdbname]['title']))
		elif opm_data[pdbname]['related_to'] in list(opm_data.keys()):
			webarchive_file.write('\t\t\t<Title> {0} </Title>\n'.format(opm_data[opm_data[pdbname]['related_to']]['title']))
			opm_data[pdbname]['title'] = opm_data[opm_data[pdbname]['related_to']]['title']
		else:
			info_filename = locations['FSYSPATH']['OPMhtmls'] + opm_data[pdbname]['related_to'] + '_info.html'
			info_file = open(info_filename, 'r')	
			info_text = info_file.read().split('\n')
			info_file.close()
			for info_line in info_text:
				if not info_line:
					continue
				if '<!-- body section -->' in line:
					activate_namesearch = True
				if struct in line and activate_namesearch:
					title = [x.strip() for x in re.findall(r'&raquo;(.*?)</h1>$', line, re.DOTALL)]
					webarchive_file.write('\t\t\t<Title> {0} </Title>\n'.format(title))
					opm_data[pdbname]['title'] = title
					break
		webarchive_file.write('\t\t</Header>\n')
		# GENERAL
		webarchive_file.write('\t\t<General>\n')
		webarchive_file.write('\t\t\t<PDBCode> {0} </PDBCode>\n'.format(pdbname))
		if cesymm_text:
			webarchive_file.write('\t\t\t<Size {0} />\n'.format(sizevalue))
		OPMpdb_dict = PDB_parser(locations['FSYSPATH']['whole'], pdbname + '_opm')
		webarchive_file.write('\t\t\t<Chains> {0} </Chains>\n'.format(len(OPMpdb_dict['CHAINS'])))
		webarchive_file.write('\t\t\t<TMChains> {0} </TMChains>\n'.format(len(opm_data[pdbname]['tmchains'])))
		for struct in sorted(names[pdbname]): #opm_data[struct]['tmchains']:
			chain = struct[5]
			webarchive_file.write('\t\t\t<Chain ID="{0}">\n'.format(chain))
			webarchive_file.write('\t\t\t\t<NTMDomains> {0} </NTMDomains>\n'.format(opm_data[pdbname][chain]['ntm']))
			webarchive_file.write('\t\t\t\t<Sequence name="{0}">\n'.format('>'+ struct +':'))
			webarchive_file.write('\t\t\t\t\t<Seq> {0} </Seq>\n'.format(opm_data[pdbname]['FASTA'][chain]))
			webarchive_file.write('\t\t\t\t</Sequence>\n')
			nid = 0
			for segment in opm_data[pdbname][chain]['segments']:
				webarchive_file.write('\t\t\t\t<TMDomain ID="{0}">\n'.format(nid))
				webarchive_file.write('\t\t\t\t\t<TMDRange> {0} - {1} </TMDRange>\n'.format(segment[0], segment[1]))
				webarchive_file.write('\t\t\t\t</TMDomain>\n'.format(nid))
				nid += 1
			webarchive_file.write('\t\t\t</Chain>\n')
		webarchive_file.write('\t\t\t<GeneralStructure_gif>\n')
		webarchive_file.write('\t\t\t\t<Path> {0} </Path>\n'.format(locations['FSYS']['wholegif'] + pdbname + '.gif'))
		webarchive_file.write('\t\t\t\t<DownloadPDBFile> {0} </DownloadPDBFile>\n'.format(locations['FSYS']['whole'] + pdbname + '_opm.pdb'))
		webarchive_file.write('\t\t\t</GeneralStructure_gif>\n')
		webarchive_file.write('\t\t\t<PDB_URL> {0} </PDB_URL>\n'.format(locations['FIXED']['pdbmainsemiurl']+pdbname))
		if opm_data[pdbname]['is_representative']:
			webarchive_file.write('\t\t\t<OPM_URL> {0} </OPM_URL>\n'.format(locations['FIXED']['opmrefsemiurl']+pdbname))
		else:
			webarchive_file.write('\t\t\t<OPM_URL> {0} </OPM_URL>\n'.format(locations['FIXED']['opmrelsemiurl']+pdbname))
		webarchive_file.write('\t\t</General>\n')
		# STRUCTCLASS
		webarchive_file.write('\t\t<StructureInformation>\n')
		for struct in sorted(names[pdbname]): #opm_data[struct]['tmchains']:
			webarchive_file.write('\t\t\t<ChainInformation ID="{0}">\n'.format(struct))
			webarchive_file.write('\t\t\t\t<Member> {0} </Member>\n'.format(struct))
			webarchive_file.write('\t\t\t\t<Class> {0} </Class>\n'.format(neighbors[struct][0]))
			webarchive_file.write('\t\t\t\t<TMdomains> {0} </TMdomains>\n'.format(neighbors[struct][1]))
			webarchive_file.write('\t\t\t\t<SequenceNeighbors> {0} </SequenceNeighbors>\n'.format(neighbors[struct][2]))
			webarchive_file.write('\t\t\t\t<StructureNeighbors> {0} </StructureNeighbors>\n'.format(neighbors[struct][3]))
			webarchive_file.write('\t\t\t\t<TotalNeighbors> {0} </TotalNeighbors>\n'.format(neighbors[struct][4]))
			webarchive_file.write('\t\t\t</ChainInformation>\n')
		webarchive_file.write('\t\t</StructureInformation>\n')
		webarchive_file.write('\t\t<Symmetry>\n')
		webarchive_file.write('\t\t\t<CE-Symm Algorithm="jCE-symm" Version="2.2" >\n')
		if cesymm_text:
			webarchive_file.write(cesymm_text)
		else:
			webarchive_file.write('Not available for this entry\n')
#		webarchive_file.write('\t\t\t</CE-Symm>')
		webarchive_file.write('\t\t\t<SymD Algorithm="SymD" Version="1.6" >\n')
		if symd_text:
			webarchive_file.write(symd_text)
		else:
			webarchive_file.write('Not available for this entry\n')
#		webarchive_file.write('\t\t\t</SymD>')
		webarchive_file.write('\t\t</Symmetry>\n')
		webarchive_file.write('\t\t<SymmetryAnalysis>\n')
		if analysis_text:
			webarchive_file.write(analysis_text)
		else:
			webarchive_file.write('Not available for this entry\n')
		webarchive_file.write('\t\t</SymmetryAnalysis>\n')
		webarchive_file.write('\t</Structure>\n')
	webarchive_file.write('</Database>\n')
	webarchive_file.close()
	
	webarchive_filename = locations['SYSFILES']['webchain']
	webarchive_file = open(webarchive_filename, 'w')
	webarchive_file.write('<Database>\n')
	for pdbname in sorted(list(names.keys())):
		for struct in sorted(names[pdbname]): #opm_data[struct]['tmchains']:
			chain = struct[5]
			chainstmp_filename = locations['FSYSPATH']['webarchive'] + 'temp/' + pdbname + '_chains.xml'
			cesymm_text = ''
			symd_text = ''
			analysis_text = ''
			transfer_text = ''
			if os.path.exists(chainstmp_filename):
				chainstmp_file = open(chainstmp_filename, 'r')
				text = chainstmp_file.read().split('\n')
				chainstmp_file.close()
				cesymm = False
				symd = False
				read_on = False
				cesymm_text = ''
				symd_text = ''
				for line in text:
					if not line:
						continue
					if struct in line:
						read_on = True
					elif read_on and '</Structure>' in line:
						read_on = False
					if not read_on:
						continue
					if cesymm:
						cesymm_text += '\t'+line+'\n'
					elif symd:
						symd_text += '\t'+line+'\n'
					elif '<Size' in line:
						sizevalue = line.split()[1]
					if '<CE-Symm>' in line:
						cesymm = True
					elif '</CE-Symm>' in line:
						cesymm = False
					elif '<SymD>' in line:
						symd = True
					elif '</SymD>' in line:
						symd = False
			else:
				print("Not available in Toni's files: {0}".format(struct))
				continue
			chainantmp_filename = locations['FSYSPATH']['webarchive'] + 'temp/analysis/' + pdbname + '_analysis_chains.xml'
			if os.path.exists(chainantmp_filename):
				chainantmp_file = open(chainantmp_filename, 'r')
				text = chainantmp_file.read().split('\n')
				chainantmp_file.close()
				analysis = False
				read_on = False
				analysis_text = ''
				for line in text:
					if not line:
						continue
					if struct in line:
						read_on = True
					elif read_on and '</Structure>' in line:
						read_on = False
					if not read_on:
						continue
					if '</Analysis>' in line:
						analysis = False
					if analysis:
						analysis_text += line+'\n'
					if '<Analysis>' in line:
						analysis = True
			chaintrtmp_filename = locations['FSYSPATH']['webarchive'] + 'temp/transfer/' + struct + '_transfer.xml'
			if os.path.exists(chaintrtmp_filename):
				chaintrtmp_file = open(chaintrtmp_filename, 'r')
				text = chaintrtmp_file.read().split('\n')
				chaintrtmp_file.close()
				transfer = False
				read_on = False
				transfer_text = ''
				for line in text:
					if not line:
						continue
					if struct in line:
						read_on = True
					elif read_on and '</Structure>' in line:
						read_on = False
					if not read_on:
						continue
					if '</Transfer>' in line:
						transfer = False
					if transfer:
						transfer_text += line+'\n'
					if '<Transfer>' in line:
						transfer = True
			topology_path = locations['FSYSPATH'][neighbors[struct][0]] + str(neighbors[struct][1]) + '/'
			strpdb_filename = topology_path + locations['TREE']['straln'] + 'str_' + struct + '_pdb.dat'

			webarchive_file.write('\t<Structure ID="{0}">\n'.format(struct))
			# HEADER
			webarchive_file.write('\t\t<Header>\n')
			webarchive_file.write('\t\t\t<Title> {0} - chain {1} </Title>\n'.format(opm_data[pdbname]['title'], chain))
			webarchive_file.write('\t\t</Header>\n')
			# GENERAL
			webarchive_file.write('\t\t<General>\n')
			webarchive_file.write('\t\t\t<PDBCode> {0} </PDBCode>\n'.format(struct))
			if cesymm_text:
				webarchive_file.write('\t\t\t<Size {0} />\n'.format(sizevalue))
			webarchive_file.write('\t\t\t<TMDomains> {0} </TMDomains>\n'.format(opm_data[pdbname][chain]['ntm']))
			webarchive_file.write('\t\t\t<Sequence name="{0}">\n'.format('>' + struct + ':'))
			webarchive_file.write('\t\t\t\t<Seq> {0} </Seq>\n'.format(opm_data[pdbname]['FASTA'][chain]))
			webarchive_file.write('\t\t\t</Sequence>\n')
			nid = 0
			for segment in opm_data[pdbname][chain]['segments']:
				webarchive_file.write('\t\t\t<TMDomain ID="{0}">\n'.format(nid))
				webarchive_file.write('\t\t\t\t<TMDRange> {0} - {1} </TMDRange>\n'.format(segment[0], segment[1]))
				webarchive_file.write('\t\t\t</TMDomain>\n'.format(nid))
				nid += 1
			webarchive_file.write('\t\t\t<GeneralStructure_gif>\n')
			webarchive_file.write('\t\t\t\t<Path> {0} </Path>\n'.format(locations['FSYS']['chainsgif'] + struct + '.gif'))
			webarchive_file.write('\t\t\t\t<DownloadPDBFile> {0} </DownloadPDBFile>\n'.format(locations['FSYS']['chains'] + struct + '.pdb'))
			webarchive_file.write('\t\t\t</GeneralStructure_gif>\n')
			webarchive_file.write('\t\t\t<PDB_URL> {0} </PDB_URL>\n'.format(locations['FIXED']['pdbmainsemiurl']+pdbname))
			if opm_data[pdbname]['is_representative']:
				webarchive_file.write('\t\t\t<OPM_URL> {0} </OPM_URL>\n'.format(locations['FIXED']['opmrefsemiurl']+pdbname))
			else:
				webarchive_file.write('\t\t\t<OPM_URL> {0} </OPM_URL>\n'.format(locations['FIXED']['opmrelsemiurl']+pdbname))
			webarchive_file.write('\t\t</General>\n')
			# STRUCTCLASS
			webarchive_file.write('\t\t<StructureInformation>\n')
			webarchive_file.write('\t\t\t<Member> {0} </Member>\n'.format(struct))
			webarchive_file.write('\t\t\t<Class> {0} </Class>\n'.format(neighbors[struct][0]))
			webarchive_file.write('\t\t\t<TMdomains> {0} </TMdomains>\n'.format(neighbors[struct][1]))
			webarchive_file.write('\t\t\t<SequenceNeighbors> {0} </SequenceNeighbors>\n'.format(neighbors[struct][2]))
			webarchive_file.write('\t\t\t<StructureNeighbors> {0} </StructureNeighbors>\n'.format(neighbors[struct][3]))
			webarchive_file.write('\t\t\t<TotalNeighbors> {0} </TotalNeighbors>\n'.format(neighbors[struct][4]))
			webarchive_file.write('\t\t\t<Images>\n')
			webarchive_file.write('\t\t\t\t<RadialDistribution> {0} </RadialDistribution>\n'.format(locations['FSYS']['polar'] + 'p_' + struct + '.png'))
			webarchive_file.write('\t\t\t\t<AlignmentsInDensityPlot> {0} </AlignmentsInDensityPlot>\n'.format(locations['FSYS']['densityscatter'] + 'ds_' + struct + '.png'))
			webarchive_file.write('\t\t\t\t<ResDistanceDistribution> {0} </ResDistanceDistribution>\n'.format(locations['FSYS']['distributions'] + 'distr_' + struct + '.png'))
			webarchive_file.write('\t\t\t\t<Topology> {0} </Topology>\n'.format(locations['FSYS']['topologies'] + 'top_' + struct + '.jpeg'))
			webarchive_file.write('\t\t\t</Images>\n')
			webarchive_file.write('\t\t\t<AdditionalFiles>\n')
			webarchive_file.write('\t\t\t\t<ResDistanceDistribution> {0} </ResDistanceDistribution>\n'.format(locations['FSYS']['distributions'] + 'distr_' + struct + '.dat'))
			webarchive_file.write('\t\t\t\t<StructureAlignments> {0} </StructureAlignments>\n'.format(locations['FSYS'][neighbors[struct][0]] + str(neighbors[struct][1]) + '/' + locations['TREE']['straln'] + 'str_' + struct + '_pdb.dat'))
			webarchive_file.write('\t\t\t</AdditionalFiles>\n')
			webarchive_file.write('\t\t</StructureInformation>\n')
			webarchive_file.write('\t\t<Symmetry>\n')
			webarchive_file.write('\t\t\t<CE-Symm Algorithm="jCE-symm" Version="2.2" >\n')
			if cesymm_text:
				webarchive_file.write(cesymm_text)
			else:
				webarchive_file.write('Not available for this entry\n')
#			webarchive_file.write('\t\t\t</CE-Symm>\n')
			webarchive_file.write('\t\t\t<SymD Algorithm="SymD" Version="1.6" >\n')
			if symd_text:
				webarchive_file.write(symd_text)
			else:
				webarchive_file.write('Not available for this entry\n')
			webarchive_file.write('\t\t</Symmetry>\n')
			webarchive_file.write('\t\t<SymmetryAnalysis>\n')
			if analysis_text:
				webarchive_file.write(analysis_text)
			else:
				webarchive_file.write('Not available for this entry\n')
#			webarchive_file.write('\t\t\t</SymD>\n')
			webarchive_file.write('\t\t</SymmetryAnalysis>\n')
			webarchive_file.write('\t\t<Transfer>\n')
			if transfer_text:
				webarchive_file.write(transfer_text)
			else:
				webarchive_file.write('No results for this entry\n')
			webarchive_file.write('\t\t</Transfer>\n')
			webarchive_file.write('\t</Structure>\n')
	webarchive_file.write('</Database>\n')
	webarchive_file.close()


def draw_topology(options, locations, neighbors, names, opm_data):
	def create_cylinder(init_p, end_p, radius, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += "draw cylinder {{{0} {1} {2}}} {{{3} {4} {5}}} radius {6} resolution 100\n".format(init_p[0], init_p[1], init_p[2], end_p[0], end_p[1], end_p[2], radius)
		return text
	
	def create_slab(init_p, end_p, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += "draw cylinder {{{0} {1} {2}}} {{{3} {4} {5}}} radius {6}\n".format(init_p[0], init_p[1], init_p[2], end_p[0], end_p[1], end_p[2], 1)
		return text
	
	def create_sphere(center_p, radius, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += "draw sphere {{{0} {1} {2}}} radius {3}\n".format(center_p[0], center_p[1], center_p[2], radius)
		return text
	
	def create_Bezier(BezierPoints, nsegments, radius, color='', gappy=False):
		print("BezierPoints", BezierPoints)
		centers = []
		text = ''
		N = len(BezierPoints) - 1
		for t in np.arange(0, 1+1/nsegments, 1/nsegments):
			BezierFormula = np.zeros_like(BezierPoints[0])
			for n in range(N+1):
				BezierFormula += ((1-t)**(N-n))*(t**n)*ncr(N,n)*BezierPoints[n]
			centers.append(BezierFormula)
#			centers.append(((1-t)**5)*BezierPoints[0] + 5*t*(1-t)**4*BezierPoints[1] + 10*(t**2)*(1-t)**3*BezierPoints[2] + 10*(t**3)*(1-t)**2*BezierPoints[3] + 5*(t**4)*(1-t)*BezierPoints[4] + (t**5)*BezierPoints[5])
		if color:
			text = "draw color {0}\n".format(color)
		for nc in range(len(centers)-1):
			text += create_sphere(centers[nc], radius, color='')
			if not gappy:
				text += create_cylinder(centers[nc], centers[nc+1], radius, color='')
		text += create_sphere(centers[len(centers)-1], radius, color='')
		return text

	def create_triangle(p1, p2, p3, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += "draw triangle {{{0} {1} {2}}} {{{3} {4} {5}}} {{{6} {7} {8}}}\n".format(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2])
		return text

	def create_rectangle(p1, p2, p3, p4, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += create_triangle(p1, p2, p3, color='')
		text += create_triangle(p3, p4, p1, color='')
		return text

	def create_convex_polygon(points, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		if len(points) < 3:
			raise NameError("For creating a polygon you have to pass at least 3 points")
		text += create_triangle(points[0], points[1], points[2], color='')
		for i in range(3, len(points)):
			text += create_triangle(points[i-1], points[i], points[0], color='')
		return text

	def create_prism(points, hvector, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		text += create_convex_polygon(points, color='')
		text += create_convex_polygon(points+hvector, color='')
		for i in range(len(points)-1):
			rectpoints = np.array([points[i], points[i+1], points[i+1]+hvector, points[i]+hvector])
			text += create_convex_polygon(rectpoints, color='')
		rectpoints = np.array([points[-1], points[0], points[0]+hvector, points[-1]+hvector])
		text += create_convex_polygon(rectpoints, color='')
		return text

	def create_arrow(init_p, end_p, normal, color=''):
		text = ''
		if color:
			text = "draw color {0}\n".format(color)
		vec = end_p - init_p
		normvec = vec/np.linalg.norm(vec)
		normal = np.array([0, 1, 0])
		widthvec = np.cross(normvec, normal)
		p1 = init_p + widthvec
		p2 = end_p - 4*normvec + widthvec
		p3 = end_p - 4*normvec - widthvec
		p4 = init_p - widthvec
		points = np.array([p1, p2, p3, p4])
		text += create_prism(points-normal/4, normal/2, color='')
		p1 = end_p - 4*normvec + widthvec*2
		p2 = end_p
		p3 = end_p - 4*normvec - widthvec*2
		points = np.array([p1, p2, p3])
		text += create_prism(points-normal/4, normal/2, color='')
		return text

	dssp_path = options['dssp_path']
	fnull = open(os.devnull, 'w')	
	for struct in sorted(list(neighbors.keys())):
		if "4bwz" not in struct:
			continue
#		if (neighbors[struct][0] == 'alpha' and neighbors[struct][1] < 12) or (neighbors[struct][0] == 'beta' and neighbors[struct][1] < 20):
#			continue
		print('struct', struct)
		struct_path = locations['FSYSPATH']['chains'] + struct + '.pdb'
		dssp_out_filename = locations['FSYSPATH']['DSSP'] + struct + '.dssp'
		top_figname = locations['FSYSPATH']['topologies'] + 'top_' + struct + '.jpeg'
		### DANGER
		if os.path.exists(top_figname):
			continue
		###
		pdbname = struct[:4]
		chain = struct[5]
		segments = opm_data[pdbname][chain]['segments']
		print(dssp_path, '-i', struct_path, '-o', dssp_out_filename)
		p = subprocess.Popen([dssp_path, '-i', struct_path, '-o', dssp_out_filename], stdout=fnull, stderr=fnull)
		p.wait()

		# Analyze DSSP file and create a handy dictionary
		# Create a DSSP_segment list with elements such as ('A', 10, 34), where 'A' = helix and 'B' = sheet. DSSP: B, E = beta, H, G, I = alpha.
		if not os.path.exists(dssp_out_filename):
			print("Struct", struct, "jumped")
			continue
		dssp_out_file = open(dssp_out_filename, 'r')
		text = dssp_out_file.read().split('\n')
		dssp_out_file.close()

		dssp_to_edo = {'H':'A', 'G':'A', 'I':'A', 'B':'B', 'E':'B'}
		residues = []
		dssp_segments = []
#		dssp_seg = []
		ss_set = set()
		ss_dict = {}
		start_next = False
		for line in text:
			if not line:
				continue
			fields = line.split()
			if start_next and fields[1] != '!':
				if line[16] in list(dssp_to_edo.keys()):
					ss_type = dssp_to_edo[line[16]]
					ss_set.add(int(line[5:10].strip()))
					ss_dict[int(line[5:10].strip())] = ss_type
#					residues.append((fields[2], int(fields[1]), fields[3], ss_type))
#					if not dssp_seg:
#						dssp_seg.append([ss_type])
#						dssp_seg.append(int(fields[1]))
#					elif dssp_seg[0] == ss_type:
#						dssp_seg.append(int(fields[1]))
#					elif dssp_seg[0] != ss_type:
#						dssp_segments.append(dssp_seg)
#						dssp_seg = [[ss_type], int(fields[1])]
#				else:
#					residues.append([fields[2], int(fields[1]), fields[3], 'L'])
#					if dssp_seg:
#						dssp_segments.append(dssp_seg)
#						dssp_seg = []
			if fields[0].strip() == '#':
				start_next = True

#		print(dssp_segments)

#		if neighbors[struct][0] == 'alpha':
		for seg in opm_data[struct[:4]][struct[5]]['segments']:
			for nres in range(int(seg[0]), int(seg[1])+1):
				ss_set.add(nres)
				if nres not in list(ss_dict.keys()):
					ss_dict[nres] = 'A'
				elif ss_dict[nres] != 'A':
					print("Inconsistency!", nres, ss_dict[nres], 'A')
		new_dssp_segments = []
		dssp_seg = []
		for nres in opm_data[struct[:4]]['FROM_PDB'][struct[5]]['RESIDS']:
			nres = int(nres)
			if nres in ss_set:
				if not dssp_seg:
					dssp_seg.append(ss_dict[nres])
					dssp_seg.append(nres)
				elif dssp_seg[0] == ss_dict[nres]:
					dssp_seg.append(nres)
				elif dssp_seg[0] != ss_dict[nres]:
					new_dssp_segments.append(dssp_seg)
					dssp_seg = []
					dssp_seg.append(ss_dict[nres])
					dssp_seg.append(nres)
			else:
				if dssp_seg:
					new_dssp_segments.append(dssp_seg)
					dssp_seg = []
		dssp_segments = new_dssp_segments

		print(dssp_segments)

		for n_dssp_seg in range(len(dssp_segments)):
			if dssp_segments[n_dssp_seg][0][0] == 'A' and len(dssp_segments[n_dssp_seg]) < 7:   # No helices with less than 6 CA
				dssp_segments[n_dssp_seg][0][0] == 'X'
			if dssp_segments[n_dssp_seg][0][0] == 'B' and len(dssp_segments[n_dssp_seg]) < 4:   # No sheets with less than 3 CA in a row
				dssp_segments[n_dssp_seg][0][0] == 'X'
		

		# Read atom positions
		struct_file = open(struct_path, 'r')
		text = struct_file.read().split('\n')
		struct_file.close()
		bb_at = ['N', 'CA', 'C']
		residue = {}
		residue_list = []
		residue_dict = {}
		for line in text:
			if not line:
				continue
			fields = line.split()
			if not (fields[0] == 'ATOM' and line[13:16].strip() in bb_at):
				continue
			if not residue:
				residue['index'] = int(line[22:26].strip())
				residue_dict[int(line[22:26].strip())] = {}
				residue_dict[int(line[22:26].strip())][line[13:16].strip()] = np.array([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
				prev_index = int(line[22:26].strip())
			elif residue and int(line[22:26].strip()) != prev_index:
				residue_list.append(residue)
				residue_dict[int(line[22:26].strip())] = {}
				residue = {'index':int(line[22:26].strip())}
				prev_index = int(line[22:26].strip())
			residue[line[13:16].strip()] = np.array([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
			residue_dict[int(line[22:26].strip())][line[13:16].strip()] = np.array([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])

		exclude = set()
		for resid in list(residue_dict.keys()):
			if len(list(residue_dict[resid].keys())) != 3:
				exclude.add(resid)
		for resid in exclude:
			del residue_dict[resid]

		# Retrieve limits
		pdbname_path = locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb'
		pdbname_file = open(pdbname_path, 'r')
		text = pdbname_file.read().split('\n')
		pdbname_file.close()
		lim_sup = ''
		lim_inf = ''
		for line in text:
			if not line or ('ATOM' not in line and 'HETATM' not in line):
				continue
			fields = line.split()
			if fields[3] == 'DUM' and fields[2] == 'O' and not lim_sup:
				lim_sup = float(line[46:54].strip())
			elif fields[3] == 'DUM' and fields[2] == 'N' and not lim_inf:
				lim_inf = float(line[46:54].strip())
			if lim_sup and lim_inf:
				break
		print("lim_sup", lim_sup, "lim_inf", lim_inf)

		# for each DSSP_segment which is 'A', take the first 4 CA and calculate the center, and the last 4 and do the same. Then calculate the distance vector, norm it and apply it to the two points so that they get translated by 1 A away from the original vector. Those are the extrema.
		# for each DSSP_segment which is 'B', take the first and last 4 bb atoms (only N, CA, C) and calculate the centers. Then calculate the distance vector, norm it and apply it to the two points so that they get translated by 1 A away from the original vector. Those are the extrema. Also record the direction 
		# if the segment crosses the membrane, also record the first and last center (4 CAs for 'A', 4 BB atoms for 'B') of that element falling inside the membrane. These are the pivots.
		new_dssp_segments = []
		for dssp_seg in dssp_segments:
			if not (dssp_seg[0][0] == 'A' or dssp_seg[0][0] == 'B'):
				continue
			if len(dssp_seg) < 3:   # If there only is one residue...
				continue
			init_point = ''
			end_point = ''
			init_mem = ''
			end_mem = ''
			if dssp_seg[0][0] == 'A':
				tail = 3
			elif dssp_seg[0][0] == 'B':
				tail = 1
			if len(dssp_seg)-tail < 3:
				continue
			for i in range(1, len(dssp_seg)-tail):
				if dssp_seg[0][0] == 'A':
					gcenter = (residue_dict[dssp_seg[i]]['CA'] + residue_dict[dssp_seg[i+1]]['CA'] + residue_dict[dssp_seg[i+2]]['CA'] + residue_dict[dssp_seg[i+3]]['CA'])/4
				elif dssp_seg[0][0] == 'B':
					gcenter = (residue_dict[dssp_seg[i]]['N'] + residue_dict[dssp_seg[i]]['CA'] + residue_dict[dssp_seg[i]]['C'] + residue_dict[dssp_seg[i+1]]['N'] + residue_dict[dssp_seg[i+1]]['CA'] + residue_dict[dssp_seg[i+1]]['C'])/6

#				if gcenter[2] > lim_sup:
#					mempos = 'above'
#				elif gcenter[2] < lim_inf:
#					mempos = 'below'
#				else:
#					mempos = 'inside'

				if init_point == '':
					init_point = np.copy(gcenter)
#					if mempos == 'inside':
#						init_mem = np.copy(gcenter)
#					mempos_old = mempos
					continue

#				if mempos_old != mempos:
#					if init_mem == '':
#						init_mem = np.copy(gcenter)
#					else:
#						end_mem = np.copy(gcenter)
			end_point = np.copy(gcenter)
#			if init_mem != '' and (end_mem == ''):
#				end_mem = np.copy(gcenter)
			# The ending points get translated in order to get to the first and last CA.
			# First, we calculate the distance vector between init and end points. Then we calculcate what's the distance between the init/end point and the plane passing for the first/last CA and whose normal is the normed distance vector. Then we translate the init/end point along the distance vector.
			print("init, end", dssp_seg[1], init_point, dssp_seg[-1], end_point)
			dist_vec = end_point - init_point
			norm_dist_vec = dist_vec/np.linalg.norm(dist_vec)
			if dssp_seg[0][0] == 'A':
				add_dist_init = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[1]]['CA'] - init_point)))
				add_dist_end = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[-1]]['CA'] - end_point)))
			elif dssp_seg[0][0] == 'B':
				add_dist_init = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[1]]['N'] - init_point)))
				add_dist_end = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[-1]]['C'] - end_point)))
			init_point = init_point - add_dist_init*norm_dist_vec
			end_point = end_point + add_dist_end*norm_dist_vec

			if not ((init_point[2] > lim_sup and end_point[2] > lim_sup) or (init_point[2] < lim_inf and end_point[2] < lim_inf)):
				if np.dot(end_point - init_point, np.array([0, 0, 1])) == 0:
					init_mem = init_point
					end_mem = end_point
				else:	
					line_param_sup = (lim_sup - init_point[2])/norm_dist_vec[2]
					point_sup = init_point + norm_dist_vec*line_param_sup
					line_param_inf = (lim_inf - init_point[2])/norm_dist_vec[2]
					point_inf = init_point + norm_dist_vec*line_param_inf
					if np.dot(end_point - init_point, point_inf - point_sup) > 0:
						init_mem = point_sup
						end_mem = point_inf
					else:
						init_mem = point_inf
						end_mem = point_sup
					i_im = init_mem - init_point
					im_em = end_mem - init_mem
					em_e = end_point - end_mem
					if np.dot(i_im, im_em) < 0:
						init_mem = init_point
					if np.dot(em_e, im_em) < 0:
						end_mem = end_point
			else:
				init_mem = ''
				end_mem = ''
			dssp_seg[0] = [dssp_seg[0][0], init_point, end_point, init_mem, end_mem]
			new_dssp_segments.append(dssp_seg)

		dssp_segments[:] = new_dssp_segments[:]

		# Perform rotaions! Mind to only rotate atoms C, CA, and N. Place first pivot at (0, 0, z). Translate all the complex. Then place second pivot at (x, 0, z) rotating all the complex. From the third on, place the pivot at (x, 0, z) rotating only residues after the last residue of the preceding loop / the corresponding ss structure (depending if pivot is initial or final). Stop at the (n-1)-th pivot.
		pivots = []
		for dssp_seg in dssp_segments:
			print(dssp_seg)
			if dssp_seg[0][3] != '' and dssp_seg[0][4] != '':
### CHANGE? dssp_seg[0][1] or dssp_seg[0][3]? UHMMMM
				if not pivots:
					pivots.append([dssp_seg[0][1], sorted(list(residue_dict.keys()))])
					pivots.append([dssp_seg[0][2], [x for x in sorted(list(residue_dict.keys())) if x > dssp_seg[-1]]])
				else:
					pivots.append([dssp_seg[0][1], [x for x in sorted(list(residue_dict.keys())) if x >= dssp_seg[1]]])
					pivots.append([dssp_seg[0][2], [x for x in sorted(list(residue_dict.keys())) if x > dssp_seg[-1]]])
		if not pivots:
			pivots.append([residue_dict[[x for x in sorted(list(residue_dict.keys()))][0]], sorted(list(residue_dict.keys()))])
			pivots.append([residue_dict[[x for x in sorted(list(residue_dict.keys()))][-1]], []])
		print("pivots", pivots)

		for resid in sorted(list(residue_dict.keys())):
			print(resid, list(residue_dict[resid].keys()))
			for atname in bb_at:
				residue_dict[resid][atname] = residue_dict[resid][atname] - pivots[0][0]
		trasl = pivots[0][0]
		for npiv in range(len(pivots)):
			pivots[npiv][0] = pivots[npiv][0] - trasl
		for dssp_seg in dssp_segments:
			for i in range(1,5):
				if dssp_seg[0][i] != '':
					dssp_seg[0][i] = dssp_seg[0][i] - trasl

		for npiv in range(len(pivots)-1):
#		for npiv in range(1):
			# Find the rotation matrix
			vec = pivots[npiv+1][0] - pivots[npiv][0]
			vec[2] = 0.0
			ax = np.array([1, 0, 0])
			vec = vec/(np.linalg.norm(vec))
			print("norm vec:", np.linalg.norm(vec))
			R =  np.zeros((3,3))
			R[0][0] = vec[0]*ax[0]+vec[1]*ax[1]
			R[0][1] = -(vec[0]*ax[1]-ax[0]*vec[1])
			R[1][0] = -R[0][1]
			R[1][1] = R[0][0]
			R[2][2] = 1.0
			print("axis:", ax, "vector:", vec, "matrix:", R)
			#print("Rotation matrix:", R)
			# Rotate all atoms after the ss structure to whom the pivot belongs
			for resid in pivots[npiv][1]:
				for atname in bb_at:
					residue_dict[resid][atname] = np.dot(R, (residue_dict[resid][atname] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0])
			# Rotate all pivots after this
			for nnpiv in range(npiv+1, len(pivots)):
#				print(npiv, nnpiv, pivots[npiv][0], pivots[nnpiv][0])
#				print("pivot before:", pivots[nnpiv][0], np.linalg.norm(pivots[nnpiv][0]))
#				print("pivot trasl:", pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]), np.linalg.norm(pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0])))
#				print("pivot rot:", np.dot(R, (pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]))), np.linalg.norm(np.dot(R, (pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0])))))
				pivots[nnpiv][0] = np.dot(R, (pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0])
#				print("pivot after:", pivots[nnpiv][0], np.linalg.norm(pivots[nnpiv][0]))
			# Rotate all endpoints after this
			for dssp_seg in dssp_segments:
				if dssp_seg[1] >= pivots[npiv][1][0]:
					for i in range(1,5):
						if dssp_seg[0][i] != '':
							dssp_seg[0][i] = np.dot(R, (dssp_seg[0][i] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0])

		# Compressing representation
		if neighbors[struct][0] == 'alpha':
			for npiv in range(2,len(pivots)-1,2):
				vec = pivots[npiv+1][0] - pivots[npiv][0]
				# The compression is just a backtranslation of all pivots after the one considered. The backtranslation is along the x axis and amounts to the double sine of the tilt angle of that segment
				backtrans = np.array([-1,0,0])*2*abs(vec[0])
				print("backtrans", backtrans, np.linalg.norm(backtrans))
				# We have to check that no pivot falls in the semi-plain before the one defined by the previous segment and the y-axis
				vecprev = pivots[npiv-1][0] - pivots[npiv-2][0]
				plainnorm = np.cross(vecprev, np.array([0,1,0]))/np.linalg.norm(np.cross(vecprev, np.array([0,1,0])))
				do_compression = True
				# Check if moved helix does not clash
				checkv_n = np.linalg.norm(pivots[npiv+1][0] + backtrans - pivots[npiv-1][0])
				checkv = (pivots[npiv+1][0] + backtrans - pivots[npiv-1][0])/checkv_n
				if not (np.dot(checkv, plainnorm)*np.dot(plainnorm, np.array([1,0,0])) > 0 and checkv_n*np.dot(checkv, plainnorm) > 5):
					print("No compression for", npiv, "due to", npiv+1, checkv_n, checkv, plainnorm, np.dot(checkv, [1,0,0]), np.dot(checkv, plainnorm), checkv_n*np.dot(checkv, plainnorm))
					do_compression = False
				# Check if other helices do not clash with moved helix
				veccompr = vec + backtrans
				plainnorm2 = np.cross(vecprev, np.array([0,1,0]))/np.linalg.norm(np.cross(vecprev, np.array([0,1,0])))
				for nnpiv in range(npiv+2, min(len(pivots), npiv+4)):
					checkv2_n = np.linalg.norm(pivots[nnpiv][0] + backtrans - pivots[npiv][0])
					checkv2 = (pivots[nnpiv][0] + backtrans - pivots[npiv][0])/checkv2_n
					if not (np.dot(checkv2, plainnorm2)*np.dot(plainnorm2, np.array([1,0,0])) > 0 and checkv_n*np.dot(checkv2, plainnorm2) > 6):
						print("No compression for", npiv, "due to", nnpiv, checkv2_n, checkv2, plainnorm2, np.dot(checkv2, [1,0,0]), np.dot(checkv2, plainnorm2), checkv2_n*np.dot(checkv2, plainnorm2))
						do_compression = False
				if do_compression:
					print("Compression!", npiv)
					for resid in pivots[npiv+1][1]:
						for atname in bb_at:
							residue_dict[resid][atname] = residue_dict[resid][atname] + backtrans
					for nnpiv in range(npiv+1, len(pivots)):
						pivots[nnpiv][0] = pivots[nnpiv][0] + backtrans
					for dssp_seg in dssp_segments:
						if dssp_seg[1] == pivots[npiv][1][0]:
							for i in [2,4]:
								if dssp_seg[0][i] != '':
									dssp_seg[0][i] = dssp_seg[0][i] + backtrans
						elif dssp_seg[1] > pivots[npiv][1][0]:
							for i in range(1,5):
								if dssp_seg[0][i] != '':
									dssp_seg[0][i] = dssp_seg[0][i] + backtrans

		# Further compression of the x axis
		if neighbors[struct][0] == 'alpha':
			cfactor = 2
		elif neighbors[struct][0] == 'beta':
			cfactor = 2
		for resid in list(residue_dict.keys()):
			for atname in bb_at:
				residue_dict[resid][atname][0] = residue_dict[resid][atname][0]/cfactor
		for nnpiv in range(npiv+1, len(pivots)):
			pivots[nnpiv][0][0] = pivots[nnpiv][0][0]/cfactor
		for dssp_seg in dssp_segments:
			for i in range(1,5):
				if dssp_seg[0][i] != '':
					dssp_seg[0][i][0] = dssp_seg[0][i][0]/cfactor

		# For each DSSP_segment, write the command to create a cylinder (A) or a slab (B) with those endpoints.
		text = ''
		text += 'color Display Background white\n'
		text += 'rotate x by -90\n'
		big = 1.5
		small = 0.2
		BezierCA = {}
		for dssp_seg in dssp_segments:
			if dssp_seg[0][0] == 'A':
				text += create_cylinder(dssp_seg[0][1], dssp_seg[0][2], big, 'orange')
			elif dssp_seg[0][0] == 'B':
				normal = np.cross(residue_dict[dssp_seg[2]]['N'] - residue_dict[dssp_seg[1]]['CA'], residue_dict[dssp_seg[2]]['CA'] - residue_dict[dssp_seg[1]]['N'])
				normal = normal/np.linalg.norm(normal)
				text += create_arrow(dssp_seg[0][1], dssp_seg[0][2], normal, 'blue2')
			vec_dist = dssp_seg[0][2] - dssp_seg[0][1]
			norm_vec_dist = vec_dist/(np.linalg.norm(vec_dist))
			init_handle = dssp_seg[0][1] - norm_vec_dist
			end_handle = dssp_seg[0][2] + norm_vec_dist
			BezierCA[dssp_seg[1]] = {-1 : init_handle, 0 : dssp_seg[0][1], 1 : ''}
			BezierCA[dssp_seg[-1]] = {-1 : '', 0 : dssp_seg[0][2], 1 : end_handle}
			
		# for each DSSP_segment loop (not counting first and last) take the endpoint (if any) of the preceding SS piece. This is the first Bezier point. Then, take the normal vector of the preceding SS piece. The point at 1 A from the first point in that direction is the first Bezier handle. Second and third handles are the C and N atoms. Second point is the CA. Fourth handle is the point 1 A away of the CA in the opposite direction of the vector from N to C of that residue. Then, build the Bezier curve. Take parameter at (0, 0.05, 0.1, ..., 1.0) and in each site put a sphere, and between each two sites put a cylinder.
		
		init = True
		inside_a_ss = False
		for resid in sorted(list(residue_dict.keys())):
			if init:
				if resid in list(BezierCA.keys()):
					inside_a_ss = True
					continue
				else:
					BezierCA[resid] = {-1 : '', 0 : residue_dict[resid]['CA'], 1 : end_handle}
				init = False
			if not inside_a_ss:
				if resid in list(BezierCA.keys()):
					if BezierCA[resid][1] == '':
						inside_a_ss = True
						continue
				res_N = residue_dict[resid]['N']
				res_C = residue_dict[resid]['C']
				init_handle = residue_dict[resid]['CA'] - (res_C - res_N)/(np.linalg.norm(res_C - res_N))
				end_handle = residue_dict[resid]['CA'] + (res_C - res_N)/(np.linalg.norm(res_C - res_N))
				BezierCA[resid] = {-1 : init_handle, 0 : residue_dict[resid]['CA'], 1 : end_handle}
			if inside_a_ss and resid in list(BezierCA.keys()) and (BezierCA[resid][-1] == ''):
				inside_a_ss = False

###
#		resids_in_BezierCA = sorted(list(BezierCA.keys()))
#		for nresid in range(len(resids_in_BezierCA)-1):
#			resid = resids_in_BezierCA[nresid]
#			if BezierCA[resid][1] != '' and  BezierCA[resids_in_BezierCA[nresid+1]][-1] != '':
#				BezierPoints = (BezierCA[resid][0], BezierCA[resid][1], residue_dict[resid]['C'], residue_dict[resids_in_BezierCA[nresid+1]]['N'], BezierCA[resids_in_BezierCA[nresid+1]][-1], BezierCA[resids_in_BezierCA[nresid+1]][0])
#				print('BezierPoints', BezierPoints)
#				text += create_Bezier(BezierPoints, 100, small, 'blue')
#				text += create_Bezier(BezierPoints, 100, small, 'blue')
####

		resids_in_BezierCA = sorted(list(BezierCA.keys()))
		init = False
		end = False
		gap = False
		count = 0
		BezierPoints = []
		for nresid in range(len(resids_in_BezierCA)):
			resid = resids_in_BezierCA[nresid]
			if BezierCA[resid][-1] == '':
				BezierPoints = []
				init = True
			if BezierCA[resid][1] == '' or count == 3 or (nresid+1 < len(resids_in_BezierCA) and BezierCA[resid][-1] != '' and np.linalg.norm(BezierCA[resids_in_BezierCA[nresid+1]][0] - BezierCA[resid][0]) > 4.5):
				end = True
			if BezierCA[resid][-1] != '' and nresid-1 > 0 and np.linalg.norm(BezierCA[resid][0] - BezierCA[resids_in_BezierCA[nresid-1]][0]) > 4.5:
				end = True
				gap = True
			if end:
				BezierPoints.append(BezierCA[resid][-1])
			BezierPoints.append(BezierCA[resid][0])
			count += 1
			if init:
				BezierPoints.append(BezierCA[resid][1])
				init = False
			if end:
				print("BezierPoints", BezierPoints)
				if gap:
					text += create_Bezier(BezierPoints, 10, small, 'gray', gappy=True)
				else:
					text += create_Bezier(BezierPoints, 100, small, 'gray')
				init = False
				end = False
				gap = False
				count = 1
				BezierPoints = []
				BezierPoints.append(BezierCA[resid][0])
				BezierPoints.append(BezierCA[resid][1])


		limright = 0
		limleft = 0
		for resid in sorted(list(residue_dict.keys())):
			if limright < residue_dict[resid]['CA'][0]:
				limright = residue_dict[resid]['CA'][0]
			if limleft > residue_dict[resid]['CA'][0]:
				limleft = residue_dict[resid]['CA'][0]

		limup = lim_sup - trasl[2]
		limdown = lim_inf - trasl[2]
#		for dssp_seg in dssp_segments:
#			for i in range(3,5):
#				if dssp_seg[0][i] != '':
##					if limup < dssp_seg[0][i][2]:
#						limup = dssp_seg[0][i][2]
#					if limdown > dssp_seg[0][i][2]:
#						limdown = dssp_seg[0][i][2]

#		text = 'display resize {0} {1}\n'.format(10000, int(10000*(abs(limup - limdown)/abs(limright - limleft)))) + text
		if int(opm_data[struct[:4]][struct[5]]['ntm']) < 13:
			text = 'display resize {0} {1}\n'.format(8000, 6000) + text
		else:
			text = 'display resize {0} {1}\n'.format(16000, 6000) + text
		text += 'display projection Orthographic\n'
		text += 'display depthcue off\n'

		p1 = np.array([limright, 50, limup])
		p2 = np.array([limright, 50, limdown])
		p3 = np.array([limleft, 50, limdown])
		p4 = np.array([limleft, 50, limup])
#		text += "material change opacity 0.3\n"
		text += create_rectangle(p1, p2, p3, p4, color='silver')

		scale = 0.01
		text += 'scale to {0}\n'.format(scale)
		text += 'translate to {0} {1} {2}\n'.format(-(limright - limleft)/2*scale, 0, -(limup + limdown)/2*scale)
		text += 'axes location off\n'

		tgafig = locations['FSYSPATH']['topologies'] + 'tmp_' + struct + '.tga'
		text += 'render TachyonInternal "{0}"\n'.format(tgafig)
		text += 'quit\n'

		vmdscript_filename = locations['FSYSPATH']['topologies'] + struct + '.tcl'
		vmdscript_file = open(vmdscript_filename, 'w')
		vmdscript_file.write(text)
		vmdscript_file.close()

		p = subprocess.Popen([options['vmd_path'], "-dispdev",  "text",  "-e", vmdscript_filename], stdout=fnull)
		p.wait()

#		tmpfig = locations['FSYSPATH']['topologies'] + 'tmp_' + struct + '.jpeg'		
#		p = subprocess.Popen(["convert", tgafig, tmpfig])
#		p.wait()
		
		im = Image.open("{0}".format(tgafig))
		pix = np.asarray(im)

		pix = pix[:,:,0:3] # Drop the alpha channel
		idx = np.where(pix-255)[0:2] # Drop the color when finding edges
		box = list(map(min,idx))[::-1] + list(map(max,idx))[::-1]

		region = im.crop(box)
		region_pix = np.asarray(region)

		fig, ax = plt.subplots()
		img = ax.imshow(region_pix)
		ax.axis('off')
		fig.savefig(top_figname, bbox_inches='tight', dpi=1000)
		os.remove(tgafig)


def webcontent_do_all(options, filters, locations, opm_data, table):
	print(options, filters, locations)
	if not (options and filters and locations):
		options, filters, locations = genfsys_opm.filesystem_info()
	neighbors, names, opm_data, table, equivalences = neighborlists_creator(options, locations, opm_data, table)
#	return
#	webcontent(locations, neighbors, names, opm_data)
#	return
#	polar(options, locations, table)
#	gif_instr_creator(options, locations, neighbors, names)
#	transfer_estimator(locations, neighbors, opm_data, equivalences)
#	scatterplots(options, locations, neighbors, table)
	webcontent(locations, neighbors, names, opm_data)
#	draw_topology(options, locations, neighbors, names, opm_data)

webcontent_do_all({}, {}, {}, {}, {})
#options, filters, locations = genfsys_opm.filesystem_info()
#neighbors, names, opm_data, table = neighborlists_creator(options, locations, {}, {})
#webcontent(locations, neighbors, names, opm_data)
