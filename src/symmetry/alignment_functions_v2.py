# Author: Antoniya A. Aleksandrova; script fragments from Jose Rodrigues
# Date: 19 Oct 2016
# Language: Python 3.5 (biopython)
from __future__ import print_function, division

import argparse
import os
import numpy as np

from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, Selection
from Bio.PDB.Polypeptide import is_aa

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.Seq import Seq
from Bio import SeqIO


def get_pdb_sequence(structure):
    """
    Retrieves the AA sequence from a PDB structure.
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
    return seq

def get_pdb_sequence_with_chains(structure):
    """
    Retrieves the AA sequence from a PDB structure. It's a list that looks like [(5, 'R', 'A'), (6, 'E', 'A'), (7, 'H', 'A'), (8, 'W', 'A'),...]
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'),r.get_parent().get_id())
    seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
    return seq

def calculate_identity(sequenceA, sequenceB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned.
    """

    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
    matches = [sa[i] == sb[i] for i in range(sl)]
    seq_id = (100 * sum(matches)) / sl

    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    gap_id = (100 * sum(matches)) / gapless_sl
    return (seq_id, gap_id)

def align_sequences(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = get_pdb_sequence(structA)
    resseq_B = get_pdb_sequence(structB)
#    print(resseq_B)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
#    print(aligned_A)
#    print(aligned_B)

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    # Gapless alignment
    # def _trimmer(sequence):
    #     """Returns indices of first and last ungapped position"""

    #     leading = [i for (i, aa) in enumerate(sequence) if aa != '-'][0]
    #     trailing = [i for (i, aa) in enumerate(sequence[::-1]) if aa != '-'][0]

    #     trailing = len(sequence) - trailing
    #     return (leading, trailing)

    # lead_A, trail_A = _trimmer(aligned_A)
    # lead_B, trail_B = _trimmer(aligned_B)

    # lead = max(lead_A, lead_B)
    # trail = min(trail_A, trail_B)
    # trim_aln_A = aligned_A[lead:trail]
    # trim_aln_B = aligned_B[lead:trail]
    # mismatch = ''.join(['+' if a!=b else ' ' for (a,b) in zip(trim_aln_A, trim_aln_B)])

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)
    # return ((trim_aln_A, trim_aln_B, mismatch), seq_id, g_seq_id, mapping)

def align_sequences_with_chains(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = get_pdb_sequence_with_chains(structA)
    resseq_B = get_pdb_sequence_with_chains(structB)
#    print(resseq_B)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
#    print(aligned_A)
#    print(aligned_B)

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[(resseq_A[aa_i_A][0],resseq_A[aa_i_A][2])] = (resseq_B[aa_i_B][0],resseq_B[aa_i_B][2])
            aa_i_A += 1
            aa_i_B += 1

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)

def align_sequences_with_chains_and_gaps(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = get_pdb_sequence_with_chains(structA)
    resseq_B = get_pdb_sequence_with_chains(structB)
#    print(resseq_B)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
    mapping = [[],[]]
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                mapping[0].append(('-','-'))
                mapping[1].append((resseq_B[aa_i_B][0],resseq_B[aa_i_B][2]))
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                mapping[0].append((resseq_A[aa_i_A][0],resseq_A[aa_i_A][2])) 
                mapping[1].append(('-','-'))
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[0].append((resseq_A[aa_i_A][0],resseq_A[aa_i_A][2]))
            mapping[1].append((resseq_B[aa_i_B][0],resseq_B[aa_i_B][2]))
            aa_i_A += 1
            aa_i_B += 1

    return ((aligned_A, aligned_B), mapping)

def parse_structure(spath):
    """Parses a PDB/cif structure"""

    if not os.path.isfile(spath):
        return IOError('File not found: {0}'.format(spath))

    if spath.endswith(('pdb', 'ent')):
        parser = PDBParser(QUIET=True)
    elif spath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(spath))

    sname = os.path.basename(spath.split('.')[0])
    return parser.get_structure(sname, spath)

def opm_sql(opm_dir, inputf):
###### Working with OPM SQL file
    sql=open(opm_dir+"OPM06-05-16.sql","r")
    flag=0
    chain=''
    name=''
    family_name=''
    opm_id=inputf
    family='not applicable'
    #  Check whether the pdb id is related to another structure or is representative
    for line in sql:
        if "-- Dumping data for table `protein`" in line:
            flag=2
        if ("-- Dumping data for table `relatedproteins" in line) and flag==2:
            flag=6
        if flag==2 and inputf in line:
            flds=line.split(',')
            family=flds[1].strip()
            flds=line.split("'")
            name='"'+flds[3]+'"'
            flag=0
        if flag==6 and inputf in line:
            flds=line.split("'")
            if inputf in flds[3]:
                opm_id=flds[1]
                break  
    sql.seek(0)
    flag=0
    #  Get subunits (chains) for the protein
    for line in sql:
        if "-- Dumping data for table `protein`" in line:
            flag=2
        if flag==2 and opm_id in line:
            flds=line.split(',')
            family=flds[1].strip()
            flds=line.split("'")
            name='"'+flds[3]+'"'
            flag=0
        if "-- Dumping data for table `subunits`" in line:
            flag=1
        if flag==1 and opm_id in line:
            flds=line.split()[3]
            chain=chain+flds[1:][:-2]+";"
        if flag==1 and '-- -' in line:
            break
    sql.seek(0)
    #  Get name and family of the protein
    flag=0
    superfml='not_applicable'
    for line in sql:
        if "-- Dumping data for table `family`" in line:
            flag=1
        if flag==1 and line.startswith("("+family+','):
            index=line.split()[1]
            superfml=index[:-5]+"'"
            flag=0
        if "-- Dumping data for table `superfamily`" in line:
            flag=2
        if flag==2 and superfml in line:
            fld=line.split("'")[3]
            family_name='"'+fld+'"'
            flag=0
    sql.close()
    return chain, name, family_name

def cesymm_related_repeats(inputf,cesymm_dir):
    axes=open(cesymm_dir+inputf+".axes","r")
    flag=0
    repeats_all=[]
    repeats_type=[]
    repeats_level=[]
    axes_per_level=[]
    order_of_level=[]
    for line in axes:
        if inputf[0:4] in line:
            flag=1
            flds=line.split()
            order=int(flds[1])
#            repeats_all.append([])
            repeats=flds[8].strip(')/(/\n')
            if ')(' in repeats:
                repeats=repeats.split(')(')
                repeats=[k.split(';') for k in repeats]
            else:
                repeats=[repeats.split(';')]
            if len(repeats_all)<order:
                repeats_all.append(repeats)
                repeats_type.append(flds[2])
                repeats_level.append(order)
                axes_per_level.append(1)
                order_of_level.append(flds[3])
            else:
                for k in repeats:
                    repeats_all[order-1].append(k)
                axes_per_level[order-1]+=1
                if flds[2]!=repeats_type[order-1]:
                    print(("Warning: Combining %s with %s symmetries!")%(repeats_type[order-1],flds[2]))
    if flag==1:
        return np.array(repeats_all), repeats_type, repeats_level, axes_per_level, order_of_level
    else:
        raise SystemExit("Error: %s has no detected symmetry\n" % inputf)

def find_tm_chains(inputf,opm_dir,flag):
# Turn flag to 0 to use Edo's dictionary str_info; Turn to 1 when the structure is not part of Edo's dictionary
    if flag==0:
        str_info=read_data(opm_dir+"opm_archive.txt")
        tm_ch=str_info[inputf]['tmchains']
    else:
        chain=""
    #    if downloaded==1:
    #      chain=chain_id(oriented)
        ppm_results=open(opm_dir+"all_opm_ppm_results_corrected.dat","r")
        for line in ppm_results:
            if line.startswith(inputf):
                flds=line.split()
                chain=flds[1]
        ppm_results.close()
        if chain=="":
            tm_ch=[]
            raise SystemExit(inputf+" has no known TM chains.")
        else:
            tm_ch=chain.split(';')[:-1]
    return(tm_ch)


def strip_tm_chains(dir,inputf,oriented_opm,chains):
	f=open(oriented_opm,'r')
	o=open(dir+"/"+inputf+"_tmp.pdb","w")
	LINELEM="{:76s}{:>2s}\n"
	for line in f:
	  if (line.startswith("ATOM") or line.startswith("TER ") or line.startswith("TER\n")) and line[21:22] in chains:
	    if line[76:78].strip()!='':
	      o.write(line)
	    else:
	      if line.startswith("TER ") or line.startswith("TER\n"):
	        o.write(line)
	      else:
	        atom=line[12:16].strip()
	        elem=atom[0]
	        o.write(LINELEM.format(line[0:76],elem))
	  if line.startswith("HETATM") and line[17:20]=='MSE' and line[21:22] in chains:
	    if line[76:78].strip()!='':
	       o.write("ATOM  "+line[6:])
	    else:
	      atom=line[12:16].strip()
	      elem=atom[0]
	      o.write("ATOM  "+line[6:76]+" "+elem+"\n")
	   
	o.write("END\n")
	f.close()
	o.close()
	
def rotation_axis_big_angle(rot,transl,c, theta):
    """ 
    This is only for cases where the angle is at least 5 degrees different from 0
    rot is the 3x3 rotation matrix
    transl is the translation vector
    c is the cosine of the rotation angle (cos(theta))
    
    We are following the method descirbed in Kim et al. 2010 (SI) and implemented in 
    CE-Symm's function calculateRotationalAxis
     
    """
    sum=0
    rot_ax=[]
    
    for i in range(0,3):
        rot_ax.append(np.sqrt(rot[i,i]-c)) # v_i=np.sqrt(R[i,i]-c)
        sum+=rot_ax[i]*rot_ax[i]
    for i in range(0,3):
        rot_ax[i]=rot_ax[i]/np.sqrt(sum)
        
    # Determine the sign
    d0 = rot[2,1]-rot[1,2] #=2u[0]*sin(theta]
    d1 = rot[0,2]-rot[2,0] #=2u[1]*sin(theta]
    d2 = rot[1,0]-rot[0,1] #=2u[2]*sin(theta]

    s12 = rot[2,1]+rot[1,2] #=2*u[1]*u[2]*(1-cos(theta]]
    s02 = rot[0,2]+rot[2,0] #=2*u[0]*u[2]*(1-cos(theta]]
    s01 = rot[1,0]+rot[0,1] #=2*u[0]*u[1]*(1-cos(theta]]
    print(d0, d1, d2)
    print(s01, s02)
    # Take the biggest d for the sign to ensure numerical stability
    if np.abs(d0)<np.abs(d1): # not d0
        if np.abs(d1) < np.abs(d2): # d2
            if d2>=0: # u[2] positive
                if s02 < 0: 
                    rot_ax[0] = -rot_ax[0];
                if s12 < 0: 
                    rot_ax[1] = -rot_ax[1];
            else: #u[2] negative
                rot_ax[2] = -rot_ax[2];
                if s02 >= 0:
                    rot_ax[0] = -rot_ax[0];
                if  s12 >= 0:
                    rot_ax[1] = -rot_ax[1];
            
        else: #d1
            if(d1>=0):#u[1] positive
                if s01 < 0:
                    rot_ax[0] = -rot_ax[0];
                if s12 < 0:
                    rot_ax[2] = -rot_ax[2];
            else: #u[1] negative
                rot_ax[1] = -rot_ax[1];
                if s01 >= 0:
                    rot_ax[0] = -rot_ax[0];
                if s12 >= 0:
                    rot_ax[2] = -rot_ax[2];
    else: # not d1
        if( np.abs(d0) < np.abs(d2) ): #d2
            if(d2>=0): #u[2] positive
                if s02 < 0:
                    rot_ax[0] = -rot_ax[0];
                if s12 < 0: 
                    rot_ax[1] = -rot_ax[1];
            else: #u[2] negative
                rot_ax[2] = -rot_ax[2];
                if s02 >= 0:
                    rot_ax[0] = -rot_ax[0];
                if s12 >= 0:
                    rot_ax[1] = -rot_ax[1];
        else: #d0
            if(d0>=0): #u[0] positive
                if s01 < 0:
                    rot_ax[1] = -rot_ax[1];
                if s02 < 0: 
                    rot_ax[2] = -rot_ax[2];
            else: #u[0] negative
                rot_ax[0] = -rot_ax[0];
                if s01 >= 0: 
                    rot_ax[1] = -rot_ax[1];
                if s02 >= 0: 
                    rot_ax[2] = -rot_ax[2];
    
    scale=np.dot(transl,np.array(rot_ax))
    screw_transl_vec=scale*np.array(rot_ax)
    perp_transl_vec=transl-screw_transl_vec
    h=np.cross(perp_transl_vec,np.array(rot_ax))/(2.0*np.tan(theta/2.0))
    axis_pos=h+0.5*perp_transl_vec
    
    return rot_ax, np.abs(scale), axis_pos


def rotation_axis_small_angle(rot, transl):
    rot_ax=transl/np.linalg.norm(transl)
    scale=np.dot(transl,np.array(rot_ax))
    axis_pos='na'
    return rot_ax, np.abs(scale), axis_pos

def get_rotation_axis(rotation,transl):
    c=(np.trace(rotation)-1)/2.0 # =cos(theta)
    print(c)
    # c is sometimes slightly out of the [-1,1] range due to numerical instabilities
    if( -1-1e-8 < c and c < -1 ):
        c = -1
    if( 1+1e-8 > c and c > 1 ): 
        c = 1
    if( -1 > c or c > 1 ):
        raise SystemExit("Not a valid rotation matrix")

    theta=np.arccos(c)
    min_angle=5*np.pi/180
    if theta < min_angle:
        rot_ax, screw_transl, axis_pos=rotation_axis_small_angle(rotation, transl)
    else:
        rot_ax, screw_transl, axis_pos=rotation_axis_big_angle(rotation, transl, c, theta)
    return theta, rot_ax, screw_transl, axis_pos


def get_rmsd_cesymm(coord_repeats):
    ### CE-Symm RMSD
    sumSqDist=0
    comparisons=0
    
    for r1 in range(0,len(coord_repeats)):
        for c in range(0,len(coord_repeats[r1])):
            refAtom=coord_repeats[r1][c]
            if refAtom==None:
                continue
            
            nonNullSqDist=0
            nonNullLength=0
            for r2 in range(r1+1,len(coord_repeats)):
                atom=coord_repeats[r2][c]
                if atom!=None:
                    nonNullSqDist +=(np.linalg.norm(refAtom-atom))**2
                    nonNullLength+=1
                
            if nonNullLength>0:
                comparisons+=1
                sumSqDist += 1.0*nonNullSqDist / nonNullLength
    return np.sqrt(sumSqDist/comparisons)
