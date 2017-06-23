# Author: Antoniya A. Aleksandrova; script patches from Jose Rodrigues
# Date: 12 Sep 2016
# Last Update: 3 Apr 2017
# Language: Python 3.5 (biopython)
# Takes a template structure and a target structure and transfers the symmetry from
# the template to the target/mobile structure using CE-Symm alignment; calculates RMSD
# and TM-score between the repeats. 
# To transfer symmetry from a run that contained just a segment of the template, execute:
# > python transfer_symmetry.py 3v8g.C_252-416.pdb 3v8g_A.pdb 


"""
Sequence-based structural alignment of two proteins.
"""

from __future__ import print_function, division

#import argparse
import os
import numpy as np
import sys
from Bio.PDB.Vector import *

from alignment_functions_v2 import *
from itertools import combinations

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)
    
def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def tm_score_calc(mobi_ca,transformed_ca,lmin):
    """Takes the lists of vector coordinates of the initial and transformed structures and calculates a TM score"""
    d0=1.24*(lmin-15)**(1/3.0)-1.8
    tm_score_elem=0

    for i in range(0,len(mobi_ca)):
        tm_score_elem+=1.0/(1+pow(np.linalg.norm(mobi_ca[i]-transformed_ca[i])/d0,2))
#        print(np.linalg.norm(mobi_ca[i]-transformed_ca[i]))
    return tm_score_elem/lmin
    
def rmsd_score(original_ca, transformed_ca):
    """Calculates the contribution to the RMSD from a given set of coordinates"""
    rmsd_elem=0
    for i in range(0,len(transformed_ca)):
        rmsd_elem=rmsd_elem + pow(np.linalg.norm(transformed_ca[i]-original_ca[i]),2)
    rmsd_elem=rmsd_elem/len(transformed_ca)
    return np.sqrt(rmsd_elem)

###### The workhorse ####
def transfer_sym(arg1,arg2):
    inputf = arg1[0:4]
    inputf2=arg2[0:4]
    if '_' in arg1:
        chain=arg1[5:6]
    else:
        chain=""
    if '_' in arg2:
        chain2=arg2[5:6]
    else:
        chain2="" 

    dir=os.getcwd()
    # Find the structure in the OPM database
    if os.path.isfile(dir+"/"+inputf+".pdb"):
      oriented_opm=dir+"/"+inputf+".pdb"
    else:
      raise SystemExit("The %s is not in the current folder."%(inputf+".pdb"))

    cesymm_dir=dir+"/"

    # Remove strange residues and extract only transmembrane chains
    flag=1 # Turn flag to 0 to use Edo's dictionary str_info; Turn to 1 when the structure is not part of Edo's dictionary
    tm_ch=find_tm_chains(inputf, dir+'/', flag)
    strip_tm_chains(dir,inputf,oriented_opm,tm_ch)
    oriented_opm=dir+"/"+inputf+"_tmp.pdb"

    if os.path.isfile(dir+"/"+inputf2+".pdb"):
      oriented_opm_2=dir+"/"+inputf2+".pdb"
    else:
      raise SystemExit("The %s is not in the current folder."%(inputf2+".pdb"))
    tm_ch=find_tm_chains(inputf2, dir+'/', flag)
    strip_tm_chains(dir,inputf2,oriented_opm_2,tm_ch)
    oriented_opm_2=dir+"/"+inputf2+"_tmp.pdb"

    # Load the template
    s_template=parse_structure(oriented_opm)
    if chain=="":
        template=s_template[0]
    else:
        try:
            template=s_template[0][chain]
        except KeyError:
            raise Exception('Chain {0} not found in template structure'.format(chain))

    # Load the target and its movable copy mobile
    s_target=parse_structure(oriented_opm_2)
    s_mobile=parse_structure(oriented_opm_2)
    if chain2=="":
        target=s_target[0]
        mobile=s_mobile[0]
    else:
        try:
            target=s_target[0][chain2]
            mobile=s_mobile[0][chain2]
        except KeyError:
            raise Exception('Chain {0} not found in target structure'.format(chain2))


    # Align sequences to get mapping between residues
    aln, seq_id, gapless_id, res_map = align_sequences_with_chains(template, target)
    print(res_map)
    print("\n")
    # Create a list of the residue ids of the target structure
    target_res_list=[]
    for res in target.get_residues():
        if is_aa(res.get_resname(),standard=True):
            target_res_list.append((res.get_id()[1],res.get_parent().get_id()))
    #print(target_res_list)
    target_sym_descriptors=[]
    transforms=[]
    #template_ca_list, target_ca_list = [], []
    #for temp_res in res_map:
    #    template_ca_list.append(s_template[0][temp_res[1]][temp_res[0]]['CA'])
    #    target_ca_list.append(s_target[0][res_map[temp_res][1]][temp_res[0]]['CA'])

    resseq_A = get_pdb_sequence_with_chains(template) # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm=SeqIO.to_dict(SeqIO.parse(cesymm_dir+arg1[:-4]+".fasta", "fasta"))
    repeats_all,repeats_type, repeats_levels, axes_per_level, order_of_level=cesymm_related_repeats(arg1[:-4],cesymm_dir)
    print(repeats_all,repeats_type)
    #repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    unique_repeats=[]
    mode=""
    for ind, symm in enumerate(repeats_all): 
        positions=sum([sum(1 for c in cesm[symm[k][1]].seq if c!='/') for k in range(0,len(symm))])
        algn_map=[[] for i in range(positions)]
        map_ind=0
        print(symm)
        for sub_symm in symm:               #sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            for k in sub_symm:              #k='1PV6.A_6-94'
                if k not in unique_repeats:
                    unique_repeats.append(k)     
                r=cesm[k]
                r_id=k.split('_')
                ch=r_id[0][-1]
                repeat=r_id[1].split('-')
            #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
                repeat_start=[i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch==sublist[2])] #finds the index of the residue that starts the repeat
            #    print(repeat_start)
                count=0
                for i,c in enumerate(r.seq): #Last edited
                    if c!='/':
                      if c!='-' and c.isupper() and resseq_A[repeat_start[0]+count][1] == c:
            #    		    print(c,resseq_A[repeat_start[0]+count][1],resseq_A[repeat_start[0]+count][0],"\n") 
                        algn_map[i+map_ind].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
                        count+=1
             #           if resseq_A[repeat_start[0]+count][2] not in chains:
             #               chains.append(resseq_A[repeat_start[0]+count][2])
                      if c=='-':
                        algn_map[i+map_ind].append(('-','-'))
                      if c.islower():
                        algn_map[i+map_ind].append(('-','-'))
                        count+=1  
            map_ind+=sum(1 for c in cesm[k].seq if c!='/')
        print(len(algn_map))
        print("The step is ",ind,"\n")
        print(algn_map,"\n")
    #print(chains)


        rotations=2
        tm_score=0
        rmsd=0
        for resi in algn_map:
            if all(elem==('-','-') for elem in resi)==False:
                rotations=len(resi)
                break
        print("rotations:", rotations)



        # Create a list of atoms to be paired during superposition for obtaining the axis of rotation   
        fixed_ca_list, mobile_ca_list = [], []
        fixed_rl=[]
        mobile_rl=[]
        mobile_coord=[]
        fixed_coord=[]

        if repeats_type[ind]=='CLOSED':
            for k in range(0,rotations):
                for i,resi in enumerate(algn_map):
                    index=(rotations-1+k)%rotations
                    if resi[k]!=('-','-') and resi[index]!=('-','-') and resi[k] in res_map and res_map[resi[k]] in target_res_list and resi[index] in res_map and res_map[resi[index]] in target_res_list:
                        fixed_ca_list.append(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'])
                        fixed_coord.append(np.array(list(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].get_vector())))
                        mobile_ca_list.append(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]]['CA'])
                        mobile_coord.append(np.array(list(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]]['CA'].get_vector())))
                        fixed_rl.append(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]].get_id()[1])
                        mobile_rl.append(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]].get_id()[1])
        elif repeats_type[ind]=='OPEN':
            for k in range(1,rotations):
                for i,resi in enumerate(algn_map):
                    index=k-1
                    if resi[k]!=('-','-') and resi[index]!=('-','-') and resi[k] in res_map and res_map[resi[k]] in target_res_list and resi[index] in res_map and res_map[resi[index]] in target_res_list:
                        fixed_ca_list.append(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'])
                        fixed_coord.append(np.array(list(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].get_vector())))
                        mobile_ca_list.append(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]]['CA'])
                        mobile_coord.append(np.array(list(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]]['CA'].get_vector())))
                        fixed_rl.append(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]].get_id()[1])
                        mobile_rl.append(s_mobile[0][res_map[resi[index]][1]][res_map[resi[index]][0]].get_id()[1])

        print("Fixed is: ",fixed_rl)
        print("Mobile is: ",mobile_rl)
        vec=mobile_ca_list[0].get_vector()
        print("Before: ",vec, "and ",fixed_ca_list[0].get_vector())
        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(fixed_ca_list, mobile_ca_list)
        print(si.rotran)
#        cs_mat=np.array([[-0.2640556515444221,0.9644870420482027,-0.006273643961747427],[-0.9383960458009067,-0.2553983792857821,0.23277570552662508],[0.22290687318401162,0.06735290327301197,0.9725102119299432]])

        T=np.array(np.transpose(si.rotran[0]))
#        T=np.array(np.transpose(cs_mat))

        transl=si.rotran[1]
#        ang,rotax=m2rotaxis(si.rotran[0])
        ang, rotax, screw_transl, axis_pos=get_rotation_axis(si.rotran[0],transl)
        print(axis_pos)
        if axis_pos=='na':
            center=sum(fixed_coord)/len(fixed_coord)
            print(center)
            proj=np.dot(center,rotax)
            center_on_axis=np.array([proj*i for i in rotax])
            axis_pos=center-center_on_axis

        
        print("Rotation angle: ",ang*180.0/np.pi,"Axis: ",rotax) 
        print("Transformation RMS is", si.rms)
#        transl_magnitude=np.linalg.norm(np.dot(transl,list(rotax))*list(rotax))
#        scale=np.dot(transl,list(rotax))
#        transl_magnitude=np.linalg.norm(scale*np.array(list(rotax)))
        axis_angle=round(angle_between(list(rotax),(0,0,1))*180/np.pi,2)
        if np.abs(np.dot(list(rotax), (0,0,1)))>0.5:
            axis_type="Parallel;"
        else:
            axis_type="Antiparallel;"
        angle=ang*180.0/np.pi
        transforms.append(T)
        centroid_mobile=sum(mobile_coord)/len(mobile_coord)

        # Apply transformation to coordinates 
        aligned_coord_res=[]
        coord=[]
        for i,resi in enumerate(algn_map):
            if np.any([resi[k]!=('-','-') for k in range(0,len(resi))]):
                aligned_coord_res.append([])
                for k in range(0,len(resi)):
                    if resi[k]!=('-','-') and resi[k] in res_map and res_map[resi[k]] in target_res_list:
                        aligned_coord_res[-1].append(res_map[resi[k]])
                        pt=list(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].get_vector())
                        #pt=pt-centroid_mobile
                        for j in range(0,k):
                            #pt=np.dot((pt-transl),T)
                            #pt=pt-centroid_mobile
                            pt=np.dot((pt-transl),T)
                            centroid_mobile=np.dot(centroid_mobile,T)
                            #pt=pt+centroid_mobile
                        #pt=pt+centroid_mobile
                        coord.append(np.array(pt))
                        s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].set_coord(pt)
                    else:
                        aligned_coord_res[-1].append(('-','-'))
        print(aligned_coord_res[0])
        print("After: ",mobile_ca_list[0].get_vector(), "and ",fixed_ca_list[0].get_vector())
    #    print("Mobile: ",s_mobile[0]['A'][6]['CA'].get_vector(),s_mobile[0]['A'][101]['CA'].get_vector(),s_mobile[0]['A'][219]['CA'].get_vector(),s_mobile[0]['A'][310]['CA'].get_vector())
        centroid=sum(coord)/len(coord)
        print(centroid)
        ax2=axis_pos+10*np.array(rotax)
        target_sym_descriptors.append([float('%.2f' % angle), float('%.2f'% screw_transl), transl, ax2, axis_pos, axis_angle])


    ### Output PDB
#    io = PDBIO()
#    io.set_structure(s_mobile)
#    io.save('super_cesymm.pdb')

    #### Calculate the RMSD & TM-score ####
    unique_repeats=list(unique_repeats)
    print(unique_repeats)
    positions=sum(1 for c in cesm[unique_repeats[0]].seq if c!='/')    #ind=1, positions=97
    algn_map=[[] for i in range(positions)]
    coord_repeats=[[] for i in range(positions)]

    # Load the alignments of all repeats of the template against each other
    for k in unique_repeats:              
        r=cesm[k]
        r_id=k.split('_')
        ch=r_id[0][-1]
        repeat=r_id[1].split('-')
    #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
        repeat_start=[i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch==sublist[2])] #finds the index of the residue that starts the repeat
    #    print(repeat_start)
        count=0
        for i,c in enumerate(r.seq): #Last edited
            if c!='/':
              if c!='-' and c.isupper() and resseq_A[repeat_start[0]+count][1] == c:
    #    		    print(c,resseq_A[repeat_start[0]+count][1],resseq_A[repeat_start[0]+count][0],"\n") 
                algn_map[i].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
                coord_repeats[i].append(s_mobile[0][resseq_A[repeat_start[0]+count][2]][resseq_A[repeat_start[0]+count][0]]['CA'].get_vector().get_array())
                count+=1
     #           if resseq_A[repeat_start[0]+count][2] not in chains:
     #               chains.append(resseq_A[repeat_start[0]+count][2])
              if c=='-':
                algn_map[i].append(('-','-'))
                coord_repeats[i].append(None)
              if c.islower():
                algn_map[i].append(('-','-'))
                coord_repeats[i].append(None)
                count+=1 

    num_repeats=len(unique_repeats)
    # Find the length of the target protein
    prot_len=0
    resseq_B=get_pdb_sequence_with_chains(target)
    prot_len=len(resseq_B)
    print("Length is: ", prot_len)
    
    rmsd=0
    tm_score=0
    n=0
    dummy=list(range(0,len(unique_repeats))) # a list of the N numbers representing the N unique repeats
    combo_list=list(combinations(dummy,2)) # a combination of 2 elements from N
    for combo in combo_list:
        k=combo[0]
        ind=combo[1]
        repeat1=[]
        repeat2=[]
        for i,resi in enumerate(algn_map):
            # Check if the template alignment pair is valid and that both residues have a correspondence in the target
            if resi[k]!=('-','-') and resi[ind]!=('-','-') and resi[k] in res_map and res_map[resi[k]] in target_res_list and resi[ind] in res_map and res_map[resi[ind]] in target_res_list: #checks that they are coordinates rather than ['-']
                rmsd =rmsd+ pow(np.linalg.norm(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].get_vector().get_array()-s_mobile[0][res_map[resi[ind]][1]][res_map[resi[ind]][0]]['CA'].get_vector().get_array()),2)
                print("[",res_map[resi[k]][0],res_map[resi[ind]][0],"], ",end="")
                repeat1.append(np.array(list(s_mobile[0][res_map[resi[k]][1]][res_map[resi[k]][0]]['CA'].get_vector())))
                repeat2.append(np.array(list(s_mobile[0][res_map[resi[ind]][1]][res_map[resi[ind]][0]]['CA'].get_vector())))   
                n+=1             
        #print(repeat1, repeat2)
        tm_score=tm_score+tm_score_calc(repeat1,repeat2,prot_len)
        print(tm_score_calc(repeat1,repeat2,prot_len),n)
    print(tm_score)
    tm_score=tm_score*num_repeats/len(combo_list)
    rmsd=np.sqrt(rmsd/n)
    print(len(combo_list))
    print(num_repeats)
#    print("\nEstimated RMSD is: ",rmsd)
    print("\nEstimated TM-score is: ",tm_score)
    
    tmp=list(map(list, zip(*coord_repeats)))
    coord_repeats=tmp
    print("Estimated RMSD is: ",get_rmsd_cesymm(coord_repeats))
    
    print("Coord repeats",tmp)
    algn_map=[[] for i in range(positions)]
    for k in unique_repeats:              
        r=cesm[k]
        r_id=k.split('_')
        ch=r_id[0][-1]
        repeat=r_id[1].split('-')
    #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
        repeat_start=[i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch==sublist[2])] #finds the index of the residue that starts the repeat
    #    print(repeat_start)
        count=0
        for i,c in enumerate(r.seq): #Last edited
            if c!='/':
              if c!='-' and c.isupper() and resseq_A[repeat_start[0]+count][1] == c:
    #    		    print(c,resseq_A[repeat_start[0]+count][1],resseq_A[repeat_start[0]+count][0],"\n") 
                algn_map[i].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
                count+=1
     #           if resseq_A[repeat_start[0]+count][2] not in chains:
     #               chains.append(resseq_A[repeat_start[0]+count][2])
              if c=='-':
                algn_map[i].append(('-','-'))
              if c.islower():
                algn_map[i].append(('-','-'))
                count+=1     
    #print(algn_map)
    target_algn_map=[]
    for i,resi in enumerate(algn_map):
        # Check if the template alignment pair is valid and that both residues have a correspondence in the target
        target_algn_map.append([])
        for k in range(0,len(resi)):
            if resi[k]!=('-','-') and resi[k] in res_map and res_map[resi[k]] in target_res_list:
                target_algn_map[i].append(res_map[resi[k]])
            else:
                target_algn_map[i].append(('-','-'))
    #print(target_algn_map)
    target_repeats=['' for i in range(0,len(target_algn_map[0]))]
    count=0
    for j in range(0,len(target_algn_map[0])):
        flag=0
        for i, resi in enumerate(target_algn_map):
            if flag==0 and resi[j]!=('-','-') and resi.count(('-','-'))<=len(resi)-2:
                target_repeats[j]=inputf2+'_'+resi[j][1]+"_"+str(resi[j][0])
                flag=1
            if flag==1 and resi[j]!=('-','-') and resi.count(('-','-'))<=len(resi)-2:
                last=resi[j]
        target_repeats[j]=target_repeats[j]+'-'+str(last[0])
    print(target_repeats)
    target_repeats_all=[]
    target_repeats_type=repeats_type
    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    for ind, symm in enumerate(repeats_all): 
        target_repeats_all.append([])
        for i, sub_symm in enumerate(symm):               #sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            target_repeats_all[ind].append([])
            for k in sub_symm:
                pos=unique_repeats.index(k)
                target_repeats_all[ind][i].append(target_repeats[pos])
    print(target_repeats_all, repeats_type, repeats_levels, axes_per_level)
    print(transforms)
    # Build tree of transformations
    s=[np.array([i]) for i in reversed(range(0, len(transforms)))]
    all_transforms=[[]]
    for t in s:
        all_transforms_tmp=[]
        for i,ax in enumerate(all_transforms):
            all_transforms_tmp.append(list(ax)+list(t))
        all_transforms=all_transforms+all_transforms_tmp   

    all_transforms=all_transforms[1:]
    print("all transforms:", all_transforms)

    levels=[]
    half=len(all_transforms)+1
    for i in range(0,len(s)):
        subs=all_transforms[0:(half-1)]
        half=int((len(subs)+1)/2)
        spacing=len(subs)+1
        tmp=[subs[half-1]]
        jump=spacing
        while half+jump<=len(all_transforms):
            tmp.append(all_transforms[half-1+jump])
            jump+=spacing
        levels.append(tmp)       
    print(levels)
    
    # Calculate axes from tree
    all_axes=[]
    for ind,l in enumerate(levels):
        axes_tmp=[]
        for j in l:
            for i,t in enumerate(j):
                if i==0:
                    pt=target_sym_descriptors[t][3]
                    pt2=target_sym_descriptors[t][4]
                if i>0:
                    transl=target_sym_descriptors[t][2]
                    T=transforms[t]
                    pt=np.dot((pt-transl),T)
                    pt2=np.dot((pt2-transl),T)
            descriptor=[target_sym_descriptors[ind][0],target_sym_descriptors[ind][1],target_sym_descriptors[ind][2],pt,pt2,target_sym_descriptors[ind][4]]
            axes_tmp.append(descriptor)
        all_axes.append(axes_tmp)
    print(len(target_sym_descriptors),target_sym_descriptors,"\n")
    print(all_axes)


    # Method 2
    all_counts=[]
    print(order_of_level)
    for repeat in range(0, len(unique_repeats)):
        counts=[0]*len(transforms)            
        for i in reversed(range(0, len(counts))):
            d=int(order_of_level[i])
            counts[i]=repeat%d
            repeat=int(repeat/d)
        all_counts.append(counts)
    print(all_counts)  
    all_counts=all_counts[1:]       

    levels=[]
    half=len(all_counts)+1
    for i in range(0,len(s)):
        subs=all_counts[0:(half-1)]
        half=int((len(subs)+1)/2)
        spacing=len(subs)+1
        tmp=[subs[half-1]]
        jump=spacing
        while half+jump<=len(all_counts):
            tmp.append(all_counts[half-1+jump])
            jump+=spacing
        levels.append(tmp)       
    print(levels)
  
                
          
    #os.remove(oriented_opm)
    #if os.path.isfile(oriented_opm_2):
    #    os.remove(oriented_opm_2)
    return tm_score, rmsd

##########################################
### Main
if __name__ == '__main__':

    if len(sys.argv) < 3:
        raise SystemExit("syntax: %s template<pdbcode_chain> target<pdbcode.chain>" % sys.argv[0])
    else:
        arg1=sys.argv[1]
        arg2=sys.argv[2]
        tmscore,rmsd=transfer_sym(arg1,arg2)
   
