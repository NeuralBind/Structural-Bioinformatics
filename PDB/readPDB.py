#!/usr/bin/python3
"""
Structural Bioinformatics assignment 1 - Phi-psi angles

The script should create an
output file containing the phi and psi angles and secondary
structure assignment of each residue. 

To run, make 'PDB' your working directory and use:
    
cd ~/VU/ 2.Structural\ bioinformatics/practical/ass1/Assignment1_Files_2021/PDB

> python3 readPDB.py pdb_filename.txt

./readPDB.py <filename>.pdb

Author: Theodoros Foskolos

"""

# Packages
from sys import argv
import os
from math import sqrt, atan2, degrees

# Vector functions that we need to calculate the angles
def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
#print(dot_product([1, 2, 3], [1, 3, 2]))

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1]*v2[2] - v1[2]*v2[1]
    j = v1[2]*v2[0] - v1[0]*v2[2]
    k = v1[0]*v2[1] - v1[1]*v2[0]
    return [i,j,k]
#print(cross_product([1, 2, 3], [1, 3, 2]))

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)
#print(magnitude([1, 2, 2]))

# PDB file parser
def readPDB(PDB_file):
    """ Reads a PDB file and stores the atom
    coordinates and amino acid types of the protein """
    # open the file
    f = open(PDB_file, 'r')

    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}

    # parse each line in the file
    for line in f:
        # remove whitespace at the end of the line
        line = line.strip()
        # only parse the lines containing atom coordinates
        if line[:4] == 'ATOM':
            # ATOM type (e.g. C-alpha)
            atom_type = line[12:16].strip()
            # AMINO ACID type (e.g. alanine)
            aa_type = line[17:20].strip()
            # residue number
            res_num = int(line[22:26])
            # Protein chain
            chain = line[21]
            # coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            # if chain does not exists create new entry
            if not chain in pdbcoord:
                pdbcoord[chain] = {}
                pdbseq[chain] = {}
            # if resnum does not exists create new entry
            if not res_num in pdbcoord[chain]:
                pdbcoord[chain][res_num] = {}

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord,ycoord,zcoord]
            # store sequence
            pdbseq[chain][res_num] = aa_type

    # close file
    f.close()

    # return dictionaries
    return pdbcoord, pdbseq
#print(readPDB('Data\\1TIM.pdb'))

### THE FOLLOWING THREE FUNCTIONS ARE THE ONES YOU NEED
### TO EDIT FOR THE ASSIGNMENT. ONLY EDIT THE INDICATED
### BLOCKS
def calculateDihedral(a1, a2, a3, a4):
    """ Calculates the normal vector of the planes
    defined by four atom coordinates """

    ### START CODING HERE
    # calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined above
    # you can also use the python math function "atan2" and "degrees"
    
    

    # calculate connecting vectors (bonds)
    
    # Plane 1
    # 1st we need to calclate the bonds in order to get the normal vectors of each plane and then calculate the angle by finding sin and cos
    
    conn_vect_1 = [a2[i] - a1[i] for i in range(3)]
    conn_vect_2 = [a3[i] - a1[i] for i in range (3)]
    
    # normal vector 1
    norm_vector_1 = cross_product(conn_vect_1, conn_vect_2)
    
    # Plane 2
  
    conn_vect_3 = [a3[i] - a2[i] for i in range (3)]
    conn_vect_4 = [a4[i] - a2[i] for i in range (3)]  
    
    #calculate unit vector
    
    u = [conn_vect_3[i]/ magnitude(conn_vect_3) for i in range (3)]
    
    
    #normal vector 2
    norm_vector_2 = cross_product(conn_vect_3, conn_vect_4)
    
    # Dot product of the 2  normal vectors of the 2 planes using the functions
    
    Dot_prod = dot_product(norm_vector_1, norm_vector_2)
    
    # Cross product of the 2 normal vectors of the 2 planes using the functions
    
    Cross_prod = cross_product(norm_vector_1, norm_vector_2)
    
    # Magnitude multiplication
    
    Magn_mult = magnitude(norm_vector_1) * magnitude(norm_vector_2)

    # Calculate the cos, sin of the angle φ
    
    cos_angle = Dot_prod / Magn_mult
    
    #for the sin i needed to calculate the cross product devided by the magnitude of the 2 norm vectors also devided by the unit one
    
    sin_angle = Cross_prod[1] / (Magn_mult * u[1])
    
        
    # Calculate the angle
    angle = atan2(sin_angle, cos_angle)
    
    # Pass it in dihedral and transform it to angle
    dihedral = degrees(angle) # replace this line with your code
    ### END CODING HERE
    return dihedral

#print(calculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))

def assign_ss(phi, psi):
    """ Assign a secondary structure type based on the phi
    and psi angles of a residue """
    ### START CODING HERE
    
    secondary_structures = None
    # we need to ensure we are not on last or 1st residue in order to have a secondary structure   
    if phi is not None and psi is not None:
        
        # According to the ramachadran plot these are some approach limits for the structures
        if (-170 < phi < -20) and (-90 < psi < 50)  : #there are some right handed alpha helixes but im not sure if needed
            secondary_structures = "alpha"
    
        elif (-180 < phi < -20) and (20 < psi < 180):
            secondary_structures = "beta"
    
        else:
            #otherwise its a loop/coil
            secondary_structures = "loop"
        
       
    # for code checking purposes use the terms "loop", "alpha" or "beta"
    ### END CODING HERE
    return secondary_structures
#print(assign_ss(60, 25))



def print_phi_psi(pdbcoord, pdbseq, outfile):
    """ given the PDB coordinates, calculate the dihedral
    angles of all the residues, assign secondary structure
    types and write them into an output file """
    f = open(outfile, 'w')

    # get the chains from the PDB file
    list_chains = sorted(pdbcoord.keys())

    for chain in list_chains:
        # get the sorted residue numbers from the pdbcoord dictionary
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in list_residue_numbers:
            # if certain residues are missing in the PDB file, you will
            # get a KeyError. Make sure your program does not crash, but
            # gives a warning when this happens
            try:
                ### START CODING HERE
                phi,psi,ss = None, None,None
                phi_flag,psi_flag = True, True
                
                # The 1st residue dont have a phi only psi because it would need to calculate the Carbon of the previous one , whereas there is no previous residue
                # So calculate only psi
                
                # The last residue wont have a psi angle only a phi because we dont have a next residue Nitrogen to calculate
                # So calculate only the phi
                
                if res_num == list_residue_numbers[0]: # 1st residue
                    phi_flag = False
                
                if res_num == len(list_residue_numbers): # last one
                    psi_flag = False
                
                
                if psi_flag:              
                    # for the psi (ψ) is the N(i),Ca(i),C(i),N(i+1) torsion angle
                    N = pdbcoord[chain][res_num]["N"]
                    CA = pdbcoord[chain][res_num]["CA"]
                    C = pdbcoord[chain][res_num]["C"] 
                    N_after =  pdbcoord[chain][res_num+1]["N"]
            
                    psi = calculateDihedral(N, CA, C, N_after)
                
                if phi_flag:
                    # for the  phi (φ) i need the C(i-1),N(i),Ca(i),C(i) torsion angle
                    C_prev = pdbcoord[chain][res_num-1]["C"] 
                    N = pdbcoord[chain][res_num]["N"]
                    CA = pdbcoord[chain][res_num]["CA"]
                    C =  pdbcoord[chain][res_num]["C"]
            
                    phi = calculateDihedral(C_prev, N, CA, C)
                
                ss = assign_ss(phi,psi)
                ### END CODING HERE
            
            except KeyError:
                print('WARNING: KeyError:', KeyError, 'in residue', chain, res_num)
                
            # get amino acid
            aa_type = pdbseq[chain][res_num]
            # write into output file
            print(chain, res_num, aa_type, phi, psi, ss, file=f)
            
    f.close()
    print('written:', outfile)

def main():
    # input PDB file
    f_in = argv[1]
    f_out = 'Output/phi_psi.txt'
    # read PDB file
    pdbcoord, pdbseq = readPDB(f_in)
    print_phi_psi(pdbcoord, pdbseq, f_out)
    
    
   
    # for testing
    #for i in ['1TIM', '3PG8']:
     #    f_in = 'student/{}.pdb'.format(i)
      #   print(f_in)
       #  f_out = 'student/output/phi_psi_{}.txt'.format(i)
    
          # read PDB file
        # pdbcoord, pdbseq = readPDB(f_in)
        # print_phi_psi(pdbcoord, pdbseq, f_out)

if __name__ == '__main__':
    main()
