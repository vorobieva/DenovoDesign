import itertools
import re,sys,os,copy
from argparse import ArgumentParser
from copy import deepcopy
import numpy as np
import os
from BuildBlueprints import *
from Blueprint import Blueprint

parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-refpdb', type=str, help="input pdb")
parser.add_argument('-prefix', type=str)
parser.add_argument('-insertbp', type=str, help="blueprint for the insert. First and last residues must be remodeled.")
args = parser.parse_args()

template_xml = args.xml

cap = 0 #
# Include fixed beta-bulge but not tested beta-bulges in beta-strand residue counts!! 
#topol = "L[2-2]E[8-8]L[3-3]E[6-7]L[2-5]H[7-10]I[28-28]E[3-3]L[2-2]E[8-8]L[1-1]"
topol =  "L[2-2]E[8-8]L[3-3]E[6-7]L[2-5]H[7-10]I[28-28]E[3-3]L[2-2]E[8-8]L[1-1]"
common_bulges = {'E1':[2],'E2':[3]}
#common_bulges = {'E2':[3]}
#For now the script does not allow for simultaneous variation of strands lenght and bulge positions! 
tested_bulges = [{},{'E3': [4]}, {'E3': [2] }]

topology = [('E1','E2'),('E3','E4'),('E1','E4')]

cwd=os.getcwd()
#refpdb = os.path.join(cwd,args.refpdb)

ss,combinations = GetCombinations(topol)

print("There are %s possible combinations of secondary structure elements lengths." %(str(len(combinations))))
print(ss, combinations)

# Change this with a short code corresponding to your motif.
motif = "something"

insert_bp_data = []
with open(args.insertbp, "r") as in_bp:
	for line in in_bp:
		data = line.strip().split()
		insert_bp_data.append(data)



for bulge,comb in itertools.product(tested_bulges, combinations):
        # Build directories and add bulge combinations
        #---------------------------------------------------------
        # Make the directory name
        comb_name = ""
        filename = motif + '-'
        strand=0 ; resnum=0
        for i,s in enumerate(ss):
                comb_name +='%s%i' %(ss[i],comb[i])

        pathname=filename+comb_name

        bulges = deepcopy(common_bulges)
# First add the sampled bulge positions to the list of beta-bulges.
        for key, val in bulge.items():
                if key in bulges.keys():
                         bulges[key]+=bulge[key]
                else:
                         bulges[key] = bulge[key]
# Now add one residue to the strands that sample and additional beta-bulge.
                n_strands = 1
                found_strand = False
                for i in range(0,len(ss)):
                        if ss[i] == "E":
                                if "E"+str(n_strands) == key:
                                        comb_list = list(comb)
                                        comb_list[i] += 1
                                        comb = tuple(comb_list)
                                        found_strand = True
                                        print("Found strand for beta-bulge %s %s." %(i, n_strands))
                                        break
                                else:
                                        n_strands += 1
                if found_strand == False:
                        print("Something is wrong, the beta-strand corresponding to the beta-bulge %s %s was not found" %(key, val))
# Modify blueprint name to account for bulges.
        if bulges:
                keys = bulges.keys() # bulged strand names
                for k,key in enumerate(keys):
                        comb_name += '-b%s.%s' %(key,bulges[key]) # this is the new filename


#        blueprint_name = comb_name+ "_bp"
#        cst_name = comb_name + "_cst"
        if not os.path.exists(comb_name):
                os.makedirs(comb_name)
#for folder in A_* ; do cp A_1/bp $folder ; done


        blueprint_name = os.path.join(comb_name, comb_name+ "_bp")
        cst_name = os.path.join(comb_name, comb_name + "_cst")


        MakeRefBlueprint(ss,comb,bulges=bulges,refblue=blueprint_name,insert_bp=insert_bp_data, mostCommonLoopABEGO=True, strandsABEGO="B", pdbfile=args.refpdb)

        blue = Blueprint(blueprint_name)
        blue.reindex_blueprint(start=1)

#index variabl
        #motif index variable
        #loop through bps??
        #def findindex:
        #motifindex=[]
        #blue_file = open(trial_bp, 'r')
        #count=0
        #for line in blue_file:
        #        count +=1     
        #        if line[0:2] == "1 ":
        #                firstres= count
        #                lastres = firstres + 27
        #                variable = str(firstres) + "A" + "-" + str(lastres) + "A"
        #                motifindex.append(variable)
        #        else:
        #                continue
        #blue_file.close()
        #print ("motifindex:", motifindex)

#add to list
        # Make directory of the topology
#        if not os.path.exists(pathname):
#                os.mkdir(pathname)
#        os.chdir(pathname)

        cst_fileout = open(cst_name,'w')
        for hairpin in topology:
                if abs(int(hairpin[1][1])-int(hairpin[0][1])) == 1:
                # bulges can not happen between non-sequential strands for now
                        n_bulges = 0
                        bulge_to_add = None
                        if hairpin[0] in bulges.keys() and hairpin[0] != "E1":
                        # for now, it is only possible to add 1 bulge per hairpin, which should be sufficient for now if we consider the G1 bulge as part of the beta-turn
                                bulges[hairpin[0]].sort(reverse=True)
                                for b in bulges[hairpin[0]]:
                                        bulgepos = blue.segment_dict[hairpin[0]].bp_data[b-1][0]
                                        if (blue.segment_dict[hairpin[0]].bp_data[-1][0] - bulgepos - n_bulges) %2 != 0:
                                                n_bulges += 1
                                                bulge_to_add = bulgepos
                        elif hairpin[1] in bulges and hairpin[1] != "E1":
                                n_bulges = 0
                                bulges[hairpin[1]].sort()
                                for b in bulges[hairpin[1]]:
                                        bulgepos = blue.segment_dict[hairpin[1]].bp_data[b-1][0]
                                        if (bulgepos - blue.segment_dict[hairpin[1]].bp_data[0][0] - n_bulges) %2 == 0:
                                                n_bulges += 1
                                                bulge_to_add = bulgepos
                        elif hairpin[0] == "E1" and hairpin[1] == "E2":
                                bulgepos = 16
                                n_bulges += 1
                                bulge_to_add = bulgepos
                        print(hairpin,n_bulges,bulge_to_add)
                        if n_bulges == 1:
                                sthb = HbondsBulgedStrand(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,bulge_position=bulge_to_add) ; cst_fileout.write(sthb)
                        elif n_bulges == 0:
                                sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue) ; cst_fileout.write(sthb)
                        elif n_bulges > 1:
                                print("Something wrong with bulge count!")
                else:
			# Add code to place bulges between non-sequential strands
                        n_bulges = 0
                        bulges_to_add = None
                        if hairpin[0] in bulges.keys() and hairpin[0] != "E1":
                                bulges[hairpin[0]].sort(reverse=True)
                                print(hairpin[0])
                                for b in bulges[hairpin[0]]:
                                        bulgespos = blue.segment_dict[hairpin[0]].bp_data[b-1][0]
                                        if (blue.segment_dict[hairpin[0]].bp_data[-1][0] - bulgepos - n_bulges) %2 != 0:
                                                n_bulges += 1
                                                bulge_to_add = bulgepos
                        elif hairpin[1] in bulges and hairpin[0] != "E1":
                                n_bulges = 0
                                bulges[hairpin[1]].sort()
                                for b in bulges[hairpin[1]]:
                                        bulgepos = blue.segment_dict[hairpin[1]].bp_data[b-1][0]
                                        if (bulgepos - blue.segment_dict[hairpin[1]].bp_data[0][0] - n_bulges) %2 == 0:
                                                n_bulges += 1
                                                bulge_to_add = bulgepos
                        elif hairpin[0] == "E1" and hairpin[1] == "E4":
                                bulgepos = 4
                                n_bulges += 1
                                bulge_to_add = bulgepos
                        print(hairpin,n_bulges,bulge_to_add)
                        if n_bulges == 1:
                                sthb = HbondsBulgedStrand(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,bulge_position=bulge_to_add,str1_offset=0) ; cst_fileout.write(sthb)
                        elif n_bulges == 0:
                                sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,str1_offset=0); cst_fileout.write(sthb)
        cst_fileout.close()


#output directory generator
        #tasklist=open("tasklist", "w")
        #tasklist.close()
        x="0123456789"
        for i in (x):
                stri = str(i)
                label= comb_name + "_" + stri
        #print(label)
                if not os.path.exists(label):
                        os.makedirs(label)
                folder_name = os.path.join(comb_name, label)

                labelstr = str(label)
                labelbp = comb_name  + "_bp"
                labelbpstr = str(labelbp)
                labelcst = comb_name + "_cst"
                labelcststr = str(labelcst)
                #print (line)
                tasklistappend = open("tasklist", "a")
                bp_path = "/work/upcorreia/users/karla/bpgenerator/" + comb_name + "/" + labelbp
                cst_path = "/work/upcorreia/users/karla/bpgenerator/" + comb_name + "/" + labelcst
                strbppath = str(bp_path)
                strcstpath = str(cst_path)
                tasklistappend.write(labelstr + " " + strbppath + " " + strcstpath + "\n")
        tasklistappend.close()

################3
#tasklist2 with motifindex
       # motifindex=[]
       # with open("tasklist2", 'w') as task2:
       #         with open("tasklist", 'r') as task1:
       #                 count1=0
       #                 for taska in task1: 
       #                         count1 += 1
       #                         bp_filepath = taska.split(" ")[1]
       #                         with open(bp_filepath, 'r') as taskbp:
       #                         #motifindex=[]
        #                                count2=0
        #                                for line in taskbp:
        #                                        count2 +=1
        #                                        if line[0:2] == "1 ":
        #                                                firstres= count2
        #                                                lastres = firstres + 27
        #                                                variable = str(firstres) + "A" + "-" + str(lastres) + "A"             
        #                                                motifindex.append(variable)
        #                                                for taska in task1:
        #                                                        task2.write(taska.rstrip(' \n') + variable + '\n')
                                               # else:
                                                #        continue
                                                #for taska in task1:
                                                #task2.write(taska.rstrip('\n') + " " +  variable + '\n')
       # print("motifindex: ", motifindex)                                      
 # else:
                                        #        continue
                                               # for task in task2:

                                                #         line = line.rstrip('\n') + variable                        
##########################
#        xml_lines = open('../%s' %(template_xml),'r').readlines()
