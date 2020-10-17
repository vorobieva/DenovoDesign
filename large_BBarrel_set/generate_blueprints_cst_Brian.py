import itertools
import re,sys,os,copy
from argparse import ArgumentParser
from copy import deepcopy
import numpy as np
import os
from BuildBlueprints import *
from Blueprint import Blueprint

#==============================
# INPUT PARAMETERS
#==============================
parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-refpdb', type=str, help="input pdb")
parser.add_argument('-prefix', type=str)
args = parser.parse_args()

template_xml = args.xml

# topol is the barrel topology only. The first beta-sheet residue is one residue of the tryptohane/tyrosine corner (Gly or Ser). The tryptophane or tyrosine of the corner is the THIRD residue of the FIRST beta-strand. It is assumed that the LYS/ARG residue on the last strand is the second residue starting from the end of the last strand.
# if there is a bulge in E2, then it should be 13 residues long. If non-bulged, then 12
# cap is the capping motif to design. Currently, the classic tryptophan corner is the only supported motif. The script should further check that there is no bulge associated to the Gly of the corner (first residue of the first strand).
# cap = 0 no cap
# cap = 1 is for classic tryptophan corner
# cap = 2 will be bulge-associated corner.
# test topol 1: "E[10,12]L[3-5]E[14,16]L[2,4,7-9]E[10-12]L[3-3]E[12-14]L[2-2]E[8-8]L[3-3]E[10-10]L[2-2]E[8-8]L[3-3]E[10-10]L[1-1]"

cap = 1 # capping motif -- 1 corresponds to trp corner
topol = "E[10-12]L[3-5]E[12-14]L[2]E[10]L[3-5]E[12]L[2]E[8-10]L[3-5]E[10-12]L[2]E[8-10]L[3-5]E[12-14]L[1-1]"
nStrand = 8
shearNo = 10
max_kinks = 8 # the cutoff must be between 5 (one kink per Cbeta strip) and 10 (2 kinks per Cbeta strip).
min_kinks = 6
common_bulges = {'E2':[-2],'E4':[-2],'E6':[-2]}
topology = [('E1','E2'),('E2','E3'),('E3','E4'),('E4','E5'),('E5','E6'),('E6','E7'),('E7','E8'),('E8','E1')]

# All bulges located in even-numbered strands. 
# The bulges are fixed in the bottom turns, at position -2 of the longer even-numbered strands. Both bulges and kinks defined from the bottom of the barrel and describe the starting model. 

cwd=os.getcwd()
refpdb = os.path.join(cwd,args.refpdb)

ss,combinations,bulges = getCombinationsForBetaBarrels(topol,common_bulges)

print "There are %s possible combinations of secondary structure elements lengths... Generating blueprints and cst files with no glycine kinks." %(str(len(combinations)))
print ss, combinations, bulges
orderedCBetaStrips = findCbetaStrips(nStrand, shearNo, topol)
possible_kink_combinations = calculate_kink_combinations(orderedCBetaStrips, max_kinks, min_kinks)
print "The number of allowed kinks combinations is " + str(len(possible_kink_combinations))

# account for length differences due to bulges 
# this is special case stuff for our b barrels that works because bulges on the "bottom" of the barrel
# are guaranteed to be in the postions (for now)
# feel free to generalize.
for cb_strip in orderedCBetaStrips:
	for i, (strandNo, pos) in enumerate(cb_strip):
		if strandNo in [1, 3, 5]:
			cb_strip[i] = (strandNo, pos - 1)
		# This is a silly fix to the mistake that we made while indexing the kinks on the odd strands on base 1 instead of 0 as in the blueprint. Hopefully I'll
		# have time to fix this later. Omit strand 0 since we shortened it by 1 residue to place the tryptophan corner. 
		if strandNo in [2,4,6]:
			cb_strip[i] = (strandNo, pos - 1)
		if strandNo == 7:
			cb_strip[i] = (strandNo, pos + 2)

for bulge, comb in zip(bulges, combinations):
	comb_name = ""
	for i,s in enumerate(ss):
		comb_name+='%s%i' %(ss[i],comb[i])
	print comb_name
	ss_with_corner, comb_with_corner = add_tryptophan_corner(ss, comb)
	blueprint_name = '_bp'	
	cst_name = comb_name+'_cst'

	MakeRefBlueprint(ss_with_corner,comb_with_corner,refblue=blueprint_name, bulges=bulge, mostCommonLoopABEGO=True, strandsABEGO="B")

	blue = Blueprint(blueprint_name)
	blue.reindex_blueprint(start=1)

	cst_fileout = open(cst_name,'w')
	for hairpin in topology:
		if abs(int(hairpin[1][1])-int(hairpin[0][1])) == 1:
			# bulges can not happen between non-sequential strands for now
			n_bulges = 0
			bulge_to_add = None
			if hairpin[0] in bulge:
				# for now, it is only possible to add 1 bulge per hairpin, which should be sufficient for now if we consider the G1 bulge as part of the beta-turn
				bulge[hairpin[0]].sort(reverse=True)
				for b in bulge[hairpin[0]]:
					bulgepos = blue.segment_dict[hairpin[0]].bp_data[b-1][0]
					if (blue.segment_dict[hairpin[0]].bp_data[-1][0] - bulgepos - n_bulges) %2 != 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			elif hairpin[1] in bulge:
				n_bulges = 0
				bulge[hairpin[1]].sort()
				for b in bulge[hairpin[1]]:
					bulgepos = blue.segment_dict[hairpin[1]].bp_data[b-1][0]
					if (bulgepos - blue.segment_dict[hairpin[1]].bp_data[0][0] - n_bulges) %2 == 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			if n_bulges == 1:
				sthb = HbondsBulgedStrand(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,bulge_position=bulge_to_add) ; cst_fileout.write(sthb)
			elif n_bulges == 0:
				sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue) ; cst_fileout.write(sthb)
			elif n_bulges > 1:
				print "Something wrong with bulge count!"

		else:
			if  cap == 1:
				sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,str1_offset=1); cst_fileout.write(sthb)


	corner = cornerConstraints(capType=cap,blueprint=blue); cst_fileout.write(corner)

	cst_fileout.close()

	# we no longer need the numbering in the blue print; remove it so it we can build the structure
	blue.remodel_all()
	blue.bp_data[0] = [1, 'A', 'L', '.']
	# this where the glycine in the trp corner goes
	cap_gly = (0, 0) if cap == 1 else None 
	for i, kink_comb in enumerate(possible_kink_combinations):
		possible_combo = True
		strand_kink_pos = []
		for cbeta_strip, pos in enumerate(kink_comb):
			strand_kink_pos += [orderedCBetaStrips[cbeta_strip][p] for p in pos]       
		if cap_gly is not None:
			strand_kink_pos = [cap_gly] + strand_kink_pos
		# make a deepcopy of the blue print
		bp = deepcopy(blue)
		for strandNo, strandPos in strand_kink_pos:
			try:
				bp.segment_dict["E" + str(strandNo + 1)].bp_data[strandPos][1:3] = ["G", "EE"]
			except IndexError:
				possible_combo = False
		# make a directory for the blueprint
		if possible_combo == True:
			dirname = os.path.join(comb_name, str(i))
			os.makedirs(dirname)
			bp.dump_blueprint(os.path.join(dirname, comb_name + "_{:04d}_bp".format(i)))
		# write XML file here??
		if possible_combo == False:
			print strand_kink_pos
	os.remove(blueprint_name)
       
#
#

