#! /usr/bin/env python
""""this script use generate lines as described in Laster's paper to build hyperboloid of one sheet"""
# in this script, the residue coordinates array is populated along the strand, starting from a reference residue


from math import sqrt, sin, cos, atan2, pi
from scipy.optimize import fsolve, fmin
import argparse
import string
import numpy as np


parser = argparse.ArgumentParser(description='Program')
parser.add_argument('n', action="store", type=int, help="Number of beta strands of beta-barrel")
parser.add_argument('S', action="store", type=int, help="Shear number of beta-barrel")
parser.add_argument('nres', action="store", type=int, help="Number of residues in each strand")
parser.add_argument('dtw', action="store", type=float, default=0, help="twist angle difference relative to the opitmal")
parser.add_argument('b', action="store", type=float, default=4.4, help="inter strand distance")
parser.add_argument('da', action="store", type=float, default=0, help="radius difference relative to the optimal")
parser.add_argument('db', action="store", type=float, default=0, help="radius difference relative to the optimal")
parser.add_argument('topology', action="store", default='1,2,5,4,3,6', help="topology")
parser.add_argument('--a', action="store_true", help="Anti-parallel")

# distance between ca along the reference vector
a = 3.3
# topology = [1,2,5,4,3,6] #Greek key
PI = pi
args = parser.parse_args()
topology_str = args.topology  # Greek key
topology = topology_str.split(',')
topology = [int(c) for c in topology]
n = args.n
S = args.S
nres = int(args.nres)
b = args.b
da = args.da
db = args.db
dtw = args.dtw

# r=sqrt((n*b)**2+(S*a)**2)/(2*n*sin(PI/n)) #the radius of beta barrel
theta = atan2((S * a), (n * b))  # the tilt angle of beta strand
# r=sqrt((n*b)**2+(S*a)**2)/(2*PI)
# theta=asin(S*a/(2*PI*r))


def twistAndCoil(x):
    return (x - 20*PI / 180) ** 2 + ((2*PI * sin(theta) - n * x) / S) ** 2 + ((2*PI * cos(theta) - S * a * x / b) / n) ** 2


tw1 = fmin(twistAndCoil, 0)[0]  # twist angle between strand
c1 = (2 * PI * sin(theta) - n * tw1) / S  # coil angle within strand
c2 = (2 * PI * cos(theta) - S * a * tw1 / b) / n  # coil angle between strand
for x in (theta, tw1, c1, c2):
    print x * 180 / pi,
# rho = (2*tw1*sin(theta)*cos(theta)/b-(c1/a)*(cos(theta)**2)-(c2/b)*(sin(theta)**2))/2
# build 2 extra residues on each strand terminii to help BBQ building the backbone
nres += 4
theta *= dtw
# C = sqrt((b ** 2 * (1 + cos(tw1)) / (4 * (1 - cos(tw1)))))
# r = sqrt(b ** 2 * (1 + cos(tw1)) / (4 * (cos(tw1) - cos(2 * PI / n))))
r = sqrt((S * a) ** 2 + (n * b) ** 2) / (2 * n * sin(PI/n))
A = r * da
B = r * db
C = (cos(theta)/sin(theta)) * sqrt(A**2 * sin(2 * PI / n)**2 + B**2 * cos(2 * PI / n)**2)


print "A=%f,B=%f,C=%f,twist=%f,tilt=%f\n" % (A, B, C, tw1 * 180 / PI, theta * 180 / PI)


def main():
    # get the metrics and references characterizing each strand
    strand_reference_point, strand_reference_vector = find_strand_reference_points(A, B, C, n, PI)
    barrel_residues_coordinates = list()
    # define the sign of the reference residue

    # find the other residues of a strand based on the position of the reference residue
    for strand in range(1, n + 1):
        dresi = find_displacement_per_residue(a, strand_reference_vector[strand - 1])
        # define the sign of the reference residue
        sign0 = define_reference_residue_sign(strand)
        reference_resi = strand_reference_point[strand-1]
        reference_vector = strand_reference_vector[strand-1]
        # create the array t store the coordinates of strands (8 + 4 residues for 1-4 and 6 + 4 residues for 5-8)
        if 1 <= strand <= 7:
            coor_strand = np.empty((8, 3))
        elif strand == 8:
            coor_strand = np.empty((10, 3))
        number_residues = coor_strand.shape[0]
        residue_index_list = get_residue_index(strand, number_residues)
        for residue in range(0, number_residues):
            index = residue_index_list[residue]
            coor_strand[residue] = arrange_residues_along_strand(index, dresi, reference_resi, reference_vector, sign0)
        barrel_residues_coordinates.append(coor_strand)
    write_pdb_file(barrel_residues_coordinates)


def pos(t, v, x0, y0, z0, reference_pt, strand):
    # coordinate of points on strand i
    xt, yt, zt = x0 + v[0] * t, y0 + v[1] * t, z0 + v[2] * t
    # the vector of previous reference point to the current vector
    vec_ref_pt = (xt - reference_pt[strand - 2][0], yt - reference_pt[strand - 2][1], zt - reference_pt[strand - 2][2])
    # inner product of vector p(i,t0)-p(i-1,t0) and direction vector
    out = vec_ref_pt[0] * v[0] + vec_ref_pt[1] * v[1] + vec_ref_pt[2] * v[2]
    return out


def get_residue_index(strand, number_residues):
    # now find the coordinates of the other residues along the strand
    residue_index_list = np.empty(number_residues)
    # place the residues coordinates in the array
    if strand == 1 or strand == 3:
        # the reference residue is at position 5 out of 8. Index the residues with the middle residue as a reference
        for j in range(0, number_residues):
            residue_index = j - 4
            residue_index_list[j] = residue_index
    elif strand == 5 or strand == 7:
        # the reference residue is at position 3 out of 8.
        for j in range(0, number_residues):
            residue_index = j - 2
            residue_index_list[j] = residue_index
    elif strand == 2 or strand == 4:
        # the reference residue is at position 4 out of 8. The direction of the strand is reversed
        for j in range(0, number_residues):
            residue_index = 3 - j
            residue_index_list[j] = residue_index
    elif strand == 6 or strand == 8:
        # the reference residue is in position 6 out of 8 or 10. the direction of the strand is reversed
        for j in range(0, number_residues):
            residue_index = 5 -j
            residue_index_list[j] =residue_index
    return residue_index_list


def find_strand_reference_points(A, B, C, n, PI):
    """

    :rtype : returns two numpy arrays. The coordinate of the reference residue for each strand and
    the coordinates of the reference vector.
    """
    strand_reference_point = np.empty([8, 3])
    strand_reference_vector = np.empty([8, 3])

    for strand in range(1, n + 1):
        gyration_angle = 2 * PI * (1 - strand) / n
        x0, y0, z0 = A * cos(gyration_angle), B * sin(gyration_angle), 0
        # direction vector of two skewer lines family
        v = [-A * sin(gyration_angle), B * cos(gyration_angle), C]
        strand_reference_vector[strand - 1] = [v[0], v[1], v[2]]
        #    u = (A*sin(gyration_angle),-B*cos(gyration_angle),C)
        # find the start points for strands other than the first one
        dt = find_displacement_per_residue(a, strand_reference_vector[strand - 1])
        if strand > 1:
            # get the ca coordinate of strand strand+1 in register with reference ca of strand strand
            ti = fsolve(pos, -0.5, args=(v, x0, y0, z0, strand_reference_point, strand))[0]
            # define the residue just below this one as the new reference for strands 2-4
            if strand == 5:
                ti -= 3 * dt
            # go three residue below for strands 5-8
            else:
                ti -= dt
            # re-define the new reference points for strands 2-8
            x0, y0, z0 = x0 + v[0] * ti, y0 + v[1] * ti, z0 + v[2] * ti
        # store all reference residues in a numpy array
        strand_reference_point[strand - 1] = [x0, y0, z0]
    return strand_reference_point, strand_reference_vector


def define_reference_residue_sign(strand):
    if strand % 2 == 0:
        sign0 = 1
    else:
        sign0 = -1
    return sign0


def find_displacement_per_residue(a, vector_coor):
    print vector_coor
    sum_vector_coor = 0
    for coordinate in vector_coor:
        sum_vector_coor += coordinate ** 2
    vector_norm = sqrt(sum_vector_coor)
    displacement_per_residue = a / vector_norm  # since (xt,yt,zt) = (x0,y0,z0)+vt
    return displacement_per_residue


def arrange_residues_along_strand(index, dresi, reference_residue, reference_vector, sign0):
    t = index * dresi
    x0, y0, z0 = reference_residue[0], reference_residue[1], reference_residue[2]
    # define ca positions along the axe based on the reference ca.
    xt, yt, zt = x0 + reference_vector[0] * t, y0 + reference_vector[1] * t, z0 + reference_vector[2] * t
    sign = sign0 * ((-1) ** abs(index))

    # move the ca position up and down in zigzag. The reference ca on strand 1 is down (sign = -1)
    vn = 2 * xt / (A ** 2), 2 * yt / (B ** 2), -2 * zt / (C ** 2)
    # the normal vector of point at the surface
    sum_vector_vn_coor = 0
    for coordinate in vn:
        sum_vector_vn_coor += coordinate ** 2
    vector_vn_norm = sqrt(sum_vector_vn_coor)
    dtn = sign * sqrt(1.9 * 1.9 - 1.65 * 1.65) / vector_vn_norm

    xn, yn, zn = xt + vn[0] * dtn, yt + vn[1] * dtn, zt + vn[2] * dtn
    coor_residue = [xn, yn, zn]
    return coor_residue


def write_pdb_file(barrel_residues_coordinates):
    outfile = 'test_ca.pdb'
    chain_ids = string.uppercase
    # chain_ids = 'A'*26
    fo = open(outfile, 'w')
    residue_id = 1
    atom_id = 2
    for i in range(len(barrel_residues_coordinates)):
        coor_strand = barrel_residues_coordinates[topology[i] - 1]
        for j in range(len(coor_strand)):
            (x, y, z) = coor_strand[j]
            fo.write("ATOM  %5d  CA  ALA %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (
            atom_id, chain_ids[i], residue_id, x, y, z))
            residue_id += 1
            atom_id += 4

    fo.close()

'''
cst_file = 'all.cst'
fo = open(cst_file,'w')
nres -= 4
mid = []
hbond = {}
sum_shift = 0
for i in range(n):
    mp = int((nres+1)/2)+nres*(topology[i]-1)
    if i<n-1:
	mpn = mp+nres*(topology[i+1]-topology[i])+shifts[i+1]*((-1)**i)
	sum_shift += shifts[i+1]*((-1)**i)
    else:
	mpn = mp+nres-sum_shift #find the registered residue of the last strand
    mid.append((mp,mpn))
    shift_d = -1*int((nres-1)/2)
    shift_u = int(nres/2)
    for j in range(shift_d,shift_u+1):
   # for j in range(-2,3):
	res1 = mp + j
	res2 = mpn - j
	if i==n-1: res2 -= n*nres #handle the last strand
	if res2<=nres*(topology[(i+1)%n]-1) or res2>nres*(topology[(i+1)%n]):
	#    print res1,res2,topology[i],topology[(i+1)%n]
	    continue
	up_b = 5.5+0.1*abs(j)
	tol = 0.5+0.1*abs(j)
	tolh = 0.1+0.1*abs(j)
#	fo.write("AtomPair CA %3d CA %3d BOUNDED 4.0 %.2f %.2f\n" %(res1,res2,up_b,tol))
	if i==0:
	    if j%2 == 0:
		fo.write("AtomPair O %3d N %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
		fo.write("AtomPair N %3d O %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
		hbond[res2]=res1
	else:
	#    print mp,j,res1,res2
	    if not hbond.has_key(mp):
		if j%2 == 0:
		    fo.write("AtomPair O %3d N %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
		    fo.write("AtomPair N %3d O %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
		    hbond[res2]=res1
	    else:
		if j%2 == 1:
		    fo.write("AtomPair O %3d N %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
                    fo.write("AtomPair N %3d O %3d BOUNDED 2.5 3.5 %.2f NOE;\n" %(res1,res2,tolh))
		    hbond[res2]=res1
fo.close()
print mid

bpfile = 'test.blueprint'
fo = open(bpfile,'w')
loop_l = 2
resi = 0
for i in range(n):
    for j in range(nres):
	resi += 1
	fo.write("%-d A E\n" %(resi,))
    if i<n-1:
	for k in range(loop_l):
	    fo.write("0 x L ALLAA\n")
fo.close()
'''
main()
