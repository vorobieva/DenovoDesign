#! /usr/bin/env python
'this script use generate lines as described in Laster paper to build hyperboloid of one sheet'
from math import sqrt, asin, sin, tan, cos, atan, pi, floor
from scipy.optimize import fsolve, fmin
import argparse
import string
import sys
import numpy as np
parser = argparse.ArgumentParser(description='Program')
parser.add_argument('n', action="store", type=int, help="Number of beta strands of beta-barrel")
parser.add_argument('S', action="store", type=int, help="Shear number of beta-barrel")
parser.add_argument('nres', action="store", type=int, help="Number of residues in each strand")
parser.add_argument('dtw', action="store", type=float,default=0, help="twist angle difference relative to the opitmal")
parser.add_argument('b', action="store", type=float,default=4.4, help="inter strand distance")
parser.add_argument('da', action="store", type=float,default=0, help="radius difference relative to the opitmal")
parser.add_argument('db', action="store", type=float,default=0, help="radius difference relative to the opitmal")
parser.add_argument('topology', action="store", default='1,2,5,4,3,6', help="toplogy")
parser.add_argument('--a', action="store_true", help="Anti-parralel")

a = 3.8
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
theta = atan(S*a/(n*b))  # the tilt angle of beta strand
# r=sqrt((n*b)**2+(S*a)**2)/(2*PI)
# theta=asin(S*a/(2*PI*r))


def twistAndCoil(x):
    return (x-20*PI/180)**2+((2*PI*sin(theta)-n*x)/S)**2+((2*PI*cos(theta)-S*a*x/b)/n)**2


tw1 = fmin(twistAndCoil, 0)[0]  # twist angle between strand
c1 = (2*PI*sin(theta)-n*tw1)/S  # coil angle within strand
c2 = (2*PI*cos(theta)-S*a*tw1/b)/n  # coil angle between strand
for x in (theta, tw1, c1, c2):
    print x*180/pi,
# rho = (2*tw1*sin(theta)*cos(theta)/b-(c1/a)*(cos(theta)**2)-(c2/b)*(sin(theta)**2))/2
nres += 4
tw1 *= dtw
C = sqrt((b**2*(1+cos(tw1))/(4*(1-cos(tw1)))))
r = sqrt(b**2*(1+cos(tw1))/(4*(cos(tw1)-cos(2*PI/n))))
A = r*da
B = r*db

print "A=%f,B=%f,C=%f,twist=%f,tilt=%f\n" % (A, B, C, tw1*180/PI, theta*180/PI)


def vecLength(vec):
    sum = 0
    for x in vec:
        sum += x**2
    return sqrt(sum)


def pos(t):
    xt, yt, zt = x0+v[0]*t, y0+v[1]*t, z0+v[2]*t  # coordinate of points on strand i
    vec_refpt = (xt-ref[i-2][0], yt-ref[i-2][1], zt-ref[i-2][2])  # the vector of previous reference point to the current vector
    out = vec_refpt[0]*v[0]+vec_refpt[1]*v[1]+vec_refpt[2]*v[2]  # inner product of vector p(i,t0)-p(i-1,t0) and direction vector
    return out

strand_reference_point = np.zeros((3, 8))
strand_reference_vector = np.zeros((3, 8))
sign0 = 1  # the sign of middle residue
# find the reference position for each strand (intersection between that strand and the central ellipse)
for i in range(1, n+1):
    omega = 2*PI*(i-1)/n
    x0, y0, z0 = A*cos(omega), B*sin(omega), 0
    strand_reference_point[i-1] = [x0, y0, z0]
    # direction vector of two skewer lines family
    v = (-A*sin(omega), B*cos(omega), C)
    strand_reference_vector[i-1] = [v[0], v[1], v[2]]
#    u = (A*sin(omega),-B*cos(omega),C)
    dt = a/vecLength(v)  # since (xt,yt,zt) = (x0,y0,z0)+vt
    # find the start points for strands other than the first one

# create a numpy array to stack the point coordinates, per strand

coord_forward = list()

for move_forward in range(0, 11):
    

    if i > 1:
	ti = fsolve(pos, -0.5)[0]  # get the value of parameter t which corresponding to the ca coordinate of strand i+1 in register with reference ca of strand i

	x0, y0, z0 = x0+v[0]*ti, y0+v[1]*ti, z0+v[2]*ti

    ref.append([x0, y0, z0])
    coor_strand = []
#    print 'sign0=',sign0
    for j in range(nres):
	resi = int(floor((j-(nres-1)/2.0)))
	t = resi*dt
	xt, yt, zt = x0+v[0]*t, y0+v[1]*t, z0+v[2]*t
	sign = sign0*((-1)**abs(resi))
    # print resi, sign
    # sign = 1 if j%2 == 0 else -1
    # '''
	vn = 2*xt/(A**2), 2*yt/(B**2), -2*zt/(C**2)  # the normal vector of point at the surface
	dtn = sign*sqrt(1.9*1.9-1.65*1.65)/vecLength(vn)
	xn, yn, zn = xt+vn[0]*dtn, yt+vn[1]*dtn, zt+vn[2]*dtn
	#print sign,da,dtn
	dd = xn-xt, yn-yt, zn-zt
	# print vecLength(dd)
	# '''
        # coor_strand.append([xt,yt,zt,sign])
        coor_strand.append([xn, yn, zn])
	# sign *= -1
    coor.append(coor_strand)

# print shifts


outfile = 'test_ca.pdb'
chain_ids = string.uppercase
# chain_ids = 'A'*26
fo = open(outfile, 'w')
resi = 1
atomi = 2
for i in range(len(coor)):
    coor_strand = coor[topology[i]-1]
    # if args.a and i%2 == 1: coor_strand.reverse()
    if args.a and i % 2 == 0:
        coor_strand.reverse()
    for j in range(len(coor_strand)):
        (x, y, z) = coor_strand[j]
        fo.write("ATOM  %5d  CA  ALA %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (atomi, chain_ids[i], resi, x, y, z))
        # fo.write("ATOM  %5d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           %1d\n" % (atomi,resi,x,y,z,sign))
        resi += 1
        atomi += 4
    # fo.write("TER   %5d      ALA %1s%4d\n" % (atomi,chain_ids[i],resi))
#    resi += 5
#    atomi += 36

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
