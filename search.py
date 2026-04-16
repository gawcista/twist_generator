from numpy import array,gcd
from twist_generator.tools import calc_angle
from twist_generator.generator import *

def search_square_100(Mcut,lattice,writeflag=False):
	for m in range(1,Mcut+1):
		for n in range(1,m):
			if gcd(m,n) == 1:
				P=array([[m,n],[m,-n]])
				vecs=P @ array([lattice.basis[0],lattice.basis[1]])
				angle=calc_angle(vecs[0],vecs[1])
				print("%d\t%d\t%.3f\t%d"%(m,n,angle,2*(m**2+n**2)*lattice.nions))
				if writeflag:
					generate_square_100(m,n,lattice,P1=P)

def search_square_210(Mcut,lattice,writeflag=False):
	for m in range(1,Mcut+1):
		for n in range(1,m):
			if gcd(m,n) == 1:
					P=array([[m,n],[n,m]])
					vecs=P @ array([lattice.basis[0],lattice.basis[1]])
					angle=calc_angle(vecs[0],vecs[1])
					print("%d\t%d\t%.3f\t%d"%(m,n,angle,2*(m**2+n**2)*lattice.nions))
					if writeflag:
						generate_square_210(m,n,lattice,P1=P)

def search_hexagonal_210(Mcut,lattice=None,writeflag=False):
	for m in range(1,Mcut+1):
		for n in range(-(m-1),2*m):
			if gcd(m,n) == 1 and m != 2*n:
				v1 = (array(lattice.basis[0])*m + array(lattice.basis[1])*n).tolist()
				v2 = (array(lattice.basis[0])*m - array(lattice.basis[1])*(n-m)).tolist()
				angle = calc_angle(v1,v2)
				print("%d\t%d\t%.3f\t%d"%(m,n,angle,(m**2+n**2-m*n)*lattice.nions))
				if writeflag:
					generate_hexagonal_210(m,n,lattice,filename='%.3f_%d.vasp'%(angle,(m**2+n**2-m*n)*3*lattice.nions))

def search_aligned(Mcut,latticeA,latticeB,strain_tol=0.001,writeflag=False):
	a = latticeA.cell_length()[0]
	b = latticeB.cell_length()[0]
	for m in range(1,Mcut+1):
		a_new = m*a
		n = round(a_new/b)
		b_new = n*b
		eta = (b_new-a_new)/a_new
		if abs(eta) < strain_tol:
			print("%d\t%d\t%.3f\t%.3f%%"%(m,n,a_new,eta*100))
			if writeflag:
				generate_aligned(latticeA,latticeB,m,n,d=4,filename="%d_%d.vasp"%(m,n))