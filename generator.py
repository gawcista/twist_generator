from twist_generator.build import build_bilayer,build_bilayer_aligned

def generate_bilayer(P1,P2,lattice,d=4,filename='CONTCAR',vacuum=None,symmetry=True,symprec=1e-5,angle_tolerance=-1.0,**kargs):
    bilayer=build_bilayer(P1,P2,d,lattice,vacuum=vacuum,**kargs)
    bilayer.print_POSCAR(filename,symmetry=symmetry,symprec=symprec,angle_tolerance=angle_tolerance)

def generate_square_100(m,n,**kargs):
    P1=[[m,-n,0],[n,m,0],[0,0,1]]
    P2=[[m,n,0],[-n,m,0],[0,0,1]]
    generate_bilayer(P1,P2,**kargs)

def generate_square_110(m,n,**kargs):
    P1=[[m,-n,0],[n,m,0],[0,0,1]]
    P2=[[n,-m,0],[m,n,0],[0,0,1]]
    generate_bilayer(P1,P2,**kargs)

def generate_hexagonal_210(m,n,**kargs):
    P1=[[m,-n,0],[n,m-n,0],[0,0,1]]
    P2=[[m,n-m,0],[m-n,n,0],[0,0,1]]
    generate_bilayer(P1,P2,**kargs)

def generate_aligned(latticeA,latticeB,m,n,d=4,filename='CONTCAR',vacuum=None,symmetry=True,symprec=1e-5,angle_tolerance=-1.0,**kargs):
    P1=[[m,0,0],[0,m,0],[0,0,1]]
    P2=[[n,0,0],[0,n,0],[0,0,1]]
    bilayer=build_bilayer_aligned(P1,P2,d,latticeA,latticeB,vacuum=vacuum,**kargs)
    bilayer.print_POSCAR(filename,symmetry=symmetry,symprec=symprec,angle_tolerance=angle_tolerance)
