from numpy import degrees,arccos,array,linalg
from copy import deepcopy

try:
	import spglib
except ImportError:
	spglib = None

def calc_angle(a,b): 
	# Calculate the angle between vectors <a,b>
	return degrees(arccos(array(a).dot(array(b))/(linalg.norm(a)*linalg.norm(b))))

def calc_mismatch(lattA,lattB,printflag=True):
	# consider strain was applied on lattB
	cellA = array(lattA.cell_length())
	cellB = array(lattB.cell_length())
	mismatch = (cellA-cellB)/cellB
	if printflag:
		if mismatch[0] - mismatch[1]<1e-8:
			print("%.4f%% biaxial strain should be applied!"%(100*mismatch[0]))
		else:
			print("%.4f%% strain should be applied on a!"%(100*mismatch[0]))
			print("%.4f%% strain should be applied on b!"%(100*mismatch[1]))
		if mismatch[2]!=0:
			print("%.4f%% mismatch was found on c!"%(100*mismatch[2]))
	return mismatch

def calc_strain(mp,nq,a,b,printflag=True):
	# old version to be updated
	a1 = (-nq*b**2/(mp))**0.5
	b1 = (-mp*a**2/(nq))**0.5
	delta_a = (a1-a)/a
	delta_b = (b1-b)/b
	if abs(delta_a)<abs(delta_b):
		if printflag:
			print("%.4f%% strain was applied on a"%(100*delta_a))
			return [a1/a,1,1]
		else:
			return 'a',a1/a-1
	else:
		if printflag:
			print("%.4f%% strain was applied on b"%(100*delta_b))
			return [1,b1/b,1]
		else:
			return 'b',b1/b-1
		
def calc_strain_hex(a1,a2,printflag=True):
	# returns the required strain on a1
	delta = (a1-a2)/a1
	if printflag:
		print("%.4f%% strain was applied on the first input layer"%(100*delta))
	return delta

def equals_vec(a,b,tol=1e-8):
	# determine whether two vectors are equivalent
	flag=True
	for i in range(len(a)):
		if abs(a[i]-b[i])>tol:
			flag=False
			break
	return flag

def spglib_cell(lattice):
	if spglib is None:
		raise ImportError("spglib is required for symmetry analysis.")
	if lattice.nions==0:
		raise ValueError("Symmetry analysis requires at least one atom.")
	scale = lattice.frac if getattr(lattice,'frac',1.0)>0 else 1.0
	cell = (array(lattice.basis)*scale).tolist()
	positions = []
	for atom in lattice.atom:
		position = _direct_position(lattice,atom)
		positions.append((position[0]%1.0,position[1]%1.0,position[2]%1.0))
	element_index = {}
	numbers = []
	for atom in lattice.atom:
		if atom.element not in element_index:
			element_index[atom.element] = len(element_index)+1
		numbers.append(element_index[atom.element])
	return (cell,positions,numbers)

def get_symmetry_info(lattice,symprec=1e-5,angle_tolerance=-1.0):
	cell = spglib_cell(lattice)
	dataset = spglib.get_symmetry_dataset(
		cell,
		symprec=symprec,
		angle_tolerance=angle_tolerance,
	)
	if dataset is None:
		return {
			'found': False,
			'symprec': symprec,
			'angle_tolerance': angle_tolerance,
		}
	return {
		'found': True,
		'number': int(dataset.number),
		'international': dataset.international,
		'hall_number': int(dataset.hall_number),
		'hall': dataset.hall,
		'pointgroup': dataset.pointgroup,
		'n_operations': len(dataset.rotations),
		'symprec': symprec,
		'angle_tolerance': angle_tolerance,
	}

def format_symmetry_info(info,label='Symmetry'):
	if not info['found']:
		return "%s: not found (symprec=%g, angle_tolerance=%g)"%(label,info['symprec'],info['angle_tolerance'])
	return (
		"%s: SG %d %s, Hall %d %s, point group %s, operations %d "
		"(symprec=%g, angle_tolerance=%g)"
	)%(
		label,
		info['number'],
		info['international'],
		info['hall_number'],
		info['hall'],
		info['pointgroup'],
		info['n_operations'],
		info['symprec'],
		info['angle_tolerance'],
	)

def print_symmetry_info(lattice,label='Symmetry',symprec=1e-5,angle_tolerance=-1.0):
	try:
		info = get_symmetry_info(lattice,symprec=symprec,angle_tolerance=angle_tolerance)
		print(format_symmetry_info(info,label=label))
		return info
	except ImportError as exc:
		print("%s: skipped (%s)"%(label,exc))
	except ValueError as exc:
		print("%s: skipped (%s)"%(label,exc))
	return None

def _check_z_aligned_cell(lattice):
	# The z placement below assumes a slab cell whose c axis is normal to xy.
	tol = 1e-8
	basis = lattice.basis
	if abs(basis[0][2])>tol or abs(basis[1][2])>tol or abs(basis[2][0])>tol or abs(basis[2][1])>tol:
		raise ValueError("combine_lattice only supports slab cells with c axis along z.")
	if basis[2][2]<=0:
		raise ValueError("combine_lattice requires a positive c lattice vector.")
	return basis[2][2]

def _direct_position(lattice,atom):
	if lattice.tag=='Cartesian':
		return (linalg.inv(array(lattice.basis).T).dot(array(atom.position))).tolist()
	return atom.position[:]

def _slab_z_profile(lattice):
	# Identify the external vacuum as the largest periodic gap in fractional z.
	c = _check_z_aligned_cell(lattice)
	if lattice.nions==0:
		raise ValueError("combine_lattice requires at least one atom in each lattice.")

	positions = [_direct_position(lattice,atom) for atom in lattice.atom]
	z_values = sorted([pos[2]%1.0 for pos in positions])
	gaps = []
	for i,z in enumerate(z_values):
		z_next = z_values[(i+1)%len(z_values)]
		if i==len(z_values)-1:
			gap = z_values[0]+1.0-z
		else:
			gap = z_next-z
		gaps.append((gap,z,z_next))

	gaps_sorted = sorted(gaps,key=lambda item:item[0],reverse=True)
	max_gap,gap_start,gap_end = gaps_sorted[0]
	if max_gap*c<=1e-8:
		raise ValueError("Could not identify external vacuum along c axis.")
	if len(gaps_sorted)>1 and abs(max_gap-gaps_sorted[1][0])*c<=1e-8:
		raise ValueError("Could not identify a unique external vacuum gap along c axis.")

	thickness_frac = 1.0-max_gap
	slab_bottom = gap_end%1.0
	rel_z = []
	for pos in positions:
		z_rel_frac = (pos[2]%1.0-slab_bottom)%1.0
		if abs(z_rel_frac-thickness_frac)*c<=1e-8:
			z_rel_frac = thickness_frac
		if z_rel_frac-thickness_frac>1e-8:
			raise ValueError("Could not place atom outside identified slab region.")
		rel_z.append(z_rel_frac*c)

	return {
		'positions': positions,
		'thickness': thickness_frac*c,
		'vacuum': max_gap*c,
		'rel_z': rel_z,
	}

def combine_lattice(latticeA,latticeB,d,vacuum=None):
	# combine two lattices and set d as the Cartesian distance between inner slab surfaces
	if d<0:
		raise ValueError("Layer distance d must be non-negative.")
	if vacuum is not None and vacuum<0:
		raise ValueError("Vacuum must be non-negative.")
	_check_z_aligned_cell(latticeA)
	_check_z_aligned_cell(latticeB)
	if equals_vec(latticeA.cell_angle(),latticeB.cell_angle()) and equals_vec(latticeA.cell_length()[:2],latticeB.cell_length()[:2]):
		profileA = _slab_z_profile(latticeA)
		profileB = _slab_z_profile(latticeB)
		vacuum_final = max(profileA['vacuum'],profileB['vacuum']) if vacuum is None else vacuum
		c_final = profileA['thickness']+d+profileB['thickness']+vacuum_final
		if c_final<=0:
			raise ValueError("Final c lattice length must be positive.")

		latticeNew = deepcopy(latticeA) # define the new lattice
		new_basis = deepcopy(latticeA.basis)
		new_basis[2] = [0.,0.,c_final]
		latticeNew.set_basis(new_basis)
		latticeNew.tag = 'Direct'

		layerA_bottom = vacuum_final/2.0
		layerB_bottom = layerA_bottom+profileA['thickness']+d
		for atom,pos,z_rel in zip(latticeNew.atom,profileA['positions'],profileA['rel_z']):
			atom.set_position([pos[0],pos[1],(layerA_bottom+z_rel)/c_final])
		for atom,pos,z_rel in zip(latticeB.atom,profileB['positions'],profileB['rel_z']):
			atom_new = deepcopy(atom)
			atom_new.set_position([pos[0],pos[1],(layerB_bottom+z_rel)/c_final])
			latticeNew.add_atom(atom_new)
		latticeNew.thickness = (profileA['thickness']+d+profileB['thickness'])/c_final
		return latticeNew
	else:
		print('WARNING: cell angles or sized does not match!! The CONTCAR may be incorrect!!' )

def is_hexagonal(lattice,tolerance=1e-5):
	# Determine whether a lattice is hexagonal
    a = lattice.cell_length()[0]
    b = lattice.cell_length()[1]
    alpha = lattice.cell_angle()[2]
    if abs(a - b) < tolerance and abs(alpha - 120) < tolerance:
        return True
    else:
        return False
	
def basis_rotation(lattice,type='hex'):
    if type == 'hex':
        a = lattice.cell_length()[0]
        c = lattice.cell_length()[2]
        newbasis = [[a,0,0],
                    [-a/2,a*(3**0.5)/2,0],
                    [0,0,c]]
        lattice.set_basis(newbasis)
        return lattice
