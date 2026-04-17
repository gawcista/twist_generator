from twist_generator.tools import calc_angle,get_symmetry_info,print_symmetry_info
from numpy import linalg,floor,array
from copy import deepcopy

class atom_unit:
	element=''
	dynamics=''
	position = []

	def __init__(self,position=[],element='',dynamics=['T','T','T']):
		self.position = position
		self.element = element
		self.dynamics = dynamics

	def read_poslines(self,pos_line):
		# read data from a position line (string): 
		# "r_a r_b r_c d_a d_b d_c" 
		# r_i is the fractional position on i-axis; 
		# d_i is a boolean value decides the states of selective dynamics
		pos_list = pos_line.split()
		self.position = [float(pos_list[0]),float(pos_list[1]),float(pos_list[2])]
		if len(pos_list) >= 6:
			self.dynamics = [pos_list[3],pos_list[4],pos_list[5]]

	def set_position(self,position):
		self.position = position

	def set_position_c(self,c):
		self.position[2] = c

	def set_element(self,element):
		self.element = element

	def set_dynamics(self,dynamics):
		self.dynamics = dynamics

	def to_Cartesian(self,basis):
		# convert fractional position to Cartesian position
		cartesian = (array(basis).T.dot(array(self.position))).tolist()
		return cartesian

class atomic_structure:
	title = ''
	frac = 1.0
	basis = []
	elements = []
	natoms = []
	selective_dynamics = 1
	tag = 'Direct'
	atom = []
	nions = 0
	thickness = 0

	def __init__(self,file='POSCAR',mode='vasp',title=None,symmetry=True,symprec=1e-5,angle_tolerance=-1.0,*args, **kwargs):
		if mode == 'vasp':
			self.generate_from_vasp(file,*args, **kwargs)
		if mode == 'lammps':
			self.generate_from_lammps(file,*args, **kwargs)
		if mode == 'empty':
			self.frac = 1.0
			self.basis = []
			self.elements = []
			self.natoms = []
			self.selective_dynamics = 1
			self.tag = 'Direct'
			self.atom = []
			self.nions = 0
			self.thickness = 0

		if title is not None:
			self.title = title
		if mode!='empty' and symmetry:
			self.print_symmetry(label='Input symmetry',symprec=symprec,angle_tolerance=angle_tolerance)

	def generate_from_vasp(self,file='POSCAR',auto_standardize=False):
		# initialize a lattice from a POSCAR file
		with open(file,'r') as fin:
			POSCAR_lines = fin.readlines()
		self.title = POSCAR_lines[0][:-2]
		self.frac = float(POSCAR_lines[1].split()[0])
		self.basis = []
		for i in range(3):
			x = float(POSCAR_lines[i+2].split()[0])
			y = float(POSCAR_lines[i+2].split()[1])
			z = float(POSCAR_lines[i+2].split()[2])
			self.basis.append([x,y,z])
		self.elements = POSCAR_lines[5].split()
		self.natoms = []
		for s in POSCAR_lines[6].split():
			self.natoms.append(int(s))
		if POSCAR_lines[7].split()[0][0] in ['S','s']:
			self.selective_dynamics = 1
		else:
			self.selective_dynamics = 0
		if POSCAR_lines[7+self.selective_dynamics].split()[0][0] in ['D','d']:
			self.tag = 'Direct'
		elif POSCAR_lines[7+self.selective_dynamics].split()[0][0] in ['C','c']:
			self.tag = 'Cartesian'
		self.nions = 0
		self.atom = []
		for i in range(len(self.elements)):
			for j in range(self.natoms[i]):
				self.nions += 1
				atom_ij = atom_unit(element=self.elements[i])
				atom_ij.read_poslines(POSCAR_lines[8+self.selective_dynamics+self.nions-1])
				self.atom.append(atom_ij)
		if auto_standardize:
			self.standardize()

	def generate_from_lammps(self,file='POSCAR',mass_dict=None,auto_standardize=False):
		if mass_dict is None:
			mass_dict = {10.811:'B',14.0067:'N',12.0107:'C'}
		with open(file,'r') as fin:
			lines_raw = fin.readlines()
			lines = [line.strip() for line in lines_raw if line.strip()]
		self.title = lines[0]
		self.frac = 1.0
		self.nions = int(lines[1].split()[0])
		self.tag = 'Cartesian'
		ntypes = int(lines[2].split()[0])
		xx = float(lines[3].split()[1])-float(lines[3].split()[0])
		yy = float(lines[4].split()[1])-float(lines[4].split()[0])
		zz = float(lines[5].split()[1])-float(lines[5].split()[0])
		xy = float(lines[6].split()[0])
		xz = float(lines[6].split()[1])
		yz = float(lines[6].split()[2])
		self.set_basis([[xx,0.,0.],[xy,yy,0.],[xz,yz,zz]])
		element_dict={}
		atom_element=[]
		for i in range(ntypes):
			index_element = int(lines[8+i].split()[0])
			mass_element = float(lines[8+i].split()[1])
			element_dict[index_element] = mass_dict[mass_element]
			if mass_dict[mass_element] not in self.elements:
				self.elements.append(mass_dict[mass_element])
				self.natoms.append(0)
				atom_element.append([])
		line_start = 9 + ntypes
		for i in range(self.nions):
			element  = element_dict[int(lines[line_start+i].split()[2])]
			x = float(lines[line_start+i].split()[4])
			y = float(lines[line_start+i].split()[5])
			z = float(lines[line_start+i].split()[6])
			atom_i = atom_unit([x,y,z],element)
			self.natoms[self.elements.index(element)] += 1
			
			atom_element[self.elements.index(element)].append(atom_i)
		self.atom = sum(atom_element,[])
		if auto_standardize:
			self.standardize()

	def reset(self):
		self.frac = 1.0
		self.basis = []
		self.elements = []
		self.natoms = []
		self.selective_dynamics = 1
		self.tag = 'Direct'
		self.atom = []
		self.nions = 0

	def set_basis(self,basis,ab_ony=False):
		if ab_ony:
			self.basis[0]= basis[0]
			self.basis[1]= basis[1]
		else:
			self.basis=basis

	def atom_in_cell(self,position):
		# Decide whether a position is in the cell
		# trans: if not in the cell, which neighbor cell it's in
		tol=1e-4
		flag = True
		trans=[0,0,0]
		for i in range(3):
			if abs(position[i])<tol:
				position[i]=0.
			if abs(position[i]-1)<tol:
				position[i]=1.0	
			if not(position[i]>=0 and position[i]<1.):
				flag = False
				trans[i] = int(floor(position[i]))
		return flag,trans

	def add_atom(self,atom):
		if atom.element not in self.elements:
			self.elements.append(atom.element)
			self.natoms.append(0)
		i = self.elements.index(atom.element)
		pos = 0
		for j in range(i+1):
			pos += self.natoms[j]
		self.nions += 1
		self.natoms[i] += 1
		self.atom.insert(pos,atom)

	def del_atom(self,i):
		self.nions -= 1
		index = self.elements.index(self.atom[i].element)
		self.natoms[index] -= 1
		del(self.atom[i])

	def cell_angle(self):
		alpha = round(calc_angle(self.basis[0],self.basis[2]),6)
		beta = round(calc_angle(self.basis[1],self.basis[2]),6)
		gamma = round(calc_angle(self.basis[0],self.basis[1]),6)
		return [alpha,beta,gamma]

	def cell_length(self):
		a = linalg.norm(self.basis[0])
		b = linalg.norm(self.basis[1])
		c = linalg.norm(self.basis[2])
		return [a,b,c]

	def get_symmetry(self,symprec=1e-5,angle_tolerance=-1.0):
		return get_symmetry_info(self,symprec=symprec,angle_tolerance=angle_tolerance)

	def print_symmetry(self,label='Symmetry',symprec=1e-5,angle_tolerance=-1.0):
		return print_symmetry_info(self,label=label,symprec=symprec,angle_tolerance=angle_tolerance)

	def strain(self,strain_list):
		for i in range(3):
			self.basis[i]=(array(self.basis[i])*strain_list[i]).tolist()
	
	def add_vaccum(self,d=None,type='add'):
		old_c = self.basis[-1][-1]
		if type == 'add':
			new_c = old_c + d
		elif type == 'set':
			new_c = d
		if self.tag=='Direct':
			for i in range(self.nions):
				self.atom[i].set_position_c(self.atom[i].position[-1]*old_c/new_c)
		self.basis[-1][-1]=new_c
		self.standardize()

	def cart2direct(self):
		if self.tag !='Direct':
			self.tag = 'Direct'
			basis_inv = linalg.inv(array(self.basis).T)
			for atom in self.atom:
				atom.position=basis_inv @ atom.position
	
	def direct2cart(self):
		if self.tag !='Cartesian':
			self.tag = 'Cartesian'
			for atom in self.atom:
				atom.position=array(self.basis).T @ array(atom.position)


	def standardize(self):
		def vec_shift(v):
			v_new = v[:]
			for i in range(3):
				flag=True
				while flag:
					if v_new[i]<-0.5:
						v_new[i] += 1
						flag=True
					elif v_new[i]>=0.5:
						v_new[i] -= 1
						flag=True
					else:
						flag=False
			return v_new
		if self.tag == 'Cartesian':
			self.cart2direct()
		c_min = 1
		c_max = -1
		for atom in self.atom:
			atom.position = vec_shift(atom.position)   

			if atom.position[-1]>c_max:
				c_max=atom.position[2]
			if atom.position[-1]<c_min:
				c_min=atom.position[2]
		if c_max-c_min>0.5:
			self.thickness = 1-(c_max-c_min) 
		else:
			self.thickness = c_max-c_min
		center = c_min+self.thickness/2
		for atom in self.atom:
			atom.set_position_c((atom.position[2]+0.5-center)%1)

	def checkcenter(self):
		if self.tag == 'Cartesian':
			self.cart2direct()
		c_min = 1
		c_max = -1
		for atom in self.atom:   
			if atom.position[-1]>c_max:
				c_max=deepcopy(atom.position[2])
			if atom.position[-1]<c_min:
				c_min=deepcopy(atom.position[2])
		if c_max-c_min>0.5:
			thickness = 1-(c_max-c_min) 
		else:
			thickness = c_max-c_min
		center = c_min+thickness/2
		print(c_min,c_max,center,thickness)

	def transform(self,P=[[1,0,0],[0,1,0],[0,0,1]],p=[0,0,0]):
		latticeNew=deepcopy(self)
		latticeNew.reset()
		if linalg.det(array(P)) == 0:
			print("-Error: The determination of the transform matrix is 0!")
		else:
			if linalg.det(array(P)) < 0:
				print("-Warning: The transform matrix changes the coordinate system from right- to left-handed.")
			if linalg.det(array(P)) != 1:
				print("-The transform matrix changes the cell volume.")
			latticeNew.set_basis((array(self.basis).T.dot(array(P))).T)
			inv_P = linalg.pinv(array(P)) #the atom position vector should dot this
			#################search atoms in the new basis#####################
			stack = [[0,0,0],[1,0,0],[0,1,0],[0,-1,0],[-1,0,0],[0,0,1],[0,0,-1]]
			search_map = []
			while len(stack)!=0:
				n_new  = 0
				for i in range(self.nions):
					newatom = deepcopy(self.atom[i])
					position = ((array(self.atom[i].position)+array(stack[0])).dot(array(inv_P).T)+array(p)).tolist()
					flag,trans = self.atom_in_cell(position)
					if flag:
						newatom.set_position(position)
						latticeNew.add_atom(newatom)
						n_new += 1
				if n_new>0:
					for nextcell in [[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,1],[0,0,-1]]:
						cell = array(nextcell)+array(stack[0])
						if cell.tolist() not in search_map and cell.tolist() not in stack:
							stack.append(cell.tolist())
				search_map.append(stack.pop(0))
		latticeNew.standardize()
		return latticeNew

	def print_POSCAR(self,file='CONTCAR',symmetry=True,symprec=1e-5,angle_tolerance=-1.0):
		with open(file,'wt') as fout:
			print(self.title,file=fout)
			print("%19.14f"%(self.frac),file=fout)
			for i in range(3):
				print(" %22.16f %22.16f %22.16f"%(self.basis[i][0],self.basis[i][1],self.basis[i][2]),file=fout)
			for val in self.elements:
				print("   %-2s"%(val),end='',file=fout)
			print(file=fout)
			for val in self.natoms:
				print("   %-2d"%(val),end='',file=fout)  #vasp POSCAR standard "%6d"
			print(file=fout)
			if self.selective_dynamics == 1:
				print('Selective Dynamics',file=fout)
			print(self.tag,file=fout)
			for i in range(self.nions):
				print(" %20.16f %20.16f %20.16f"%(self.atom[i].position[0],self.atom[i].position[1],self.atom[i].position[2]),end='',file=fout)
				if self.selective_dynamics == 1:
					print(" %s %s %s"%(self.atom[i].dynamics[0],self.atom[i].dynamics[1],self.atom[i].dynamics[2]),file=fout)
				else:
					print(file=fout)
		if symmetry:
			self.print_symmetry(label='Output symmetry',symprec=symprec,angle_tolerance=angle_tolerance)
