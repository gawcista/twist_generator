from twist_generator.tools import calc_angle,combine_lattice,calc_strain_hex

def build_bilayer(P1,P2,d,lattice,use_axis='A',translation=[0.,0.,0.],vacuum=None):
	latticeA=lattice.transform(P1)
	latticeB=lattice.transform(P=P2,p=translation)
	twist_angle=calc_angle(latticeA.basis[0],latticeB.basis[0])
	if use_axis=='A':
		latticeB.set_basis(latticeA.basis,ab_ony=True)
	elif use_axis=='B':
		latticeA.set_basis(latticeB.basis,ab_ony=True)
	print("The rotation angle is %.3f deg"%(twist_angle))
	bilayer=combine_lattice(latticeA,latticeB,d,vacuum=vacuum)
	return bilayer

def build_trilayer_ABC_hex(P1,P2,P3,d,lattA,lattB,lattC,autocorrection=True):
	#autocorrection: set the lattice constant as lattice B by default
	latticeA=lattA.transform(P1)
	latticeB=lattB.transform(P2)
	latticeC=lattC.transform(P3)
	print("The rotation angle is %.3f deg between layer 1&2"%(calc_angle(latticeA.basis[0],latticeB.basis[0])))
	if autocorrection and latticeA.cell_length()!=latticeB.cell_length():
		calc_strain_hex(latticeA.cell_length()[0],latticeB.cell_length()[0])
		latticeA.set_basis(latticeB.basis)
	bilayer=combine_lattice(latticeA,latticeB,d)
	print("The rotation angle is %.3f deg between layer 2&3"%(calc_angle(latticeA.basis[0],latticeC.basis[0])))
	if autocorrection and bilayer.cell_length()!=latticeC.cell_length():
		calc_strain_hex(bilayer.cell_length()[0],latticeC.cell_length()[0])
		latticeC.set_basis(bilayer.basis)
	trilayer=combine_lattice(bilayer,latticeC,d)
	return trilayer

def build_trilayer_ABC_hex(P1,P2,P3,d,lattA,lattB,lattC,autocorrection=True):
	#autocorrection: set the lattice constant as lattice B by default
	latticeA=lattA.transform(P1)
	latticeB=lattB.transform(P2)
	latticeC=lattC.transform(P3)
	print("The rotation angle is %.3f deg between layer 1&2"%(calc_angle(latticeA.basis[0],latticeB.basis[0])))
	if autocorrection and latticeA.cell_length()!=latticeB.cell_length():
		calc_strain_hex(latticeA.cell_length()[0],latticeB.cell_length()[0])
		latticeA.set_basis(latticeB.basis)
	bilayer=combine_lattice(latticeA,latticeB,d)
	print("The rotation angle is %.3f deg between layer 2&3"%(calc_angle(latticeA.basis[0],latticeC.basis[0])))
	if autocorrection and bilayer.cell_length()!=latticeC.cell_length():
		calc_strain_hex(bilayer.cell_length()[0],latticeC.cell_length()[0])
		latticeC.set_basis(bilayer.basis)
	trilayer=combine_lattice(bilayer,latticeC,d)
	return trilayer

def build_bilayer_aligned(P1,P2,d,lattA,lattB,use_axis='A',translation=[0.,0.,0.],vacuum=None):
	latticeA=lattA.transform(P1)
	latticeB=lattB.transform(P=P2,p=translation)
	if latticeA.cell_length()!=latticeB.cell_length():
		eta=calc_strain_hex(latticeA.cell_length()[0],latticeB.cell_length()[0],printflag=False)
		if use_axis=='A':
			latticeB.set_basis(latticeA.basis,ab_ony=True)
			print("%.4f%% strain was applied on layer A"%(100*eta))
		elif use_axis=='B':
			latticeA.set_basis(latticeB.basis,ab_ony=True)
			print("%.4f%% strain was applied on layer B"%(-100*eta))
	bilayer=combine_lattice(latticeA,latticeB,d,vacuum=vacuum)
	return bilayer
