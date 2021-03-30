import numpy as np
import networkx as nx
from rdkit import Chem

def MolFromGraphs(G):
	'''
	Function that takes as input the networkx graph (each node must have an atom property H,N,C,O etc) and return the mol (rdkit) object
	the function dont discriminate between different type of bond, it care only about the connectivity
	'''
	#https://stackoverflow.com/questions/51195392/smiles-from-graph
	# extract the adjacency matrix and the list of atoms
	adjacency_matrix = nx.adjacency_matrix(G).todense().tolist()
	node_list = []
	for i,node in enumerate(G):
		node_list.append(G.nodes[i]['atom'])
	# create empty editable mol object
	mol = Chem.RWMol()
	# add atoms to mol and keep track of index
	node_to_idx = {}
	for i in range(len(node_list)):
		a = Chem.Atom(node_list[i])
		molIdx = mol.AddAtom(a)
		node_to_idx[i] = molIdx
	# add bonds between adjacent atoms
	for ix, row in enumerate(adjacency_matrix):
		for iy, bond in enumerate(row):
			# only traverse half the matrix
			if iy >= ix:
				break
			# add relevant bond type (there are many more of these)
			if bond == 0:
				continue
			else:
				bond_type = Chem.rdchem.BondType.SINGLE
				mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
			#elif bond == 2:
			#    bond_type = Chem.rdchem.BondType.DOUBLE
			#    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)          
	Chem.SanitizeMol(mol)
	return mol

def chem_formula(form):
	charge = ''
	if form[-1] == '+' or form[-1] == '-':
		charge = form[-1]
		form = form[:-1]

	num = ['0','1','2','3','4','5','6','7','8','9']
	atoms = ['H','C','N','O']
	n_a = [0,0,0,0]
	name = ''
	for j in range(len(form)):
		if form[j] in atoms:
			if j == len(form) - 1:
				for k in range(len(atoms)):
					if form[j] == atoms[k]:
						n_a[k] = n_a[k] + 1
			else:       
				if form[j+1] in num:
					tmp_n = int(form[j+1])
					for k in range(len(atoms)):
						if form[j] == atoms[k]:
							n_a[k] = n_a[k] + tmp_n
				else:
					for k in range(len(atoms)):
						if form[j] == atoms[k]:
							n_a[k] = n_a[k] + 1
		elif form[j] in num:
			continue 					        
		else:
			name = 'ERROR'
			break

	if name != 'ERROR':		
		for i,ix in enumerate(atoms):
			if n_a[i] != 0:
				name = name + ix + str(n_a[i])
		name  = name + charge		
	return(name)