from itertools import combinations, chain, permutations
import networkx as nx
import matplotlib.pyplot as plt
import random 

random.seed(0)

class Atom:
	def __init__(self,info):
		self.atom           = info[0]
		self.electrons      = info[1]
		self.atomic_numb    = info[2]
		self.valence        = info[3]
		self.orbitals		= info[4]

def powerset(someset):
	"""
	Function that gives it back all the possible sub-set of a set (powerset)
	"""
	_empty_powerset = ((), )
	try:
		someset.isdisjoint
	except AttributeError:
		raise TypeError(
			f"{powerset.__name__} accepts only a set-like object as parameter"
		) from None
	size = len(someset)
	combs = (combinations(someset, k) for k in range(1, size+1))
	return chain(_empty_powerset, *combs)

def create_list_graph(num_nodes):
	'''
	Function that gives it back all possible graph, in networkx object,  with n-nodes
	'''
	# Compute all possible edge like possible combination between two numbers
	edges = set(combinations(range(num_nodes), 2))
	# Check that all the possibilities must be equal to 0.5*n*(n-1)
	assert len(edges) == num_nodes*(num_nodes - 1)/2
	# Compute all the possible graphs, from the empty graph to the completed
	possible_graphs = powerset(edges)
	list_of_graphs = []
	for graph in powerset(edges):
		G = nx.Graph()
		G.add_nodes_from(range(num_nodes))
		for edge in list(graph):
			G.add_edge(*edge)
		list_of_graphs.append(G)     
	return list_of_graphs	

def mol_graph_image(G):
	'''
	Draw graph like molecule
	'''
	# extract nodes with specific setting of the attribute
	oxigen_nodes = [n for (n,at) in \
		nx.get_node_attributes(G,'atom').items() if at == 'O']
	carbon_nodes = [n for (n,at) in \
		nx.get_node_attributes(G,'atom').items() if at == 'C']
	nitrogen_nodes = [n for (n,at) in \
		nx.get_node_attributes(G,'atom').items() if at == 'N']
	hydrogen_nodes = [n for (n,at) in \
		nx.get_node_attributes(G,'atom').items() if at == 'H']	
	# and find all the remaining nodes.
	other_nodes = list(set(G.nodes()) - set(oxigen_nodes) - set(carbon_nodes) - set(nitrogen_nodes) - set(hydrogen_nodes))

	# now draw them in subsets  using the `nodelist` arg
	pos = nx.spring_layout(G)
	nx.draw_networkx_nodes(G, pos, nodelist=oxigen_nodes, \
		node_color='red', node_shape='o')
	nx.draw_networkx_nodes(G, pos, nodelist=carbon_nodes, \
		node_color='black', node_shape='o')
	nx.draw_networkx_nodes(G, pos, nodelist=nitrogen_nodes, \
		node_color='blue', node_shape='o')
	nx.draw_networkx_nodes(G, pos, nodelist=hydrogen_nodes, \
		node_color='grey', node_shape='o')	
	nx.draw_networkx_nodes(G, pos, nodelist=other_nodes, \
		node_color='purple', node_shape='o')
	
	edges = G.edges()
	try: 
		weights = [G[u][v]['bond'] for u,v in edges]
		for i,item in enumerate(weights):
			if item == 2:
				weights[i] = 'dashed'
			elif item == 3:
				weights[i] = 'dotted'
			else:
				weights[i] = 'solid'	
		nx.draw_networkx_edges(G, pos, width=2.0, style=weights, alpha=0.5) 
	except:
		nx.draw_networkx_edges(G, pos, width=2.0, style='solid', alpha=0.5) 
	
	plt.savefig("graph.png", format="PNG")
	plt.show()
	plt.clf()
	return

def atom_valence(G):
	'''
	Function that form a graph gives it back if all the atom respect the valence, return(true/false, list of index-node)
	'''
	tmp_bol = True
	n_atom = []
	for k,atom in enumerate(G):
		if G.degree[k] > G.nodes[k]['valence']:
			tmp_bol = False
			n_atom.append(k)
	return(tmp_bol,n_atom)

def hierarchy_pos(G, root=None, width=20., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5):

    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.  
    Licensed under Creative Commons Attribution-Share Alike 

    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.

    G: the graph (must be a tree)

    root: the root node of current branch 
    - if the tree is directed and this is not given, 
      the root will be found and used
    - if the tree is directed and this is given, then 
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given, 
      then a random choice will be used.

    width: horizontal space allocated for this branch - avoids overlap with other branches

    vert_gap: gap between levels of hierarchy

    vert_loc: vertical location of root

    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=20., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos


    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)

def del_iso_graph(tmp,prop):
	'''
	delete graphs in a list (tmp) that are isomorph respect to an graph attribute (prop) or a list of attributes
	'''
	if prop == None:
		for i, graph_1 in enumerate(tmp):
			for j, graph_2 in enumerate(tmp[i+1:]):
				if nx.is_isomorphic(graph_1,graph_2):
					tmp.remove(graph_2)
	else:
		nm = nx.algorithms.isomorphism.categorical_node_match(prop,prop)
		for i, graph_1 in enumerate(tmp):
			for j, graph_2 in enumerate(tmp[i+1:]):
				if nx.is_isomorphic(graph_1,graph_2,node_match=nm):
					tmp.remove(graph_2)
	return(tmp)

def del_iso_graph_bond(tmp):
	'''
	delete graphs in a list (tmp) that are isomorph respect to an graph attrobute (prop) 
	'''
	prop = 'bond'
	nm = nx.algorithms.isomorphism.categorical_edge_match(prop,prop)
	for i, graph_1 in enumerate(tmp):
		for j, graph_2 in enumerate(tmp[i+1:]):
			if nx.is_isomorphic(graph_1,graph_2,edge_match=nm):
				tmp.remove(graph_2)
	return(tmp)

def del_iso_mol(tmp):
	'''
	delete molecules in a list (tmp) that are isomorph respect to bonds and electronic structures 
	'''
	prop = ['atom','non_b_electrons']
	nm_a = nx.algorithms.isomorphism.categorical_node_match(prop,prop)
	nm_b = nx.algorithms.isomorphism.categorical_edge_match('bond','bond')
	for i, graph_1 in enumerate(tmp):
		for j, graph_2 in enumerate(tmp[i+1:]):
			if nx.is_isomorphic(graph_1,graph_2,edge_match=nm_b,node_match=nm_a):
				tmp.remove(graph_2)
	return(tmp)

def H_addiction(G,dict_element):
	'''
	Function that add in all possible position an hydrogen to the graph, create a list of allowed molecules by valence and isomerism
	'''
	H = nx.Graph()
	H.add_node(len(G))
	list_gr = []
	H.nodes[len(G)].update(dict_element[0])
	for i,node in enumerate(G):
		tmp_gr = H.copy()
		tmp_gr.add_edge(node,len(G))
		tmp_gr = nx.compose(G,tmp_gr).copy()
		list_gr.append(tmp_gr)
	tmp_l = del_iso_graph(list_gr,'atom').copy()
	list_gr = []
	for i,item in enumerate(tmp_l):
		if atom_valence(item)[0]:
			list_gr.append(item)
	return(list_gr)

def total_bonds(G,node):
	'''
	Function that take as input a moleculer-graph and the label of a specific node and then it gives out the number of bonds of that specific node/atom
	'''
	tot_bonds = 0
	for (u, v) in G.edges(node):
		tot_bonds = tot_bonds + G[u][v]['bond']
	return(tot_bonds)

def mol_total_bonds(G):
	'''
	Function that take as input a moleculer-graph and then it gives out the list of the number of bonds of that molecules (single,double,triple)
	'''
	bond = [0,0,0]

	for (u, v) in G.edges:
		if G[u][v]['bond'] == 1:
			bond[0] = bond[0] + 1
		elif G[u][v]['bond'] == 2:
			bond[1] = bond[1] + 1
		else:
			bond[2] = bond[2] + 1
	return(bond)

def possible_bonds(G):
	'''
	Function that takes as input a moleculer-graph and give you back a list of all possibile, if is possible, edges/sites, type tuple, available for double or triple bonds
	'''
	tmp_edge = []
	for (u, v) in G.edges:
		if total_bonds(G,u) < G.nodes[u]['valence'] and total_bonds(G,v) < G.nodes[v]['valence'] and G[u][v]['bond'] < 3:
			tmp_edge.append([u,v])
	return(tmp_edge)
	
def all_bonds_molecules(G):
	'''
	Function that takes as input a moleculer-graph and give you back a list of all possibile, if is possible, molecules with double or triple bonds
	'''
	tmp_edge = possible_bonds(G)
	tmp_mol = []
	for i,edge in enumerate(tmp_edge):
		tmp_g = G.copy()
		tmp_g[edge[0]][edge[1]]['bond'] = tmp_g[edge[0]][edge[1]]['bond'] + 1
		tmp_mol.append(tmp_g.copy())
	return(tmp_mol)

def possible_electron(list_electron,av_electron):
	'''
	Function that takes as input a list, in which each element rapresent an atom avialable to have at least one electron and the element list[i] rapresent the maximum electrion avialable for that atom, 
	and the total number of electrons avialable in the molecule to put as non-bond electron 
	'''
	rng = []
	for i,item in enumerate(list_electron):
		rng = rng + list(range(item + 1))
	#rng = list(range(av_electron + 1)) * list_electron
	tmp = set(i for i in permutations(rng, len(list_electron)) if sum(i) == av_electron)
	tmp_l = []
	for i,item in enumerate(tmp):
		a = True
		for k,elem in enumerate(item):
			if elem > list_electron[k]:
				a = False
				break
		if a == True:
			tmp_l.append(item)		
	return(tmp_l)

def list_molecules_bonds(G):
	'''
	Function that takes in input a graph and gives it back all the possible allowed bond-isomer molecules
	'''
	#search and store all avialable molecules with different bonds
	list_mol = [G.copy()]
	tmp_i = 0
	while True:
		tmp_list = all_bonds_molecules(list_mol[tmp_i])
		list_mol = list_mol + tmp_list
		tmp_i = tmp_i + 1
		if len(list_mol) == tmp_i:
			break
	#delete isomorph graphs by bonds type 
	list_mol = del_iso_graph_bond(list_mol)
	#update number electron in each atoms and delete molecules that have more electron in bonds then total number electron
	tmp_l = []  
	for k,mol in enumerate(list_mol):
		tmp_e = 0
		for i,node in enumerate(mol):
			mol.nodes[i]['bonds'] = total_bonds(mol,i)
			mol.nodes[i]['electrons'] = 2 * mol.nodes[i]['bonds']
			tmp_e = tmp_e + mol.nodes[i]['electrons'] // 2
		if tmp_e <= G.graph['tot_electrons']:
			mol.graph['electrons'] = tmp_e
			tmp_b = mol_total_bonds(mol)
			mol.graph['N_bond_single'] = tmp_b[0]
			mol.graph['N_bond_double'] = tmp_b[1]
			mol.graph['N_bond_triple'] = tmp_b[2]
			tmp_l.append(mol.copy())
	list_mol = tmp_l.copy()
	return(list_mol)

def list_resonance_structure(G):
	'''
	Function that takes in input a graph and gives it back all the possible allowed way to put non-bonds electrons in the molecules, as a list.
	The number of lone pairs and radicals on each element of the list respect the "pauli principle" of max molteplity of spin
	'''	
	av_electron = G.graph['tot_electrons'] - G.graph['electrons'] #number of total electrons avialble to put in the molecule
	list_ris = []
	if av_electron != 0:
		tmp_l = [] # atom-label of avialable atoms to have more electron
		tmp_e = [] # max electron avialable for the atoms index by label
		for i,node in enumerate(G):
			tmp_el = (G.nodes[i]['orbitals']-total_bonds(G,i))*2
			if tmp_el != 0:
				tmp_l.append(i)
				tmp_e.append(tmp_el)
		tmp_c = possible_electron(tmp_e,av_electron)
		for k,item in enumerate(tmp_c):
			tmp_g = G.copy()
			for j in range(len(item)):
				tmp_g.nodes[tmp_l[j]]['non_b_electrons'] = item[j]
				tmp_g.nodes[tmp_l[j]]['electrons'] = tmp_g.nodes[tmp_l[j]]['electrons'] + item[j]
			#put the electron in max molteplicity way (pauli principle)
			for i,node in enumerate(tmp_g):
				tmp_nb = tmp_g.nodes[i]['non_b_electrons']
				tmp_or  = tmp_g.nodes[i]['orbitals']-total_bonds(tmp_g,i)
				if tmp_nb != 0 and tmp_or != 0:
					if (tmp_nb / tmp_or) <= 1:
						tmp_g.nodes[i]['radicals'] = tmp_nb
						tmp_g.nodes[i]['lone_pair'] = 0
					else:
						tmp_g.nodes[i]['radicals'] = tmp_or * 2 - tmp_nb
						tmp_g.nodes[i]['lone_pair'] = tmp_or - tmp_g.nodes[i]['radicals']
				tmp_g.nodes[i]['formal_charge'] = tmp_g.nodes[i]['valence_electron'] - tmp_g.nodes[i]['bonds'] - tmp_g.nodes[i]['non_b_electrons']
				tmp_g.graph['radicals'] = tmp_g.graph['radicals'] + tmp_g.nodes[i]['radicals']		
				tmp_g.graph['lone_pairs'] = tmp_g.graph['lone_pairs'] + tmp_g.nodes[i]['lone_pair']
				tmp_g.graph['abs_total_formal_charge'] = tmp_g.graph['abs_total_formal_charge'] + abs(tmp_g.nodes[i]['formal_charge'])
			list_ris.append(tmp_g.copy())
		list_ris = del_iso_mol(list_ris)	
	else:
		list_ris.append(G)
	return(list_ris)

def isomers_generator(n_h,n_c,n_o,n_n):
	'''
	Function that take as input the number of: H,C,O,N and gives it back the tree of the generation of all the possibile molecules. This last one are the leafs of the tree. 
	'''
	H = Atom(("H",1,1,1,1))
	C = Atom(("C",4,6,4,4))
	N = Atom(("N",5,7,4,4))
	O = Atom(("O",6,8,2,4))
	tmp_h = {'atom': H.atom, 'valence_electron': H.electrons, 'atomic_numb': H.atomic_numb, 'valence': H.valence, 'orbitals': H.orbitals, 'lone_pair': 0,'bonds': 0, 'electrons': 0, 'radicals': 0, 'formal_charge': 0, 'non_b_electrons': 0}
	tmp_c = {'atom': C.atom, 'valence_electron': C.electrons, 'atomic_numb': C.atomic_numb, 'valence': C.valence, 'orbitals': C.orbitals, 'lone_pair': 0,'bonds': 0, 'electrons': 0, 'radicals': 0, 'formal_charge': 0, 'non_b_electrons': 0}
	tmp_n = {'atom': N.atom, 'valence_electron': N.electrons, 'atomic_numb': N.atomic_numb, 'valence': N.valence, 'orbitals': N.orbitals, 'lone_pair': 0,'bonds': 0, 'electrons': 0, 'radicals': 0, 'formal_charge': 0, 'non_b_electrons': 0}
	tmp_o = {'atom': O.atom, 'valence_electron': O.electrons, 'atomic_numb': O.atomic_numb, 'valence': O.valence, 'orbitals': O.orbitals, 'lone_pair': 0,'bonds': 0, 'electrons': 0, 'radicals': 0, 'formal_charge': 0, 'non_b_electrons': 0}
	dict_element = [tmp_h.copy(),tmp_c.copy(),tmp_n.copy(),tmp_o.copy()]
	del tmp_h, tmp_c, tmp_n, tmp_o

	NA     = n_h + n_c + n_o + n_n
	NA_con = NA - n_h
	list_atom = []
	list_con_atom = []
	tot_el = n_h * H.electrons + n_c * C.electrons + n_o * O.electrons + n_n * N.electrons
	#list of connectivity atom
	for i in range(n_h):
		list_atom.append(H.atom)
	for i in range(n_c):
		list_atom.append(C.atom)
		list_con_atom.append(C.atom)
	for i in range(n_o):
		list_atom.append(O.atom)
		list_con_atom.append(O.atom)
	for i in range(n_n):
		list_atom.append(N.atom)
		list_con_atom.append(N.atom)

	#generate tree
	tree = nx.Graph()
	root = nx.Graph()
	root.add_nodes_from(range(NA_con))
	tmp_attr = {'graph': root, 'block': 'root'}
	tree.add_node(0)
	tree.nodes[0].update(tmp_attr.copy())
	del root, tmp_attr

	# generate all graphs
	list_of_graphs = create_list_graph(NA_con)

	# connected graph filter
	graphs_connected = []	
	for g in list_of_graphs:
		if nx.is_connected(g):
			graphs_connected.append(g.copy())
	del list_of_graphs, g   

	# isomorphism graph filter
	graphs_connected_iso = graphs_connected.copy()
	graphs_connected_iso = del_iso_graph(graphs_connected_iso,None)           
	del graphs_connected	

	#tree generator 1 child
	tmp_i = len(tree)
	for i,item in enumerate(graphs_connected_iso):
		i = i + tmp_i
		tmp_attr = {'graph': item, 'block': 'structure'}
		tree.add_node(i)
		tree.nodes[i].update(tmp_attr.copy())
		tree.add_edge(0,i)
	del tmp_i, tmp_attr, item

	#create all the possible permutation without repetition, set -> delete repetation
	if len(list_con_atom) == 1:
		tmp = (list_con_atom[0])
		list_color = list(tmp)
	else:    
		list_color = list(set(permutations(list_con_atom)))

	#coloring graph, in the allowed way by valence, and add info to each atoms
	list_graphs_connected_iso = graphs_connected_iso.copy()
	for k,struc in enumerate(list_graphs_connected_iso):
		tmp_col = []
		#adding attribute/color(atom) to each node for each possibile permutation of color and structure
		for j,col in enumerate(list_color):
			tmp_G = list_graphs_connected_iso[k].copy()
			for i,node in enumerate(tmp_G):
				tmp_attr = {'atom': col[i]}
				tmp_G.nodes[i].update(tmp_attr.copy())
				#adding all info to each node/atom
				for l in range(len(dict_element)):
					if col[i] == dict_element[l]['atom']:
						tmp_G.nodes[i].update(dict_element[l].copy())
			#delete graph in witch at least one atom doesn't respet the valence          
			if atom_valence(tmp_G)[0]:                   
				tmp_col.append(tmp_G.copy())
		#delete isomorphic colored-graphs generated by the k-structure 
		tmp = tmp_col.copy()
		tmp = del_iso_graph(tmp,'atom')
		#list_graph_connectted_iso is a list of list of colored-graph    
		list_graphs_connected_iso[k] = tmp
	del tmp_col, tmp_G, tmp_attr, tmp, graphs_connected_iso, col, node, struc, list_color, list_con_atom

	#tree generator 2 child
	tmp_k = len(tree) - len(list_graphs_connected_iso)
	tmp_num = len(tree)
	for j,item1 in enumerate(list_graphs_connected_iso):
		tmp_i = len(tree)
		tmp_j = j + tmp_k
		for i,item in enumerate(list_graphs_connected_iso[j]):
			i = i + tmp_i
			tmp_attr = {'graph': item, 'block': 'atom'}
			tree.add_node(i)
			tree.nodes[i].update(tmp_attr.copy())
			tree.add_edge(tmp_j,i)
	del tmp_k, tmp_attr, tmp_j, tmp_i, item1, item, list_graphs_connected_iso, 

	#Add Hs to the graph in all possibile way allowed
	if n_h > 0:
		tmp_nf = len(tree)
		for j,item in enumerate(tree):
			if j >= tmp_num:
				list_gr = H_addiction(tree.nodes[j]['graph'].copy(),dict_element)
				for k in range(n_h-1):
					tmp_a = []
					for i,tmp in enumerate(list_gr):
						tmp_a = tmp_a + H_addiction(tmp,dict_element)
					list_gr = tmp_a.copy()
				tmp_i = list_gr.copy()
				tmp_i = del_iso_graph(tmp_i,'atom')
				tmp_k = len(tree)
				tmp_g = nx.Graph()
				for i,tmp in enumerate(tmp_i):
					tmp_attr = {'graph': tmp, 'block': 'hydrogen'}
					tmp_g.add_node(tmp_k)
					tmp_g.nodes[tmp_k].update(tmp_attr.copy())
					tmp_g.add_edge(tmp_k,j)       
					tmp_k = tmp_k + 1
				tree = nx.compose(tree,tmp_g)
				del tmp_k, tmp_g, tmp_i, list_gr, tmp, i
		del item, j, tmp_num, tmp_attr
	else:
		tmp_nf = len(tree)
		tmp_k = len(tree)
		tmp_g = nx.Graph()
		for i,tmp in enumerate(tree):
			if i >= tmp_num:
				tmp_attr = {'graph': tree.nodes[i]['graph'].copy(), 'block': 'hydrogen'}
				tmp_g.add_node(tmp_k)
				tmp_g.nodes[tmp_k].update(tmp_attr)
				tmp_g.nodes[tmp_k].update()
				tmp_g.add_edge(tmp_k,i)       
				tmp_k = tmp_k + 1
		tree = nx.compose(tree,tmp_g)

	#generate all possibile bonds in molecules
	tmp_k = len(tree)
	tmp_nnf = len(tree)
	for i,item in enumerate(tree):
		if i >= tmp_nf:
			g = tree.nodes[i]['graph'].copy()
			#inizialization multi-bond
			for (u, v) in g.edges:
				g[u][v]['bond'] = 1
			#inizialization graph properties    
			g.graph['tot_electrons'] = tot_el
			g.graph['electrons'] = tot_el
			g.graph['radicals'] = 0 
			g.graph['lone_pairs'] = 0
			g.graph['N_bond_single'] = 0
			g.graph['N_bond_double'] = 0
			g.graph['N_bond_triple'] = 0
			g.graph['abs_total_formal_charge'] = 0
			#generate all possibile bonds in molecules 
			list_mol = list_molecules_bonds(g)
			#update tree
			tmp_g = nx.Graph()
			for j,tmp in enumerate(list_mol):
				tmp_attr = {'graph': tmp, 'block': 'Bond'}
				tmp_g.add_node(tmp_k)
				tmp_g.nodes[tmp_k].update(tmp_attr.copy())
				tmp_g.add_edge(tmp_k,i)       
				tmp_k = tmp_k + 1
			tree = nx.compose(tree,tmp_g)

	#generate all possibile electronics structure in molecules with bonds
	tmp_k = len(tree)
	for i,item in enumerate(tree):
		if i >= tmp_nnf:
			g = tree.nodes[i]['graph'].copy()
			list_ris = list_resonance_structure(g)
			#update tree
			tmp_g = nx.Graph()
			for j,tmp in enumerate(list_ris):
				tmp_attr = {'graph': tmp, 'block': 'electron'}
				tmp_g.add_node(tmp_k)
				tmp_g.nodes[tmp_k].update(tmp_attr.copy())
				tmp_g.add_edge(tmp_k,i)       
				tmp_k = tmp_k + 1
			tree = nx.compose(tree,tmp_g)
	return(tree)