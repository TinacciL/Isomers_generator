# Isomers Generator

Isomers Generator is a python program for generate all the possible molecules starting from the molecular formula.
The porpouse of this open source program is not to provide a optimatize tool but a frindly tool to work and generate all the possibile molecules.    

## Installation

This program used a python3 interface, to run this code you must install on your machine this list of packages:

* ```matplotlib```
* ```networkx```
* ```itertools```
* ```random```

The program is inside the ```functions.py``` file, inside this file there are some file that can you help to visualize the molecules and show the properties of those.

## Usage

**Main code:**

```python
import networkx as nx
from IG_lib import isomers_generator

n_H = 1 # the number (int) of hydrogen atoms
n_C = 1 # the number (int) of carbon atoms
n_N = 1 # the number (int) of nitrogen atoms
n_O = 1 # the number (int) of oxygen atoms

tree = isomers_generator(n_H,n_C,n_N,n_O) # tree is a tree in which each node are a molecule in the process of creation, the leaf are the all possibile molecules generated
```

**Other functions:**

Function that print the tree.

```python
from IG_lib import hierarchy_pos

pos = hierarchy_pos(tree,0)    
nx.draw(tree, pos=pos, with_labels=True) 
```

Function that print one molecule from the tree:

```python
from IG_lib import mol_graph_image

i = 4 # the i-node of the tree
g = tree.nodes[i]['graph']
mol_graph_image(g)
```

Function that print all the info of the atoms in the molecule:

```python
from IG_lib import atoms_property

atoms_property(g)
```

## Documentation

* [Isomers Generator](https://github.com/TinacciL/Isomers_generator/blob/master/IG_documentation.pdf)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
