import glob, os
from Bio import Phylo
def prune_tree(tree, leaf1, leaf2, name): 
    ''' 
    Prune a clade into one single node
    '''
    nodes_to_remove = []

    found = False
    for clade in tree.find_clades():
        if clade.name == leaf1 or clade.name == leaf2 and not found:
            found = True
            nodes_to_remove.append(clade)
            continue
        if found:    
            nodes_to_remove.append(clade)
        if clade.name == leaf1 or clade.name == leaf2 and found:
            break
    for node in nodes_to_remove:
        tree.collapse(node)
    for leaf in tree.get_terminals():
        if leaf.name == None:
            leaf.name = name

    return tree

path_name = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/chunks_names/' # Path to the folder with the files of the different subtrees with the original name
name_paths = glob.glob(os.path.join(path_name, '*.fasta'))
path = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/trees/' # Path to the tree files of the subtrees, with as name only the number of the tree
reference_tree = Phylo.read('l0.2_s3_4_1500_o1.0_a0_constr_localpair.tre', 'newick')
tree_paths = glob.glob(os.path.join(path, '*.tre'))
i = 0
for tree_path in tree_paths:
    name_path = name_paths[i]
    split_orde = name_path.split('o__')
    try:
        split_fam = split_orde[-1].split('_f__')
        name = split_fam[-1].split('.')
    except:
        name = split_orde[-1].split('.')
    naam = []
    name_tree = Phylo.read(tree_path, 'newick')
    for leaf in name_tree.get_terminals():
        naam.append(leaf.name)
    

    pruned_tree = prune_tree(reference_tree, naam[0], naam[-2], name[0])
    reference_tree = pruned_tree
    i += 1
Phylo.write(pruned_tree, 'chunks_rt_name.nwk', 'newick')
pruned_tree = Phylo.read('chunks_rt_name.nwk', 'newick')
Phylo.draw(pruned_tree)



