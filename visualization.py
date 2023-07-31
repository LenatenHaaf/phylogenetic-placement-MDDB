import copy, glob, os
from Bio import Phylo, SeqIO
def rename_clade(tree, leaf1, leaf2, name):
    '''
    Rename the nodes in the tree that do not have a name into ' ' or the name of the subtree
    ''' 
    common_ancestor = tree.common_ancestor(leaf1, leaf2)
    for clade in tree.find_clades():
        if clade == common_ancestor:
            clade.name = name
            for clade_subtree in clade.find_clades():
                if clade_subtree.name == None:
                    clade_subtree.name = ' '
            break

    return tree

path_name = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/chunks_names/'
name_paths = glob.glob(os.path.join(path_name, '*.fasta'))
path = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/trees/'
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

    rename_tree = rename_clade(reference_tree, naam[0], naam[-2], name[0])
    reference_tree = rename_tree
    i += 1
for clade in reference_tree.find_clades():
    if clade.name == None:
        clade.name = ' '
Phylo.write(reference_tree, 'visualization_rt_names.nwk', 'newick')
