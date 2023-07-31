'''
Assign colours to chunks in pruned rt
'''
from Bio import Phylo
import pandas as pd
import glob, os
from ete3 import Tree, CircleFace, faces, TextFace, ImgFace, ProfileFace
from ete3.treeview import TreeStyle, NodeStyle
from scipy import stats
path = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/trees/'
path_name = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/chunks_names/'
reference_tree = 'chunks_rt_name.nwk'
subtrees = pd.read_csv('Test3/distance_test3.csv')

def colour_branch_one_sample(path, path_name, sample, treefile):
    '''
    Colour the subtrees that contain queries. Different colour for the number queries there are in the subtrees.
    '''
    name_number = {}
    name_paths = glob.glob(os.path.join(path_name, '*.fasta'))
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
        # len_path = len(tree_path)
        number = tree_path[-7:-4]
        name_number[number] = name[0]
        i+=1
    count_chunks = []
    subtrees = sample['tree'].astype(str).to_list()
    for subtree in subtrees:
        padded_number = str(subtree).zfill(3)
        name_chunk = name_number.get(padded_number)
        count_chunks.append(name_chunk)
    count_name = {}
    for chunkname in count_chunks:
        number = count_chunks.count(chunkname)
        count_name[chunkname] = number

    tree_file = treefile
    tree = Tree(tree_file)
     
    for n in tree.iter_leaves():
        name = n.name
        no_in_ch = count_name.get(name)
        if no_in_ch == 1: # Change these numbers if the number of queries in a subtree is very different
            nstyle = NodeStyle()
            nstyle["bgcolor"] = "#DC267F"
            nstyle["size"] = 10
            n.set_style(nstyle)
        if no_in_ch == 2: # Change these numbers if the number of queries in a subtree is very different
            nstyle = NodeStyle()
            nstyle["bgcolor"] = "#785EF0"
            nstyle["size"] = 15
            n.set_style(nstyle)
        if no_in_ch >= 3: # Change these numbers if the number of queries in a subtree is very different
            nstyle = NodeStyle()
            nstyle["bgcolor"] = "#FFB000"
            nstyle["size"] = 20
            n.set_style(nstyle)
    ts = TreeStyle()
    ts.show_leaf_name = True 
    ts.scale = 2000
    tree.render("tree_test3.png", tree_style=ts)
# colour_branch_one_sample(path=path, path_name=path_name, sample=subtrees, treefile=reference_tree)



def colour_branch_three_samples(path, path_name, sample1, sample2, sample3, treefile):
    '''
    Add dots to each leaf when there is a query in the subtree. Compare three different samples. The size of the dot corresponds to the number fo queries in the subtrees.
    '''
    subtrees_1 = sample1['tree'].astype(str).to_list()
    subtrees_2 = sample2['tree'].astype(str).to_list()
    subtrees_3 = sample3['tree'].astype(str).to_list()
    name_number = {}
    name_paths = glob.glob(os.path.join(path_name, '*.fasta'))
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
        # len_path = len(tree_path)
        number = tree_path[-7:-4]
        name_number[number] = name[0]
        i+=1
    full_subtrees_1 = []
    for subtree in subtrees_1:
        padded_number = str(subtree).zfill(3)
        name_chunk = name_number.get(padded_number)
        full_subtrees_1.append(name_chunk)
    full_subtrees_2 = []
    for subtree in subtrees_2:
        padded_number = str(subtree).zfill(3)
        name_chunk = name_number.get(padded_number)
        full_subtrees_2.append(name_chunk)
    full_subtrees_3 = []
    for subtree in subtrees_3:
        padded_number = str(subtree).zfill(3)
        name_chunk = name_number.get(padded_number)
        full_subtrees_3.append(name_chunk)
    count_name1 = {}
    for chunkname in full_subtrees_1:
        number = full_subtrees_1.count(chunkname)
        count_name1[chunkname] = number
    count_name2 = {}
    for chunkname in full_subtrees_2:
        number = full_subtrees_2.count(chunkname)
        count_name2[chunkname] = number
    count_name3 = {}
    for chunkname in full_subtrees_3:
        number = full_subtrees_3.count(chunkname)
        count_name3[chunkname] = number
    tree_file = 'chunks_rt_name.nwk'
    tree = Tree(tree_file)
    ts = TreeStyle()
    ts.show_leaf_name = True 
    ts.scale = 2000 

    for n in tree.iter_leaves():
        if n.name in full_subtrees_1 or n.name in full_subtrees_2 or n.name in full_subtrees_3:
            nstyle = NodeStyle()
            nstyle['bgcolor'] = '#DADADA'
            nstyle["size"] = 10
            n.set_style(nstyle)
        if n.name in full_subtrees_1:
            if count_name1.get(n.name) == 1: # Change these numbers if the number of queries in a subtree is very different
                circle_face = CircleFace(radius=20, color='#DC267F', style='circle')
            if count_name1.get(n.name) == 2: 
                circle_face = CircleFace(radius=30, color='#DC267F', style='circle')
            if count_name1.get(n.name) >= 3:
                circle_face = CircleFace(radius=40, color='#DC267F', style='circle')
            circle_face.margin_left = 5
            circle_face.margin_bottom = 5
            n.add_face(circle_face, column=0, position="aligned")
        if n.name in full_subtrees_2:
            if count_name2.get(n.name) == 1: # Change these numbers if the number of queries in a subtree is very different
                circle_face = CircleFace(radius=20, color='#FFB000', style='circle')
            if count_name2.get(n.name) == 2:
                circle_face = CircleFace(radius=30, color='#FFB000', style='circle')
            if count_name2.get(n.name) >= 3:
                circle_face = CircleFace(radius=40, color='#FFB000', style='circle')
            circle_face.margin_left = 5
            circle_face.margin_bottom = 5
            n.add_face(circle_face, column=1, position="aligned")
        if n.name in full_subtrees_3:
            if count_name3.get(n.name) == 1: # Change these numbers if the number of queries in a subtree is very different
                circle_face = CircleFace(radius=20, color='#785EF0', style='circle')
            if count_name3.get(n.name) == 2:
                circle_face = CircleFace(radius=30, color='#785EF0', style='circle')
            if count_name3.get(n.name) >= 3:
                circle_face = CircleFace(radius=40, color='#785EF0', style='circle')
            circle_face.margin_left = 5
            circle_face.margin_bottom = 5
            n.add_face(circle_face, column=2, position="aligned")
 
    tree.show()
    tree.render("tree_test_double3.png", tree_style=ts, dpi=300)

tree_file = 'chunks_rt_name.nwk' # Path to the pruned reference tree file

# Paths to the files with results of different tests
sample3 = pd.read_csv('Test3/distance_test3.csv') 
sample2 = pd.read_csv('Test2/Resultaten/distance_test2.csv')
sample1 = pd.read_csv('Test1/Resultaten/distance_test1.csv')
# colour_branch_three_samples(path=path, path_name=path_name, sample1=sample1, sample2=sample2, sample3=sample3, treefile=tree_file

