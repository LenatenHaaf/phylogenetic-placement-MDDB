import pandas as pd
from Bio import SeqIO
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
import os, glob
from collections import Counter
import random
import json
from ete3 import Tree
import time
start_time = time.time()

'''
Output of the different tests converted into a dataframe
'''
# -------- Output test 1 ---------
output_test1 = pd.read_csv('Test1/test1_output.csv')
columns = ['Query', 'Subject', 'Pident', 'Evalue', 'Length', 'Qcovhsp']
output_test1.columns = columns
output = output_test1[output_test1["Query"].str.contains("#") == False]
# --------------------------------

# ------- Output test 2 ----------
# output_test2 = pd.read_csv('Test2/test2_output.csv')
# columns = ['Query', 'Subject', 'Pident', 'Evalue', 'Length', 'Qcovhsp']
# output_test2.columns = columns
# output = output_test2[output_test2["Query"].str.contains("#") == False]
# --------------------------------


# # ------- Output test 3 ----------
# output_test3 = pd.read_csv('Test3/test3_output.csv')
# columns = ['Query', 'Subject', 'Pident', 'Evalue', 'Length', 'Qcovhsp']
# output_test3.columns = columns
# output = output_test3[output_test3["Query"].str.contains("#") == False]
# # --------------------------------


def random_test_set():
    '''
    Get a random test set of size n from the data set
    '''
    data = list(SeqIO.parse('l0.2_s3_4_1500_o1.0_a0_constr_localpair.fasta', "fasta")) #Specifiy the path to the fasta file of the whole tree
    seqid = []
    rec_list = []
    for seq in data:
        seqid.append(seq.id)
        rec = SeqRecord(seq.seq, id=seq.id)
        rec_list.append(rec)

    n = 25
    random_n = random.sample(seqid, n)
    # ----- These variables are used to add duplicates to the last test -----
    # double = ['SH1145193.08FU_EU222979_reps', 'SH1189121.08FU_MG920147_reps', 'SH1152267.08FU_KR017144_reps', 'SH2614834.08FU_UDB0680020_reps', 'SH1191717.08FU_MF488989_refs', 'SH1190303.08FU_MH926037_reps']
    # random_n += double
    # random_n_set = set(random_n)
    # print(len(random_n_set))
    # print(random_n)
    # -----------------------------
    data = []
    test = []
    for rec in rec_list:
        if rec.id in random_n:
            test.append(rec)
        else:
            data.append(rec)
    # print(len(test))
    # print(len(data))
    SeqIO.write(test, "Test3/test_3.faa", "fasta-2line") #Specify the path of the file the test set should be written to
    SeqIO.write(data, "Test3/data_3.faa", "fasta-2line") #Specify the path of the file the data set should be written to

def search_subtree(path, test_path):
    ''' 
    Determine the subtree for each query. Input is the path to the folder containing the fasta files of the subtrees (path) and the path to the fasta file with 
    queries that are tested (test_path).
    Dataframe with the query, the original tree, the length of the original tree, the subtrees with hits of BLAST, the number of BLAST
    hits in these subtrees, the length of these subtrees, if the correct subtree is found
    Creates a csv file of the dataframe

    The subtree names should be the treenumber

    Return: list of dictionaries with for each query the determined subtree, dictionary with for each query the original subtree
    '''
    q_sub = {} #{query: {subtree: number of hits}}
    s_sub = {} #{query: original subtree}
    length = {}
    len_path = len(path)
    columns = ['query', 'original_tree', 'length_original', 'subtrees', 'hits_blast', 'length_trees', 'correct']
    subtrees = pd.DataFrame(columns=columns)
    n = 3 #number of decimals for no subtree

    seq_list = list(SeqIO.parse(test_path, "fasta")) 

    seqid = []
    for seq in seq_list:
        seqid.append(seq.id)
    
    fasta_paths = glob.glob(os.path.join(path, '*.fasta'))
    for seq in seqid: #Loop over all test queries
        sub_trees = []
        hits_st = []
        len_subtrees = []
        for fasta_path in fasta_paths: #loop over all subtrees
            subtree = fasta_path[len_path: len_path+n]
            id = []
            for seq_record in SeqIO.parse(fasta_path, 'fasta'):
                id.append(seq_record.id)
            length[fasta_path[len_path: len_path+n]] = len(id)
            length_st = len(id)
            count_q = 0
            subtree_dict = {}
            if seq in id:
                s_sub[seq] = fasta_path[len_path: len_path+n] #Add the subtree in which they originally are
                original_tree = fasta_path[len_path: len_path+n]
                length_o_st = len(id)
            subject = output[output['Query'] == seq]
            for i in subject['Subject']:
                if i in id:
                    count_q += 1
            
            if count_q != 0:
                sub_trees.append(subtree)
                hits_st.append(count_q)
                len_subtrees.append(length_st)
                subtree_dict[fasta_path[len_path: len_path+n]] = count_q
            if subtree_dict != {}:
                if seq in q_sub:
                    q_sub[seq].update(subtree_dict)
                else:
                    q_sub[seq] = subtree_dict

        max_hits = max(hits_st)
        location = hits_st.index(max_hits)
        most_hits = sub_trees[location]
        if original_tree == most_hits:
            correct_tree = 'Yes'
        else:
            correct_tree = 'No'
        
        new_line = pd.DataFrame([[seq, original_tree, length_o_st, sub_trees, hits_st, len_subtrees, correct_tree]], columns=columns)
        subtrees = pd.concat([subtrees, new_line], ignore_index=True)

    maximum = [] #Will be a list containing dicts: {query: (subtree, number of hits)}
    for key, val in q_sub.items():
        count = Counter(val)
        maximum.append({key: count.most_common(1)[0]})
    

    # Information about the correct and wrong queries
    # correct = []
    # wrong = []
    # for i in range(len(maximum)):
    #     keys = [key for key in maximum[i]]
    #     query = maximum[i][keys[0]]
    #     sub = s_sub[keys[0]]
    #     if query[0] == str(sub):
    #         correct.append(keys[0])
    #     else:
    #         wrong.append(keys[0])

    #Print information about the results
    # print('Original subtree:')
    # print(s_sub)
    # print('BLAST:')
    # print(q_sub)
    # print('length:')
    # print(length) 
    # print('Number correct:')      
    # print(len(correct))
    # print('wrong:')
    # print(wrong)
    # print(maximum)

    # Write the results to a csv file
    # if len(subtrees) == len(seqid):
    #         subtrees.to_csv('Test3/subtrees_timetest.csv', index=False) #Print dataframe to csv
    return maximum, s_sub
        
def prune_tree(maximum):
    '''
    Prune the query/queries from the tree
    Input: the list of dicts from search_subtree(path, test_path)
    Return a newick file with the query/queries pruned
    '''
    no_tree_count = []
    for j in range(len(maximum)):
        keys = [key for key in maximum[j]]
        no_tree = maximum[j].get(keys[0])[0]
        no_tree_count.append(no_tree)
        if no_tree_count.count(no_tree) == 1: #Prune different queries from the same tree
            tree = Phylo.read('trees/' + no_tree + '.tre', 'newick')
        else: 
            tree = Phylo.read('pruned/pruned_tree_' + no_tree + '.nwk', 'newick')
        # print(keys[0])
        try: 
            prune = tree.common_ancestor(keys)
        except:
            print(keys[0] + 'is not in tree' + no_tree)
            continue
        tree.prune(prune)
        Phylo.write(tree, 'pruned/pruned_tree_' + no_tree + '.nwk', 'newick')

def remove_from_fasta(path, maximum, removed_path, query_path):
    '''
    Remove the queries from the fasta files of the subtrees
    Input: path to the folder that contains the fasta files of the subtrees (path),
    the list of dicts from search_subtree(path, test_path) (maximum),
    the path to the file the remaining data should be written to (removed_path),
    the path to the file where the query/queries should be written to (query_path)
    Return: a file with the sequence(s) of the query/queries and a file with the other sequences from the subtree

    '''
    len_path = len(path)
    n = 3
    fasta_paths = glob.glob(os.path.join(path, '*.fasta'))
    # Rename the subtree to only the number, if this is not done already
    for fasta in fasta_paths:
        no = fasta[len_path: len_path+n]
        new_name = path + no + '.fasta'
        os.rename(fasta, new_name)

    no_tree_count = []
    for j in range(len(maximum)):
        keys = [key for key in maximum[j]]
        no_tree = maximum[j].get(keys[0])[0]
        no_tree_count.append(no_tree)
        if no_tree_count.count(no_tree) == 1:
            query_list = list(SeqIO.parse(path + no_tree + '.fasta', "fasta"))
        else:
            query_list = list(SeqIO.parse(removed_path + no_tree + '.fasta', "fasta"))
        rec_list = []
        query = []
        for seq in query_list:
            if seq.id != keys[0]:
                rec = SeqRecord(seq.seq, id=seq.id)
                rec_list.append(rec)
                SeqIO.write(rec_list, removed_path + no_tree + ".fasta", "fasta-2line")
            else:
                rec = SeqRecord(seq.seq, id=seq.id)
                query.append(rec)
                try:
                    with open(query_path + no_tree + ".fasta", 'a') as file:
                        SeqIO.write(query, file, "fasta-2line")
                except:
                    SeqIO.write(query, query_path + no_tree + ".fasta", "fasta-2line")
              
def number_of_nodes(leaf1, leaf2, tree):
    '''
    Calculate the distance and the number of nodes between two leaves (using ete3)
    Input: two leaves in a tree, the tree
    Return: the number of nodes and the distance
    '''
    distance = tree.get_distance(leaf1, leaf2)
    total = tree.get_distance(leaf1, leaf2, topology_only=True) #Number of nodes
    return total, distance

def add_likelihood(number, name, part_of_path, mode):
    '''
    Retrieve LWR, number of possible positions, posterior probability and likelihood from output pplacer
    Part of path is the path to the folder with the .jplace and .tog results of pplacer
    '''
    path = part_of_path + "query_" + number + ".jplace"
    with open(path, "r") as file:
        jplace_data = json.load(file)
    placements = jplace_data["placements"]
    for placement in placements:
        place = placement['p']
        lwr = place[0][2]
        likelihood = place[0][3]
        sh = placement['nm']
        name_query = sh[0][0]
        no_of_placements = len(place)
        if mode == 'bay':
            post_prob = place[0][-1]
        else:
            post_prob = 'No posterior probability'
        if name_query == name:
            return lwr, no_of_placements, post_prob, likelihood
        
    return 'Not find', 'Not find'

def distance_ete(maximum, s_sub, part_of_path):
    '''
    Create a csv file with the difference in distance and in nodes for each query 
    Input: the list of dicts from search_subtree(path, test_path) (maximum),
    The dictionary with for each query the original subtree from search_subtree(path, test_path),
    The path to the folder where the results of pplacer are in (part_of_path)
    Return: a dataframe with: the query, branch length of the query in pplacer and in the original, the subtree number, the neighbour, 
    branch length of the neighbour in pplacer and in the original, the distance in pplacer, the number of nodes in pplacer, 
    the distance in the original tree, the number of nodes in the original tree, difference in distance, LWR, likelihood,
    number of possible placements, (posterior probability), difference in nodes
    '''
    columns = ['query', 'branch_len_q', 'tree', 'neighbour', 'branch_len_n', 'pplacer', 'nodes_pplacer', 'original', 'nodes_original', 'difference', 'lwr', 'likelihood', 'no_of_placements', 'difference nodes']
    columns_bay = ['query', 'branch_len_q', 'tree', 'neighbour', 'branch_len_n', 'pplacer', 'nodes_pplacer', 'original', 'nodes_original', 'difference', 'lwr', 'likelihood', 'no_of_placements', 'post_prob', 'difference nodes']
   
    distance = pd.DataFrame(columns=columns)
    length = 0 
    for i in range(len(maximum)):
        keys = [key for key in maximum[i]]
        placer = maximum[i][keys[0]]
        no_orig = s_sub[keys[0]]
        branch_len_n = []
        branch_len_q = []
        try: 
            tree_file_pplacer = part_of_path + 'query_' + placer[0] + '.tog.tre'
            pplacer_tree = Tree(tree_file_pplacer)
            
        except:
            print('No pplacer tree found for ', placer[0]) #Not the correct subtree was found in search_subtree, so pplacer is not executed
            length += 1
            continue
        pplacer_tree.set_outgroup('OUTGROUP')
        tree_file_original = 'trees/' + no_orig + '.tre'
        original_tree = Tree(tree_file_original)
        target_leaf_name = keys[0]
        targets = []
        target_leaf = original_tree.search_nodes(name=target_leaf_name)[0]
        targets.append(target_leaf.name)
        sisters = target_leaf.get_sisters()
        # Find a neighbour for the target leaf
        if sisters[0].name == '':
            sister1 = sisters[0].children[0]
            if sister1.name == '':
                sister2 = sister1.children[0]
                if sister2.name == '':
                    sister3 = sister1.children[1]
                    if sister3.name == '':
                        sister4 = sister2.children[0]
                        if sister4.name == '':
                            sister5 = sister2.children[1]
                            if sister5.name == '':
                                targets.append(sister5.name)
                                branch_len_n.append(sister5.dist)
                            else:
                                targets.append(sister5.name)
                                branch_len_n.append(sister5.dist)
                        else:
                            targets.append(sister4.name)
                            branch_len_n.append(sister4.dist)
                    else:
                        targets.append(sister3.name)
                        branch_len_n.append(sister3.dist)
                else:
                    targets.append(sister2.name)
                    branch_len_n.append(sister2.dist)
            else:
                targets.append(sister1.name)
                branch_len_n.append(sister1.dist)
        else:
            targets.append(sisters[0].name)
            branch_len_n.append(sisters[0].dist)

        k = 0
        for tree in [pplacer_tree, original_tree]:
            nodes = number_of_nodes(targets[0], targets[1], tree)[0]
            dist = number_of_nodes(targets[0], targets[1], tree)[1]
            max_length = 0
            for node in tree.iter_descendants():
                if node.name == 'OUTGROUP':
                    continue
                if node.dist > max_length:
                    max_length = node.dist
            if k == 0:
                nodes_pplacer = nodes
                pplacer = dist
                pplacer_normalized = dist/max_length # Normalize the distance in the subtree by dividing it through the longest branch in the subtree
            else: 
                nodes_original = nodes
                original = dist
                original_normalized = dist/max_length # Normalize the distance in the subtree by dividing it through the longest branch in the subtree
            k += 1
        neighbour = targets[1]
        query_original = original_tree.search_nodes(name=target_leaf_name)[0]
        branch_len_q.append(query_original.dist)
        query_pplacer = pplacer_tree.search_nodes(name=target_leaf_name)[0]
        branch_len_q.append(query_pplacer.dist)
        if neighbour != '':
            neighbour_pplacer = pplacer_tree.search_nodes(name=neighbour)[0]
            branch_len_n.append(neighbour_pplacer.dist)
        difference = abs(pplacer_normalized - original_normalized)
        difference_nodes = abs(nodes_original - nodes_pplacer)
        if 'bay' in part_of_path:
            mode = 'bay'
            lwr = add_likelihood(placer[0], target_leaf_name, part_of_path, mode)
            new_line = pd.DataFrame([[target_leaf_name, branch_len_q, placer[0], neighbour, branch_len_n, pplacer, nodes_pplacer, original, nodes_original, difference, lwr[0], lwr[3], lwr[1], lwr[2], difference_nodes]], columns=columns_bay)
            distance = pd.concat([distance, new_line], ignore_index=True)
        else: 
            mode = 'ml'
            lwr = add_likelihood(placer[0], target_leaf_name, part_of_path, mode)
            new_line = pd.DataFrame([[target_leaf_name, branch_len_q, placer[0], neighbour, branch_len_n, pplacer_normalized, nodes_pplacer, original_normalized, nodes_original, difference, lwr[0], lwr[3], lwr[1], difference_nodes]], columns=columns)
            distance = pd.concat([distance, new_line], ignore_index=True)
            
        if len(distance) == len(maximum) - length:
            distance = distance.sort_values('difference nodes')
            distance.to_csv('Test2/distance_test3_bay.csv', index=False)
            return distance
        
def trees_likelihood(subtrees_path, likelihood_path):
    '''
    Create a dataframe with: the query, the number of subtrees from BLAST, the likelihood score, LWR, length of the original subtree 
    and the tree number
    '''
    i = 0
    columns = ['query', 'no of trees', 'likelihood', 'lwr', 'length', 'tree']
    likelihood_df = pd.DataFrame(columns=columns)
    for path in subtrees_path:
        data_subtrees = pd.read_csv(path)
        data_likelihood = pd.read_csv(likelihood_path[i])
        for index1, row1 in data_subtrees.iterrows():
            for index2, row2 in data_likelihood.iterrows():
                if row1['query'] == row2['query']:
                    query = row1['query']
                    likelihood = abs(row2['likelihood'])
                    lwr = row2['lwr']
                    subtrees = row1['subtrees']
                    try:
                        subtrees_list = subtrees.split(',')
                        no_of_trees = len(subtrees_list)
                    except:
                        no_of_trees = 1
                    length = row1['length_original']
                    tree = row2['tree']
                    newline = pd.DataFrame([[query, no_of_trees, likelihood, lwr, length, tree]], columns=columns)
                    likelihood_df = pd.concat([likelihood_df, newline], ignore_index=True)
        i += 1
    likelihood_df.to_csv('no_of_trees_likelihood.csv', index=False)

def visualization_file(subtrees_path, likelihood_path, name_number):
    '''
    Combine all results from all three tests into one file for visualization
    Generates a dataframe with: query, subtree (name) found by BLAST, is it the correct subtree?, posterior probability of the query, 
    difference in number of nodes, difference in distance
    '''
    i = 0
    columns = ['query', 'subtree', 'correct', 'post prob', 'number of nodes', 'distance']
    visualization_df = pd.DataFrame(columns=columns)
    for path in subtrees_path:
        data_subtrees = pd.read_csv(path)
        data_likelihood = pd.read_csv(likelihood_path[i])
        for index1, row1 in data_subtrees.iterrows():
            if row1['correct'] == 'No':
                correct = 'No'
                query = row1['query']
                tree = row1['original_tree']
                tree_number = str(tree).zfill(3)
                tree_name = name_number.get(tree_number)
                newline = pd.DataFrame([[query, tree_name, correct, None, None, None]], columns=columns)
                visualization_df = pd.concat([visualization_df, newline], ignore_index=True) 
                continue
            for index2, row2 in data_likelihood.iterrows():
                if row1['query'] == row2['query']:
                    query = row1['query']
                    posterior = row2['post_prob']
                    subtrees = row1['subtrees']
                    try:
                        subtrees_list = subtrees.split(',')
                        no_of_trees = len(subtrees_list)
                    except:
                        no_of_trees = 1
                    if row1['correct'] == 'Yes' and no_of_trees == 1:
                        correct = 'yes'
                    elif row1['correct'] == 'Yes' and no_of_trees > 1:
                        correct = 'maybe'
                    if row1['correct'] == 'No':
                        correct = 'No'
                    tree = row2['tree']
                    tree_number = str(tree).zfill(3)
                    tree_name = name_number.get(tree_number)
                    nodes_ppl = row2['nodes_pplacer']
                    nodes_orig = row2['nodes_original']
                    difference_nodes = abs(nodes_orig - nodes_ppl)
                    difference_distance = row2['difference']
                    newline = pd.DataFrame([[query, tree_name, correct, posterior, difference_nodes, difference_distance]], columns=columns)
                    visualization_df = pd.concat([visualization_df, newline], ignore_index=True)
        i += 1
    visualization_df.to_csv('visualization.csv', index=False)



        


path_o1 = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Subtrees_queries/' #Path to the folder with the fasta files of all subtrees (name should be only the number of the subtree)
path_o2 = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test2/Subtrees_queries/'
path_o3 = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test3/Subtrees_queries/'

test_path_3 = 'Test3/test_3.faa'
test_path_2 = 'Test2/test_2.faa'
test_path_1 = 'Test1/test.faa'
# search_subtree(path_o3, test_path_3)

# prune_tree(maximum=search_subtree(path)[0])

removed_path = 'Test1/removed_fasta1/removed_' 
query_path = 'Test1/queries1/query_'
# remove_from_fasta(path=path_o1, maximum=search_subtree(path_o1)[0], removed_path, query_path)

part_of_path = 'Test3/output_test3_bay_best/' # Path to the folder with the .jplace and .tog files of pplacer
# distance_ete(search_subtree(path_o3, test_path_3)[0], search_subtree(path_o3, test_path_3)[1], part_of_path)

# calculate_distances(search_subtree(path_o1)[0], search_subtree(path_o1)[1], part_of_path)

subtrees_path = ['C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test1/Resultaten/subtrees_test1.csv', 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test2/Resultaten/subtrees_test2.csv', 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test3/subtrees_test3.csv']
likelihood_path = ['C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test1/distance_test1_best_ete_norm_likelihood.csv', 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test2/distance_test2_best_ete_norm_likelihood.csv', 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/Test3/distance_test3_best_ete_norm_likelihood.csv']
# trees_likelihood(subtrees_path=subtrees_path, likelihood_path=likelihood_path)

# Dictionary of the number of the subtree and the corresponding name
name_number = {'001': 'Glomerales', '002': 'Diversisporales', '003': 'Gigasporales', '004': 'Archaeosporales', '005': 'Paraglomerales', '006': 'GS24', '007': 'Thelephorales', '008': 'Gomphales', '009': 'Hygrophoraceae', '010': 'Cortinariaceae', '011': 'Inocybaceae', '012': 'Amanitaceae', '013': 'Lycoperdaceae', '014': 'Agaricaceae', '015': 'Typhulaceae', '016': 'Clavariaceae', '017': 'Hydnangiaceae', '018': 'Tricholomataceae', '019': 'Marasmiaceae', '020': 'Mycenaceae', '021': 'Psathyrellaceae', '022': 'Strophariaceae', '023': 'Callistosporiaceae', '025': 'Omphalotaceae', '026': 'Cyphellaceae', '027': 'Entolomataceae', '028': 'Pluteaceae', '029': 'Lyophyllaceae', '030': 'Pleurotaceae', '031': 'Pterulaceae', '032': 'Bolbitiaceae', '033': 'Catathelasmataceae', '034': 'Stephanosporaceae', '035': 'Cystostereaceae', '036': 'Hymenogastraceae', '037': 'Schizophyllaceae', '038': 'Agaricales_fam_Incertae_sedis', '039': 'Crepidotaceae', '041': 'Physalacriaceae', '042': 'Nidulariaceae', '045': 'Radulomycetaceae', '046': 'Pseudoclitocybaceae', '048': 'Hymenochaetales', '049': 'Polyporales', '050': 'Boletales', '051': 'Russulales', '052': 'Corticiales', '054': 'Auriculariales', '055': 'Geastrales', '056': 'Trechisporales', '057': 'Phallales', '058': 'Gloeophyllales', '059': 'Sebacinales', '060': 'Hysterangiales', '061': 'Atheliales', '062': 'Amylocorticiales', '063': 'Agaricomycetes_ord_Incertae_sedis', '064': 'Tremellodendropsidales', '065': 'GS28', '067': 'Lepidostromatales', '070': 'Tremellales', '071': 'Cystofilobasidiales', '072': 'Trichosporonales', '073': 'Holtermanniales', '074': 'Filobasidiales', '075': 'Cystobasidiomycetes_ord_Incertae_sedis', '076': 'Erythrobasidiales', '077': 'Cyphobasidiales', '078': 'Cystobasidiales', '080': 'Sporidiobolales', '081': 'Microbotryomycetes_ord_Incertae_sedis', '082': 'Microbotryales', '083': 'Leucosporidiales', '084': 'Kriegeriales', '086': 'Agaricostilbales', '087': 'Septobasidiales', '088': 'Pucciniales', '089': 'Platygloeales', '090': 'Helicobasidiales', '091': 'Exobasidiales', '092': 'Entylomatales', '093': 'Tilletiales', '094': 'Georgefischeriales', '095': 'Microstromatales', '099': 'Tritirachiales', '100': 'Geminibasidiales', '101': 'Ustilaginales', '102': 'Urocystidales', '104': 'Atractiellales', '105': 'Spiculogloeales', '107': 'Dacrymycetales', '108': 'Dacrymycetes_ord_Incertae_sedis', '110': 'Malasseziales', '111': 'Wallemiales', '116': 'Dothideomycetes_ord_Incertae_sedis', '117': 'Capnodiales', '118': 'Pleosporales', '119': 'Acrospermales', '120': 'Tubeufiales', '121': 'Botryosphaeriales', '122': 'Dothideales', '123': 'Venturiales', '124': 'Myriangiales', '125': 'Strigulales', '127': 'Abrothallales', '128': 'Mytilinidales', '129': 'Patellariales', '130': 'Mytilinidiales', '131': 'Trypetheliales', '132': 'Hysteriales', '133': 'Jahnulales', '135': 'Stigmatodiscales', '137': 'Valsariales', '139': 'Minutisphaerales', '140': 'Helotiales', '141': 'Erysiphales', '142': 'Rhytismatales', '143': 'Thelebolales', '144': 'Triblidiales', '146': 'Phacidiales', '149': 'Aspergillaceae', '150': 'Trichocomaceae', '151': 'Thermoascaceae', '152': 'Elaphomycetaceae', '153': 'Onygenales', '154': 'Chaetothyriales', '155': 'Verrucariales', '156': 'Phaeomoniellales', '157': 'Mycocaliciales', '158': 'Sclerococcales', '159': 'Coryneliales', '160': 'Pyrenulales', '161': 'Glomerellales', '162': 'Sordariomycetes_ord_Incertae_sedis', '163': 'Microascales', '164': 'Diaporthales', '165': 'Coniochaetales', '166': 'Sordariales', '167': 'Hypocreales', '168': 'Xylariales', '169': 'Magnaporthales', '170': 'Chaetosphaeriales', '171': 'Melanosporales', '172': 'Phyllachorales', '173': 'Pleurotheciales', '174': 'Myrmecridiales', '175': 'Branch06', '176': 'Ophiostomatales', '177': 'Conioscyphales', '178': 'Hypoceales', '179': 'Boliniales', '181': 'Calosphaeriales', '182': 'Annulatascales', '183': 'Togniniales', '184': 'Xenospadicoidales', '185': 'Coronophorales', '186': 'Pararamichloridiales', '187': 'Trichosphaeriales', '188': 'Lulworthiales', '189': 'Phomatosporales', '190': 'Falcocladiales', '191': 'Savoryellales', '196': 'Ostropales', '197': 'Lecanorales', '198': 'Caliciales', '199': 'Rhizocarpales', '200': 'Peltigerales', '201': 'Umbilicariales', '202': 'Acarosporales', '203': 'Pertusariales', '204': 'Arctomiales', '206': 'Trapeliales', '207': 'Teloschistales', '208': 'Lecanoromycetes_ord_Incertae_sedis', '209': 'Leprocaulales', '210': 'Lecideales', '211': 'Baeomycetales', '212': 'Candelariales', '213': 'Sarrameanales', '214': 'GS36', '218': 'Orbiliales', '219': 'GS33', '221': 'Taphrinales', '222': 'Saccharomycetales', '223': 'GS34', '224': 'Symbiotaphrinales', '227': 'Coniocybales', '228': 'Geoglossales', '229': 'Laboulbeniales', '230': 'Pyxidiophorales', '231': 'Archaeorhizomycetales', '232': 'GS31', '233': 'Arthoniales', '234': 'Lichenostigmatales', '236': 'Sareales', '237': 'Lichinales', '241': 'GS05', '242': 'GS08', '243': 'GS07', '244': 'GS03', '245': 'GS11', '246': 'Branch01', '247': 'GS06', '249': 'GS10', '250': 'GS04', '252': 'Branch03', '254': 'Spizellomycetales', '255': 'Rhizophydiales', '256': 'Lobulomycetales', '259': 'Rhizophlyctidales', '260': 'Synchytriales', '262': 'Basidiobolales', '263': 'Endogonales', '264': 'GS21', '265': 'GS22', '267': 'Mucorales', '268': 'GS23', '269': 'Umbelopsidales', '270': 'Blastocladiales', '271': 'GS15', '272': 'Mortierellales', '273': 'Neocallimastigales', '274': 'GS16', '277': 'Olpidiales', '278': 'Monoblepharidales', '280': 'Sanchytriales', '281': 'Zoopagales', '283': 'Kickxellales', '284': 'Entorrhizales'}
# visualization_file(subtrees_path=subtrees_path, likelihood_path=likelihood_path, name_number=name_number)

end_time = time.time()
run_time = end_time - start_time 
# print(run_time)

