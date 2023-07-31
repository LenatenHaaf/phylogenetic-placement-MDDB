import glob, os
import pandas as pd
# data = pd.read_csv('Test3/distance_test3_best_ete.csv')
path_name = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/chunks_names/' # Path to the folder with the subtree names
name_paths = glob.glob(os.path.join(path_name, '*.fasta'))
path = 'C:/Users/Lena/Documents/School/Jaar 4/Bachelorklas/Bachelor thesis/Data/trees/' # Path to the folder with tree numbers
tree_paths = glob.glob(os.path.join(path, '*.tre'))
i = 0 

# Create a dictionary with {tree number: tree name}
name_number = {}
for tree_path in tree_paths:
    name_path = name_paths[i]
    split_orde = name_path.split('o__')
    try:
        split_fam = split_orde[-1].split('_f__')
        name = split_fam[-1].split('.')
    except:
        name = split_orde[-1].split('.')
    number = tree_path[-7:-4]
    name_number[number] = name[0]
    i+=1
# print(name_number)



def colour_dict():   
    '''
    Create a dictionary with {colourname: (rgb values)} to colour the leaves in the visualization
    '''
    dict_colour = {}
    colour_path = 'colours.txt' # File with the latex colours
    colour_names = []
    with open(colour_path, 'r') as file:
        for line in file:
            row = line.split('{')
            name_colour = row[1].replace('}', '')
            colour_names.append(name_colour)
            rgb = row[-1].split('}')
            value = rgb[0].split(',')
            rgb_value = ''
            for i in value:
                val = round((float(i) * 255), 0)
                rgb_value += str(int(val)) + ','
            rgb_value = rgb_value[:-1]
            dict_colour[name_colour] = '(' + rgb_value + ')'
    # print(dict_colour)
    return colour_names, dict_colour
# colour_dict()

def annotation(data, begin, end, path):
    '''
    Generate the data for the annotation files for ITOL
    Needs a dataframe with: query, subtree (name) found by BLAST, is it the correct subtree?, posterior probability of the query, 
    difference in number of nodes, difference in distance
    '''
    colour_names = colour_dict()[0]
    dict_colour = colour_dict()[1]
    output_post = []
    output_nodes = []
    output_dist = []
    output_colour  = []
    leaf_name = []
    output_yes, output_maybe, output_no = [], [], []
    output_yes2, output_maybe2, output_no2, output_post2, output_nodes2, output_dist2 = [], [], [], [], [], []
    output_yes3, output_maybe3, output_no3, output_post3, output_nodes3, output_dist3 = [], [], [], [], [], []
    output_yes4, output_maybe4, output_no4, output_post4, output_nodes4, output_dist4 = [], [], [], [], [], []
    chunk_colour = {}
    test_data = data[begin:end]
    for index, row in test_data.iterrows():
        name_tree = row['subtree']
        leaf_name.append(name_tree)
        if name_tree not in chunk_colour.keys():
            kleur = colour_names[0]
            rgb_value = dict_colour.get(kleur)
            colour_names = colour_names[1:]
            chunk_colour[name_tree] = rgb_value
        else:
            rgb_value = chunk_colour.get(name_tree)
        correct = row['correct']
        post_prob = row['post prob']
        
        distance = row['distance']
        if leaf_name.count(name_tree) == 1:
            if correct == 'yes':
                output_yes.append(name_tree + ',1')
            if correct == 'maybe':
                output_maybe.append(name_tree + ',1')
            if correct == 'No':
                output_no.append(name_tree + ',1')
                continue
            nodes = int(row['number of nodes'])
            output_colour.append(name_tree + ' label_background ' + 'rgb' + rgb_value)
            output_post.append(name_tree + ' ' + str(post_prob))
            output_nodes.append(name_tree + ',' + str(nodes) + ',-1,#000000,normal,1,0')
            output_dist.append(name_tree + ' ' + str(distance))
        if leaf_name.count(name_tree) == 2:
            if correct == 'yes':
                output_yes2.append(name_tree + ',1')
            if correct == 'maybe':
                output_maybe2.append(name_tree + ',1')
            if correct == 'No':
                output_no2.append(name_tree + ',1')
                continue
            nodes = int(row['number of nodes'])
            output_post2.append(name_tree + ' ' + str(post_prob))
            output_nodes2.append(name_tree + ',' + str(nodes) + ',-1,#000000,normal,1,0')
            output_dist2.append(name_tree + ' ' + str(distance))
        if leaf_name.count(name_tree) == 3:
            if correct == 'yes':
                output_yes3.append(name_tree + ',1')
            if correct == 'maybe':
                output_maybe3.append(name_tree + ',1')
            if correct == 'No':
                output_no3.append(name_tree + ',1')
                continue
            nodes = int(row['number of nodes'])
            output_post3.append(name_tree + ' ' + str(post_prob))
            output_nodes3.append(name_tree + ',' + str(nodes) + ',-1,#000000,normal,1,0')
            output_dist3.append(name_tree + ' ' + str(distance))
        if leaf_name.count(name_tree) == 4:
            if correct == 'yes':
                output_yes4.append(name_tree + ',1')
            if correct == 'maybe':
                output_maybe4.append(name_tree + ',1')
            if correct == 'No':
                output_no4.append(name_tree + ',1')
                continue
            nodes = int(row['number of nodes'])
            output_post4.append(name_tree + ' ' + str(post_prob))
            output_nodes4.append(name_tree + ',' + str(nodes) + ',-1,#000000,normal,1,0')
            output_dist4.append(name_tree + ' ' + str(distance))
    output = [output_yes, output_maybe, output_no, output_dist, output_nodes, output_post, output_colour]
    output2 = [output_yes2, output_maybe2, output_no2, output_dist2, output_nodes2, output_post2]
    output3 = [output_yes3, output_maybe3, output_no3, output_dist3, output_nodes3, output_post3]
    output4 = [output_yes4, output_maybe4, output_no4, output_dist4, output_nodes4, output_post4]
    file_path_yes = path + "yes.txt"
    file_path_maybe = "maybe.txt"
    file_path_no = path + "no.txt"
    file_path_post = path + "post_prob1.txt"
    file_path_nodes = path + "nodes1.txt"
    file_path_dist = path + "dist1.txt"
    file_path_yes2 = path + "yes2.txt"
    file_path_maybe2 = path + "maybe2.txt"
    file_path_no2 = path + "no2.txt"
    file_path_post2 = path + "post_prob2.txt"
    file_path_nodes2 = path + "nodes2.txt"
    file_path_dist2 = path + "dist2.txt"
    file_path_yes3 = path + "yes3.txt"
    file_path_maybe3 = path + "maybe3.txt"
    file_path_no3 = path + "no3.txt"
    file_path_post3 = path + "post_prob3.txt"
    file_path_nodes3 = path + "nodes3.txt"
    file_path_dist3 = path + "dist3.txt"
    file_path_yes4 = path + "yes4.txt"
    file_path_maybe4 = path + "maybe4.txt"
    file_path_no4 = path + "no4.txt"
    file_path_post4 = path + "post_prob4.txt"
    file_path_nodes4 = path + "nodes4.txt"
    file_path_dist4 = path + "dist4.txt"
    file_path_colour = path + "colour1.txt"
    files = [file_path_yes, file_path_maybe, file_path_no, file_path_dist, file_path_nodes, file_path_post, file_path_colour]
    files2 = [file_path_yes2, file_path_maybe2, file_path_no2, file_path_dist2, file_path_nodes2, file_path_post2]
    files3 = [file_path_yes3, file_path_maybe3, file_path_no3, file_path_dist3, file_path_nodes3, file_path_post3]
    files4 = [file_path_yes4, file_path_maybe4, file_path_no4, file_path_dist4, file_path_nodes4, file_path_post4]
    j = 0
    for file_path in files:
        with open(file_path, 'w') as file:
            for line in output[j]:
                file.write(line + '\n')
        j+=1
    m = 0
    for file_path in files2:
        if len(output2[m]) == 0:
                m+=1
                continue
        else:
            with open(file_path, 'w') as file:
                for line in output2[m]:
                    file.write(line + '\n')
        m+=1
    n = 0
    for file_path in files3:
        if len(output3[n]) == 0:
                n+=1
                continue
        else:
            with open(file_path, 'w') as file:
                for line in output3[n]:
                    file.write(line + '\n')
        n+=1
    p = 0
    for file_path in files4:
        if len(output4[p]) == 0:
                p+=1
                continue
        else:
            with open(file_path, 'w') as file:
                for line in output4[p]:
                    file.write(line + '\n')
        p+=1

data = pd.read_csv('visualization.csv') # Path to the file with all results genreated using the definition visualization_file in Search_subtree.py
begin_test1 = 0
end_test1 = 25 # Data of test 1
begin_test2 = 25
end_test2 = 50 # Data of test two
begin_test3 = 50 
end_test3 = 81 # Data of test 3
path = 'Test3/Annotation/' # Path to the destination folder of all annotation files
# annotation(data=data, begin=begin_test3, end=end_test3, path=path)