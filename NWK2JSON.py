from Bio import Phylo
import json

def tree_to_dict(node):
    data = {"name": node.name}
    if node.clades:
        data["children"] = [tree_to_dict(child) for child in node.clades]
    return data

tree = Phylo.read('visualization_rt_names.nwk', 'newick')
root = tree.clade
json_tree = tree_to_dict(root)

with open('visualization_rt_names_json.json', 'w') as f:
    json.dump(json_tree, f)



