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






# import sys
# from ete3 import Tree
# import random

# def get_json(node):
#     # Read ETE tag for duplication or speciation events
#     if not hasattr(node, 'evoltype'):
#         dup = random.sample(['N','Y'], 1)[0]
#     elif node.evoltype == "S":
#         dup = "N"
#     elif node.evoltype == "D":
#         dup = "Y"

#     node.name = node.name.replace("'", '')
        
#     json = { "name": node.name, 
#              "display_label": node.name,
#              "duplication": dup,
#              "branch_length": str(node.dist),
#              "common_name": node.name,
#              "seq_length": 0,
#              "type": "node" if node.children else "leaf",
#              "uniprot_name": "Unknown",
#              }
#     if node.children:
#         json["children"] = []
#         for ch in node.children:
#             json["children"].append(get_json(ch))
#     return json


# if __name__ == '__main__':
#     t = Tree('chunks_rt_name.nwk')
#     # if len(sys.argv) > 1:
#     #     t = Tree(sys.argv[1])
        
#     # else:
#     #     # create a random example tree
#     #     t = Tree()

#     #     t.populate(100, random_branches=True)
#     string_json = str(get_json(t))               
#     # TreeWidget seems to fail with simple quotes
#     print(string_json.replace("'", '"'))