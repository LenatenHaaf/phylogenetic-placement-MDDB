# phylogenetic-placement-MDDB
BSc project at Leiden University for performing phylogenetic placement on the backbone in MDDB using BLAST and pplacer. By Lena ten Haaf. The code in this repository is made for three purposes:
  1. Find the correct subtree based on the results of BLAST
  2. Additional code for running and analyzing the project
  3. Visualization of the results

Supervisors:
+ Prof. Dr. Ir. F. Verbeek (LIACS)
+ Dr. R. Vos (Naturalis)

## Find the correct subtree based on the results of BLAST
The definition search_subtree(path, test_path) in the [Search_subtree.py](Search_subtree.py) file is the definition that is used to find the subtree using the results of BLAST. 

## Additional code for running and analyzing the project
All other definitions in [Search_subtree.py](Search_subtree.py) are used to prepare files for the process or analyze the results. The definition distance_ete(maximum, s_sub, part_of_path) is used to calculate the number of nodes and the distance between two nodes in the tree. [run_taxid.py](run_taxid.py), [run_pplacer.sh](run_pplacer.sh), [run_raxml.sh](run_raxml.sh) and [guppy.sh](guppy.sh) are scripts to automate the running of the taxid create, pplacer, RAxML and guppy. 

## Visualization of the results
The following files are used to visualize the results using the ETE Toolkit or ITOL, or prepare the files to be able to use in D3.js for the interactive visualization that can be found [here](https://observablehq.com/d/1c69d26ecff13759):
+ [annotations.py](annotations.py) is used to generate the data for the annotation files for ITOL.
+ [prune_rt.py](prune_rt.py) is used to collapse all subtrees into one single node, to make the reference tree smaller
+ [visualization_ete.py](visualization_ete.py) is used to make a static visualization using the collapsed reference tree and ETE
+ [visualization.py](visualization.py) is used to change the names of the internal nodes, to be able to make the visualization in D3.js interactive
+ [NWK2JSON.py](NWK2JSON.py) is used to convert a .nwk file into a .json file to be able to upload it into D3.js
